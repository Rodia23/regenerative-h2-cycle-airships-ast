# -*- coding: utf-8 -*-
# [EN] Clean-burn simulation of the base case with the closed-form preconditioning (Section 4, Appendix A).
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# t_apendice.py
# Simulación con parámetros de preacondicionamiento calculados analíticamente
# según las derivaciones del Apéndice A (apendice.tex).
#
# Reemplaza los valores hardcodeados de t2_gain_scheduled.py:
#   ANTES: dm_H2_proactivo = -10 kg/s  (fijo)
#          t_anticipo      =  2.0 s    (fijo)
#   AHORA: dm_pre = η · ṁ_max          [A.1: L_eff + flujo proporcional]
#          t_anticipo = τ · t_max       [A.2: fórmula cúbica]
#
# Parámetros operacionales que el usuario fija:
#   DM_MAX    — capacidad del actuador [kg/s]
#   ETA       — fracción de DM_MAX usada en el precondicionamiento (0 < η ≤ 1)
#   D_LIM     — hundimiento máximo permitido antes del drop [m]
#   TAU       — fracción de t_max que se usa efectivamente (0 < τ ≤ 1)
#   DIRIGIBLE — plataforma a simular
# =============================================================================

import os
from datetime import datetime
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import plot_style
import parametros as p_base
import inercia
import minimax as mm
import rutas

# Las figuras se guardan en resultados/figuras/ con nombre fijo (rutas.py).

# =============================================================================
# Selección de dirigible y parámetros operacionales
# =============================================================================

# Plataforma por defecto: parametros.py (el main de _perdrop recorre las 4)
DIRIGIBLE = {
    'nombre':        p_base.nombre_dirigible,
    'Area':          p_base.NAVE_BASE['Area'],
    'm_vacio':       p_base.NAVE_BASE['m_vacio'],
    'V_casco':       p_base.NAVE_BASE['V_casco'],
    'k':             p_base.NAVE_BASE['k'],
    'm_carga':       p_base.m_carga,
}

DM_MAX  = float(p_base.max_flujo_H2)  # [kg/s]  capacidad del actuador (esquema estándar, comparación)
D_LIM   = p_base.D_lim_oper      # [m]     hundimiento máximo operacional antes del drop (D_lim)
TAU     = 1.0     # [—]     fracción de t_max realmente usada (1.0 = usar todo t_max)
H_MAX_PREC = p_base.h_max_precision  # [m] umbral de precisión (|Δh|_fine, §2.4); la quema limpia no sube
# El pulso proactivo corre a la capacidad instalada COMPLETA (η = 1): en el óptimo
# la quema limpia es una sola a ṁ*, así que no hay razón para acotar el pulso.
ETA     = 1.0

# =============================================================================
# A.1 — Factor de Sustentación Efectiva
# =============================================================================

def calcular_L_eff(rho_air, rho_H2):
    """
    L_eff = ρ_air/ρ_H2 + 8
    Cada kg de H2 quemado genera L_eff·g N de fuerza neta hacia abajo:
      - ρ_air/ρ_H2: pérdida de flotación
      - 8:          ganancia de masa por O2 absorbido (1 kg H2 + 8 kg O2 → 9 kg H2O)
    """
    return rho_air / rho_H2 + p_base.ratio_O2_H2


# =============================================================================
# A.2-A.6 — la forma cerrada (quema limpia) vive SOLO en minimax.py (minimax.limpio)
# =============================================================================


# =============================================================================
# Motor de física (adaptado de t2_gain_scheduled.py)
# =============================================================================

class _Params:
    def __init__(self, overrides=None):
        for k, v in vars(p_base).items():
            if not k.startswith('_'):
                setattr(self, k, v)
        if overrides:
            for k, v in overrides.items():
                setattr(self, k, v)


def correr(n_entregas=1, toneladas_por_drop=1, graficar=True, verbose=True):
    p = _Params({'m_carga': n_entregas * toneladas_por_drop * 1000.0,
                 'Area':         DIRIGIBLE['Area'],
                 'm_vacio':      DIRIGIBLE['m_vacio'],
                 'nombre_dirigible': DIRIGIBLE['nombre']})

    # Siempre 30 s — intervalo de 1 s mutilaba t_anticipo (≈3 s) en la gráfica base.
    # El xlim del gráfico log-X se ajusta por separado para la vista visual.
    intervalo = p.intervalo_multidrop
    t_final   = 160.0 if n_entregas == 1 else (n_entregas + 3) * intervalo

    # ── Condiciones iniciales ─────────────────────────────────────────────────
    rho_h0     = p.rho_aire_std   # densidad del aire en la altitud de referencia h0
    rho_H2_ini = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_ini   = (p.m_vacio + p.m_carga) / (rho_h0 / rho_H2_ini - 1)   # gas justo para flotar
    M_ini      = p.m_vacio + p.m_carga + m_H2_ini
    V_casco    = DIRIGIBLE['V_casco']                       # volumen total del casco (dato; constante)
    k_nave     = DIRIGIBLE['k']                             # masa añadida de esta nave
    M_I_ini    = inercia.inercial_casco(V_casco, k_nave, rho_h0)   # inercia en equilibrio (inercia.py)

    # ── Parámetros analíticos A.1-A.6 (orden correcto) ───────────────────────
    L_eff     = calcular_L_eff(rho_h0, rho_H2_ini)
    m_drop_kg = toneladas_por_drop * 1000.0

    # ── Quema limpia (minimax.limpio; hallazgo C1, 2026-09-27): una sola quema a ṁ* que cruza la
    # liberación, hunde D_lim antes de soltar y deja la nave en reposo con el H2 de equilibrio ─
    ex         = mm.limpio(M_I_ini, m_drop_kg, D_LIM, L_eff)
    alpha_star = ex['alpha']                            # = 1/(1+s)
    dm_opt     = ex['dm_star']                          # = K_pre·α*^(3/2)

    # El pulso proactivo corre a la capacidad instalada completa (η = 1): la quema es una sola.
    dm_pre = dm_opt

    # ── Por-drop: t_anticipo recalculado por masa (pulso a dm_opt = ṁ*) ──────
    # t_anticip = (6·M·D_lim/(L_eff·g·ṁ*))^(1/3): el pulso a ṁ* llega exacto a D_lim
    # y logra compensación α* (verificado: comp = α* analíticamente).
    t_drops   = [(k + 1) * intervalo for k in range(n_entregas)]

    t_anticipos = []
    for k in range(n_entregas):
        m_carga_k = p.m_carga - k * m_drop_kg
        M_k       = M_I_ini          # inercia antes de cada drop: equilibrio (inercia.py)
        # t_max con dm_pre (= dm_opt): (6·M·D_max/(L_eff·g·dm_pre))^(1/3)
        t_ant_k   = TAU * (6 * M_k * D_LIM / (L_eff * p.g * dm_pre)) ** (1.0/3.0)
        t_ant_k   = min(t_ant_k, intervalo * 0.95)   # no solapar con drop anterior
        t_anticipos.append(t_ant_k)

    t_pre_starts = [t_drops[k] - t_anticipos[k] for k in range(n_entregas)]

    if verbose:
        print("\n" + "="*66)
        print(f"  {DIRIGIBLE['nombre']} — Parámetros analíticos (Apéndice A)")
        print("="*66)
        print(f"  L_eff                = {L_eff:.3f}  (ρ_air/ρ_H2 + 8)")
        print(f"  D_lim operacional    = {D_LIM*100:.1f} cm  |  h_max (precisión) = {H_MAX_PREC} m")
        print(f"  s = √(1 − m/M_I)     = {ex['s']:.5f}   hundimiento total (1+s)·D_lim = {ex['hundimiento_total']:.3f} m")
        print(f"  α*                   = {alpha_star:.4f}  = 1/(1+s)  (fracción compensada antes de soltar)")
        print(f"  ṁ*                   = {dm_opt:.3f} kg/s  = K_pre·α*^(3/2)")
        print(f"  pulso proactivo      = ṁ* (η = 1)")
        print(f"\n  t_anticipo por drop (recalculado por masa):")
        for k in range(n_entregas):
            M_k = M_I_ini
            print(f"    Drop {k+1}: M_inercial≈{M_k:.0f} kg  "
                  f"→ t_anticipo={t_anticipos[k]:.3f} s  "
                  f"(precond. inicia t={t_pre_starts[k]:.2f} s)")
        print(f"\n  → Comparación vs hardcoded:")
        print(f"     dm_H2_proactivo : -10.0 kg/s  →  -{dm_pre:.2f} kg/s")
        print(f"     t_anticipo      :   2.0 s     →    {t_anticipos[0]:.2f} s (drop 1)")
        print("="*66 + "\n")

    # ── Simulación ODE ────────────────────────────────────────────────────────
    fase2               = [False]
    fase2_at_drop       = [0]
    ultimo_drop_anun    = -1
    precond_notif       = [False] * n_entregas   # evitar prints repetidos por paso ODE

    def dinamica(t, y):
        nonlocal ultimo_drop_anun
        v, h, m_H2, T_H2, Q_con, Q_gen, m_agua = y
        m_H2   = max(m_H2, 0.1)
        m_agua = max(m_agua, 0.0)

        drops_hechos  = min(int(t // intervalo), n_entregas)
        m_carga_actual = p.m_carga - drops_hechos * (toneladas_por_drop * 1000.0)
        # Lastre de agua (O2 atmosferico absorbido) participa en la inercia y en el
        # equilibrio de flotabilidad — igual que t2_gain_scheduled.py, segun §2.3.
        M_eff        = p.m_vacio + m_carga_actual + m_agua
        M_total_nave = M_eff + m_H2

        if drops_hechos > ultimo_drop_anun:
            if verbose and drops_hechos > 0:   # ignorar el "drop 0" en t=0
                print(f"  [t={t:6.1f}s] DROP #{drops_hechos}  "
                      f"carga restante: {m_carga_actual/1000:.1f} t")
            ultimo_drop_anun = drops_hechos

        rho_aire = p.rho_aire_std * np.exp(-(h - p.h0) / p.H_escala)
        V_H2     = (m_H2 * p.R_H2 * T_H2) / p.P_atmos
        V_eq     = (M_eff * p.R_H2 * p.T1) / \
                   (p.P_atmos * (rho_h0 / rho_H2_ini - 1))
        M_I      = inercia.inercial(M_total_nave, rho_aire, V_casco, V_H2, k_nave)

        # Controlador PD con ganancia adaptativa: Kd con la inercia real (zeta de diseño)
        zeta    = p.zeta
        # (autor, 2026-10-09) ANTES: Kd_auto = 2 * zeta * np.sqrt(p.Kp * M_I / (rho_aire * p.g)) con zeta = 1.35, que
        # con la rigidez real (rho_aire - rho_H2) g daba un amortiguamiento de 1.302; ahora un solo zeta, el real.
        Kd_auto = 2 * zeta * np.sqrt(p.Kp * M_I / ((rho_aire - rho_H2_ini) * p.g))
        V_cmd   = V_eq + p.Kp * (p.h0 - h) - Kd_auto * v
        V_error = V_cmd - V_H2

        # Gain scheduling: Fase 1 (emergencia) → Fase 2 (trim)
        if fase2[0] and drops_hechos > fase2_at_drop[0]:
            fase2[0] = False
        if not fase2[0] and drops_hechos > 0:
            carga_rest  = (n_entregas - drops_hechos) * toneladas_por_drop * 1000.0
            m_H2_eq_now = (p.m_vacio + carga_rest + m_agua) / (rho_h0 / rho_H2_ini - 1)
            if m_H2 <= m_H2_eq_now:
                fase2[0]        = True
                fase2_at_drop[0] = drops_hechos

        K_mass_eff = p.K_mass_trim if fase2[0] else p.K_mass_control
        dm = K_mass_eff * V_error
        if drops_hechos == 0:
            dm = 0.0

        # Pulso proactivo — ventanas de tiempo puro (sin depender de drops_hechos)
        # Las ventanas [t_pre_starts[k], t_drops[k]) son no solapantes por diseño,
        # por lo que el loop activa a lo sumo un pulso por instante de tiempo.
        for k in range(n_entregas):
            if t_pre_starts[k] <= t < t_drops[k]:
                dm = -dm_pre
                if verbose and not precond_notif[k]:
                    print(f"  [t={t:6.2f}s] PRECOND #{k+1} activado  "
                          f"dm={-dm_pre:.1f} kg/s  t_anticipo={t_anticipos[k]:.3f}s  "
                          f"→ drop en t={t_drops[k]:.1f}s")
                    precond_notif[k] = True
                break

        # Saturar en dm_opt, NO en DM_MAX — coherencia con el Minimax [A.6]:
        # dm_opt es la capacidad mínima instalada que balancea proactivo y reactivo.
        # Usar DM_MAX=30 kg/s haría que el controlador reactivo ignorase la optimización
        # y saturase ahí, desconectando la simulación de la teoría del apéndice.
        dm = np.clip(dm, -dm_opt, dm_opt)

        F_emp = rho_aire * p.g * V_H2
        F_arr = 0.5 * p.Cd * DIRIGIBLE['Area'] * rho_aire * v**2 * np.sign(v)
        dv    = (F_emp - M_total_nave * p.g - F_arr) / M_I
        dQgen  = abs(dm) * p.Q_comb_H2 * p.rend_gen if dm < -0.05 else 0.0
        d_agua = -(1 + p.ratio_O2_H2) * min(dm, 0.0)   # 1 kg H2 + 8 kg O2 -> 9 kg H2O (lastre)

        return [dv, v, dm, (p.T1 - T_H2)/5.0, 0.0, dQgen, d_agua]

    y0  = [0.0, p.h0, m_H2_ini, p.T1, 0.0, 0.0, 0.0]
    # OJO (M-12 g de la auditoría 3, revisado 2026-09-29): esta dinámica cambia estado DENTRO de la función (fase2,
    # drops_hechos), así que el resultado depende de los pasos de prueba del integrador: con tolerancias finas el agua
    # final deja de ser 8.94 x el H2 quemado. Se deja RK45 por defecto (como estaba) y NO se usa para cifras del
    # paper: de minimax_base.json solo se leen dm_opt y t_anticipo, que salen de la forma cerrada. Las simulaciones del
    # paper son las de drop_flujo_optimo_perdrop.py (por tramos, con eventos, VER-02).
    sol = solve_ivp(dinamica, [0, t_final], y0, method='RK45', max_step=0.05)

    t_arr = sol.t
    v_arr, h_arr, m_H2_arr = sol.y[0], sol.y[1], sol.y[2]
    Q_gen_arr  = sol.y[5]
    m_agua_arr = sol.y[6]
    a_arr     = np.gradient(v_arr, t_arr)
    dm_arr    = np.gradient(m_H2_arr, t_arr)

    if graficar:
        if n_entregas == 1:
            _graficar_base(t_arr, h_arr, v_arr, a_arr, dm_arr, m_H2_arr,
                           Q_gen_arr/3.6e6, p, t_lib=intervalo,
                           dm_pre=dm_pre, t_anticipo=t_anticipos[0],
                           dm_opt=dm_opt, alpha_star=alpha_star)
        else:
            _graficar_multi(t_arr, h_arr, v_arr, a_arr, dm_arr, m_H2_arr,
                            Q_gen_arr/3.6e6, p, n_entregas=n_entregas,
                            intervalo=intervalo, tons_p=toneladas_por_drop,
                            dm_pre=dm_pre, t_anticipo=t_anticipos[0],
                            t_pre_starts=t_pre_starts)

    dh_peak     = float(np.max(np.abs(h_arr[t_arr >= intervalo] - p.h0)))
    m_H2_neto   = float(m_H2_arr[0] - m_H2_arr[-1])
    m_H2_bruto  = float(m_H2_arr[0] - np.min(m_H2_arr))
    m_agua_fin  = float(m_agua_arr[-1])
    if verbose:
        print(f"  H2 neto = {m_H2_neto:.1f} kg   H2 bruto = {m_H2_bruto:.1f} kg   "
              f"agua (lastre O2) = {m_agua_fin:.1f} kg")
    return {'Energia_kWh': Q_gen_arr[-1]/3.6e6,
            'M_total': M_ini, 'M_inercial': M_I_ini, 'L_eff': L_eff, 's': ex['s'], 'hundimiento_total_m': ex['hundimiento_total'],
            't_reactiva': ex['t_reactiva'],
            'dm_pre': dm_pre, 't_anticipo': t_anticipos[0],
            'alpha_star': alpha_star, 'dm_opt': dm_opt, 'dh_peak_m': dh_peak,
            'm_H2_neto': m_H2_neto, 'm_H2_bruto': m_H2_bruto, 'm_agua_fin': m_agua_fin}


# =============================================================================
# Gráfica
# =============================================================================

def _graficar_base(t, h, v, a, dm, m_h2, E, p, t_lib,
                   dm_pre, t_anticipo, dm_opt, alpha_star):

    FIGSIZE  = (10, 12)
    FS_TITLE = 14
    FS_LEG   = 10
    LW       = 2.2

    COL_ALT   = 'steelblue'
    COL_VEL   = 'forestgreen'
    COL_ACC   = 'indianred'
    COL_FLUJO = 'firebrick'
    COL_MASA  = 'darkorchid'
    COL_ENERG = 'darkorange'
    COL_LIB   = 'black'

    plt.rc('axes', labelsize=11)

    fig, axes = plt.subplots(3, 1, figsize=FIGSIZE, sharex=True,
                             gridspec_kw={'height_ratios': [2.5, 1.5, 1]})
    plt.subplots_adjust(hspace=0.15)

    # Panel 1: Cinemática
    ax_alt = axes[0]
    ax_vel = ax_alt.twinx()
    ax_acc = ax_alt.twinx()
    ax_acc.spines["right"].set_position(("axes", 1.12))
    ax_acc.spines["right"].set_visible(True)

    ax_alt.plot(t, h, color=COL_ALT, lw=LW,       label='Altitude [m]')
    ax_vel.plot(t, v, color=COL_VEL, lw=LW, ls='--', label='Velocity [m/s]')
    ax_acc.plot(t, a, color=COL_ACC, lw=1.8, ls='-.', alpha=0.7, label='Acceleration [m/s²]')

    ax_alt.set_xscale('linear')
    ax_alt.set_xlim(0, t[-1])   # simulación completa desde t=0
    ax_alt.axhline(p.h0, color=COL_ALT, ls=':', alpha=0.5)
    ax_alt.axvline(t_lib, color=COL_LIB, ls='--', lw=1.5, label=f'Release ({t_lib} s)')

    subtitle = (
        f'$\\dot{{m}}_{{\\mathrm{{pre}}}}={dm_pre:.1f}$ kg/s'
        f'  |  $t_{{\\mathrm{{anticip}}}}={t_anticipo:.2f}$ s'
        f'  |  $\\alpha^*={alpha_star:.3f}$'
        f'  |  $\\dot{{m}}_{{\\mathrm{{opt}}}}={dm_opt:.2f}$ kg/s'
        f'  [Appendix A]'
    )
    _ = subtitle  # (sin título dentro de la figura: va en el caption del manuscrito)
    ax_alt.set_ylabel('Altitude [m]',      color=COL_ALT, weight='bold')
    ax_vel.set_ylabel('Velocity [m/s]',    color=COL_VEL, weight='bold')
    ax_acc.set_ylabel('Acceleration [m/s²]',      color=COL_ACC, weight='bold')
    ax_alt.grid(True, which='both', alpha=0.2)

    lines  = (ax_alt.get_legend_handles_labels()[0] +
              ax_vel.get_legend_handles_labels()[0] +
              ax_acc.get_legend_handles_labels()[0])
    labels = (ax_alt.get_legend_handles_labels()[1] +
              ax_vel.get_legend_handles_labels()[1] +
              ax_acc.get_legend_handles_labels()[1])
    ax_alt.legend(lines, labels, loc='upper right', fontsize=FS_LEG)

    # Panel 2: Flujo + Masa
    ax_f = axes[1]
    ax_m = ax_f.twinx()
    ax_f.plot(t, dm,   color=COL_FLUJO, lw=LW,  label='H₂ Flow [kg/s]')
    ax_m.plot(t, m_h2, color=COL_MASA,  lw=1.8, ls=':', marker='.', markersize=2,
              markevery=15, alpha=0.5, label='H₂ Mass [kg]')
    ax_f.axvline(t_lib, color=COL_LIB, ls='--', lw=1.2)
    ax_f.set_ylabel('H₂ Flow [kg/s]',  color=COL_FLUJO, weight='bold')
    ax_m.set_ylabel('H₂ Mass [kg]',    color=COL_MASA,  weight='bold')
    ax_f.grid(True, which='both', alpha=0.2)
    ax_f.legend([ax_f.lines[0], ax_m.lines[0]],
                ['H₂ Flow [kg/s]', 'H₂ Mass [kg]'],
                loc='lower right', fontsize=FS_LEG)

    # Panel 3: Energía
    axes[2].plot(t, E, color=COL_ENERG, lw=LW)
    axes[2].axvline(t_lib, color=COL_LIB, ls='--', lw=1.2)
    axes[2].set_ylabel('Energy [kWh]',  weight='bold')
    axes[2].set_xlabel('Time [s]', weight='bold')
    axes[2].grid(True, which='both', alpha=0.2)

    plt.tight_layout()
    nombre = DIRIGIBLE['nombre'].replace(' ', '')
    path   = rutas.figura(f"minimax_base_{nombre}")
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f"  Figura guardada: {path}")
    rutas.mostrar()


def _graficar_multi(t, h, v, a, dm, m_h2, E, p, n_entregas, intervalo,
                    tons_p, dm_pre, t_anticipo, t_pre_starts=None):
    """5-panel linear-scale figure for multi-drop missions."""

    FIGSIZE   = (12, 16)
    LW        = 2.2
    COL_ALT   = 'steelblue'
    COL_VEL   = 'forestgreen'
    COL_ACC   = 'indianred'
    COL_FLUJO = 'firebrick'
    COL_MASA  = 'darkorchid'
    COL_ENERG = 'darkorange'
    COL_LIB   = 'black'

    plt.rc('axes', labelsize=12)
    plt.rc('axes', titlesize=14)

    fig, axes = plt.subplots(5, 1, figsize=FIGSIZE, sharex=True)

    # Líneas verticales en cada drop
    t_drops = [intervalo * (k + 1) for k in range(n_entregas)]

    # Panel 1: Altitud
    axes[0].plot(t, h, COL_ALT, lw=LW)
    axes[0].axhline(p.h0, color=COL_ALT, ls=':', alpha=0.5)
    axes[0].set_ylabel('Altitude [m]', weight='bold')
    # (sin título dentro de la figura)

    # Panel 2: Velocidad + Aceleración
    ax_v = axes[1]
    ax_a = ax_v.twinx()
    l_v, = ax_v.plot(t, v, COL_VEL, lw=LW)
    l_a, = ax_a.plot(t, a, COL_ACC, lw=2.0, ls='--')
    ax_v.set_ylabel('Velocity [m/s]',    color=COL_VEL, weight='bold')
    ax_a.set_ylabel('Acceleration [m/s²]',      color=COL_ACC, weight='bold')
    vl = max(abs(ax_v.get_ylim()[0]), abs(ax_v.get_ylim()[1]))
    al = max(abs(ax_a.get_ylim()[0]), abs(ax_a.get_ylim()[1]))
    ax_v.set_ylim(-vl*1.1, vl*1.1)
    ax_a.set_ylim(-al*2.2, al*2.2)
    ax_v.legend([l_v, l_a], ['Velocity [m/s]', 'Acceleration [m/s²]'],
                loc='upper right', fontsize=9)

    # Panel 3: Flujo H2 — con marcadores de inicio de precondicionamiento
    axes[2].plot(t, dm, COL_FLUJO, lw=LW)
    axes[2].set_ylabel('H₂ Flow [kg/s]', weight='bold')
    if t_pre_starts:
        for k, tps in enumerate(t_pre_starts):
            axes[2].axvline(tps, color='royalblue', ls=':', lw=1.4, alpha=0.8)
            axes[2].annotate(f'P{k+1}\n{tps:.1f}s', xy=(tps, 0),
                             xytext=(tps + 0.3, -dm_pre * 0.35),
                             fontsize=7, color='royalblue', alpha=0.9)

    # Panel 4: Masa H2
    axes[3].plot(t, m_h2, COL_MASA, lw=LW)
    axes[3].set_ylabel('H₂ Mass [kg]', weight='bold')

    # Panel 5: Energía
    axes[4].plot(t, E, COL_ENERG, lw=LW)
    axes[4].set_ylabel('Energy [kWh]',   weight='bold')
    axes[4].set_xlabel('Time [s]',       weight='bold')

    # Líneas de drop en todos los paneles
    for ax in axes:
        for td in t_drops:
            ax.axvline(td, color=COL_LIB, ls='--', lw=1.2, alpha=0.7)
        ax.grid(True, alpha=0.2)

    fig.tight_layout()
    nombre = DIRIGIBLE['nombre'].replace(' ', '')
    path   = rutas.figura(f"minimax_multi{n_entregas}_{nombre}_dmfijo")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Figura guardada: {path}")
    rutas.mostrar()


# =============================================================================
# Entrada principal
# =============================================================================

if __name__ == '__main__':
    # Misión base (1 drop)
    resultado = correr(n_entregas=1, toneladas_por_drop=1, graficar=True)
    print(f"\n  Energía generada: {resultado['Energia_kWh']:.1f} kWh"
          f"  |  pico |Δh| = {resultado['dh_peak_m']:.2f} m")
    rutas.guardar_datos('minimax_base', {'plataforma': DIRIGIBLE['nombre'], 'D_lim_m': D_LIM,
                                         'h_max_m': H_MAX_PREC, **resultado})

    # Descomenta para probar multi-drop:
    # correr(n_entregas=4, toneladas_por_drop=1, graficar=True)
