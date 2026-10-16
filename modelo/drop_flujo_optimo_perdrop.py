# -*- coding: utf-8 -*-
# [EN] Vertical-dynamics simulator of a release: burn, controller, heat-rejection unit and hold-down device (Sections 3-5).
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# drop_flujo_optimo_perdrop.py
# Igual que drop_flujo_optimo.py salvo que dm_opt (flujo óptimo Minimax A.6)
# se recalcula para cada drop usando la masa real M_k en lugar de M_ini fija.
#
# Cambios respecto al original:
#   - Bloque A.3-A.6 + loop t_anticipos fusionados → un solo loop por drop
#   - dm_pres[k], dm_opts[k], alpha_stars[k]: listas indexadas por drop
#   - dinamica: pulso proactivo usa dm_pres[k]; saturación usa dm_opts[0]
#     (actuador dimensionado para el peor caso = drop 1, masa máxima)
#   - Todo lo demás: idéntico al original
# =============================================================================

import os
from datetime import datetime
import numpy as np
from types import SimpleNamespace
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import plot_style
import parametros as p_base
import condensador
import minimax as mm
import inercia
import rutas

# Las figuras se guardan en resultados/figuras/ con nombre fijo (rutas.py).

# =============================================================================
# Selección de dirigible y parámetros operacionales
# =============================================================================

# Plataforma por defecto: parametros.py (el main de _perdrop recorre las 4). usar_nave() la cambia.
DIRIGIBLE = {}


def usar_nave(nombre, m_carga=None):
    """Deja la nave `nombre` (clave de parametros.PLATAFORMAS) en DIRIGIBLE, con m_carga [kg] a bordo."""
    d = p_base.PLATAFORMAS[nombre]
    DIRIGIBLE.update({'nombre': nombre, 'Area': d['Area'], 'm_vacio': d['m_vacio'], 'V_casco': d['V_casco'],
                      'k': d['k'], 'm_carga': p_base.m_carga if m_carga is None else m_carga})


usar_nave(p_base.nombre_dirigible)

DM_MAX  = float(p_base.max_flujo_H2)  # [kg/s]  capacidad del actuador (esquema estándar, comparación)
D_LIM   = p_base.D_lim_oper      # [m]     hundimiento máximo operacional antes del drop (D_lim)
TAU     = 1.0     # [—]     fracción de t_max realmente usada (1.0 = usar todo t_max)
# Integración (VER-02, verificación numérica): tolerancia relativa y escala de cada estado para la
# absoluta: v [m/s], h [m], m_H2 [kg], T_H2 [K], Q_con [J], Q_gen [J], m_agua [kg] (+ P_c [W],
# agua venteada [kg], vapor retenido [kg] con condensador)
RTOL = 1e-6
METODO = 'LSODA'   # integrador (rígido: la ganancia de emergencia hace un lazo de ~1 ms; ver correr())
ESCALAS_ESTADO = (1.0, 1.0, 1e3, 1e2, 1e9, 1e9, 1e2)
ESCALAS_COND = (1e8, 1e2, 1e2)
H_MAX_PREC = p_base.h_max_precision  # [m] umbral de precisión (|Δh|_fine, §2.4); la quema limpia no sube
ETA     = 1.0     # el pulso proactivo corre a la capacidad instalada completa (η = 1)

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


def correr(n_entregas=1, toneladas_por_drop=1, graficar=True, verbose=True,
           h_max=None, D_lim=None, dm_opt=None, intervalo=None, t_post=130.0,
           inyeccion=True, cond=None, devolver_serie=False, max_step=0.05, rtol=RTOL, atol=None,
           metodo=None, precond=True, retener=True, helices=None):
    """Simula la misión con la quema limpia (minimax.limpio).

    Por defecto (sin argumentos opcionales) reproduce exactamente las corridas del
    artículo: D_lim = 0.10 m, ṁ* de la quema limpia recalculado por entrega; h_max
    (1 m) solo es el umbral de las métricas de asentamiento.
    Para el mapa de requisitos (1 entrega): h_max, D_lim y dm_opt se pasan
    explícitos; intervalo = instante de la liberación, que debe dejar espacio a la
    ventana de anticipación; t_post = s simulados después. inyeccion=False: el
    controlador solo quema (no devuelve H2 de los tanques). precond=False: sin quema previa, el
    controlador responde solo después de soltar (control REACTIVO puro, con el agua retenida al instante;
    reactivo_vs_limpia.py, revisión 2026-10-01); dm_opt es entonces la saturación del controlador.
    Opciones del dispositivo de sujeción (freno.py, sección 5.2 del paper):
      retener=False: el escape se ventea (sin agua a bordo): L = rho_aire/rho_H2 - 1 en la quema previa;
      helices = {'F_max': N, 'tau': s, 'wn': rad/s, 'zeta': -, 'A': m2, 'FM': -, 'bidireccional': bool}: las hélices
      vectoriales de la nave, activas desde la primera liberación. Fuerza vertical F (positiva hacia abajo) con
      retardo de primer orden tau; la orden cancela el desbalance estático (flotación - peso) y corrige la altitud con
      un lazo de segundo orden (wn, zeta), saturada en F_max (y en -F_max si son bidireccionales).
      helices['anticipar'] = {'rampa': N/s, 't_previa': s} (revisión 2026-10-06, freno.py): el freno ANTICIPA como la
      quema limpia. Su orden sube en rampa desde t_drop - t_previa - tau (adelantada un tau para el retardo) y, tras
      soltar, la parte estática de la orden no pasa de la rampa. Sin esta clave, todo queda como antes.
    Con los valores por defecto el resultado es el del artículo.
    """
    h_max = H_MAX_PREC if h_max is None else h_max
    D_lim = D_LIM if D_lim is None else D_lim
    metodo = METODO if metodo is None else metodo
    p = _Params({'m_carga': n_entregas * toneladas_por_drop * 1000.0,
                 'Area':         DIRIGIBLE['Area'],
                 'm_vacio':      DIRIGIBLE['m_vacio'],
                 'nombre_dirigible': DIRIGIBLE['nombre']})

    intervalo = p.intervalo_multidrop if intervalo is None else intervalo
    # Ventana: último drop en n_entregas*intervalo; con el lastre de agua (§2.3) la
    # nave se asienta ~30-40 s después → +2 intervalos alcanzan para mostrar el
    # regreso completo a (h0, 0) con una cola corta. (1 entrega: 30 + 130 = 160 s)
    t_final   = intervalo + t_post if n_entregas == 1 else (n_entregas + 2) * intervalo

    # ── Condiciones iniciales ─────────────────────────────────────────────────
    rho_h0     = p.rho_aire_std   # densidad del aire en la altitud de referencia h0
    rho_H2_ini = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_ini   = (p.m_vacio + p.m_carga) / (rho_h0 / rho_H2_ini - 1)   # gas justo para flotar
    M_ini      = p.m_vacio + p.m_carga + m_H2_ini
    V_casco    = DIRIGIBLE['V_casco']                       # volumen total del casco (dato; constante)
    k_nave     = DIRIGIBLE['k']                             # masa añadida de esta nave (Lamb, su esbeltez)
    if m_H2_ini * p.R_H2 * p.T1 / p.P_atmos > V_casco * (1.0 + 1e-9):
        raise ValueError(f"{DIRIGIBLE['nombre']}: {p.m_carga:.0f} kg no flotan en {V_casco:.0f} m3")
    I_eq       = inercia.inercial_casco(V_casco, k_nave, rho_h0)   # inercia en equilibrio (inercia.py)

    L_eff     = calcular_L_eff(rho_h0, rho_H2_ini) if retener else rho_h0 / rho_H2_ini - 1.0
    if cond is not None and not retener:
        raise ValueError('retener=False es sin unidad de rechazo: no combinar con cond')
    i_hel = 7 + (3 if cond is not None else 0)                 # índice del estado de la fuerza de las hélices
    m_drop_kg = toneladas_por_drop * 1000.0

    t_drops = [(k + 1) * intervalo for k in range(n_entregas)]

    # ── Quema limpia (minimax.limpio) — recalculada por drop con M_k real ────
    # α* = 1/(1+s) ≈ 1/2 en todos los drops; ṁ*_k cambia con M_k.

    dm_opts     = []
    alpha_stars = []
    dm_pres     = []
    t_anticipos = []
    t_ant_teoricos = []

    for k in range(n_entregas):
        m_carga_k    = p.m_carga - k * m_drop_kg
        M_k          = I_eq          # inercia antes de cada drop: equilibrio, (1 + k)·rho·V_casco

        # quema limpia (minimax.limpio; después de soltar acelera M_k - m_drop; hallazgo C1, 2026-09-27)
        ex_k         = mm.limpio(M_k, m_drop_kg, D_lim, L_eff)
        alpha_star   = ex_k['alpha']                            # = 1/(1+s)
        dm_opt_k     = ex_k['dm_star']                          # = K_pre·α*^(3/2)
        if dm_opt is not None and k == 0:
            dm_opt_k = dm_opt                                   # ṁ impuesto (mapa, verificación)

        dm_opts.append(dm_opt_k)
        alpha_stars.append(alpha_star)
        dm_pres.append(dm_opt_k)   # pulso a ṁ* (η = 1)

        t_ant_k = TAU * (6 * M_k * D_lim / (L_eff * p.g * dm_opt_k)) ** (1/3)
        t_ant_teoricos.append(t_ant_k)
        t_ant_k = min(t_ant_k, intervalo * 0.95)
        t_anticipos.append(t_ant_k)

    if not precond:                                    # reactivo puro: sin ventana de anticipación
        t_anticipos = [0.0] * n_entregas
    t_pre_starts = [t_drops[k] - t_anticipos[k] for k in range(n_entregas)]
    if cond is not None and cond.get('prearranque'):
        # operación del paper (sec. 2.3): la unidad de rechazo recibe la orden n_tau_arranque·tau antes
        # de la primera quema previa (condensador.config)
        cond = dict(cond, t_arranque=t_pre_starts[0] - p.n_tau_arranque * cond['tau'])

    if verbose:
        print("\n" + "="*74)
        print(f"  {DIRIGIBLE['nombre']} — Parámetros analíticos por drop (Apéndice A)")
        print("="*74)
        print(f"  L_eff = {L_eff:.3f}  |  D_lim = {D_lim*100:.1f} cm  |  h_max = {h_max} m  |  "
              f"α* = 1/(1+s) = {alpha_star:.4f}")
        print(f"\n  {'Drop':>4}  {'M_k [kg]':>10}  {'α*':>7}  "
              f"{'ṁ_opt [kg/s]':>13}  {'t_ant [s]':>10}  {'t_pre [s]':>10}")
        print(f"  {'-'*4}  {'-'*10}  {'-'*7}  {'-'*13}  {'-'*10}  {'-'*10}")
        for k in range(n_entregas):
            m_carga_k = p.m_carga - k * m_drop_kg
            M_k       = I_eq
            print(f"  {k+1:>4}  {M_k:>10.0f}  {alpha_stars[k]:>7.4f}  "
                  f"{dm_opts[k]:>13.3f}  {t_anticipos[k]:>10.3f}  {t_pre_starts[k]:>10.2f}")
        print(f"\n  → Comparación vs hardcoded (drop 1):")
        print(f"     dm_H2_proactivo : -10.0 kg/s  →  -{dm_pres[0]:.2f} kg/s")
        print(f"     t_anticipo      :   2.0 s     →    {t_anticipos[0]:.2f} s (drop 1)")
        print("="*74 + "\n")

    # ── Simulación ODE, por tramos (M13 de la auditoría, 2026-09-27) ──────────
    # El modo del controlador (entregas hechas, pulso previo activo, fase de trim) es CONSTANTE dentro de
    # cada tramo. Los cambios de instante conocido (inicio de cada pulso previo, cada liberación, arranque
    # de la unidad de rechazo) son bordes de tramo, y el paso a trim es un evento de solve_ivp. Antes el
    # modo se cambiaba DENTRO del lado derecho, que el integrador evalúa también en etapas y pasos de
    # prueba; con la ganancia de emergencia (lazo de volumen de ~1 ms) RK45 podía oscilar entre quemar y
    # reinyectar dentro de un paso y dejar agua de más (hasta +12 % de H2 por tonelada en EXP-10, según el
    # tamaño exacto del escalón). Por eso el método por defecto es LSODA, que maneja la rigidez.
    modo = {'drops': 0, 'pre': None, 'fase2': False}

    def m_H2_equilibrio(m_agua, drops):
        carga_rest = (n_entregas - drops) * toneladas_por_drop * 1000.0
        return (p.m_vacio + carga_rest + m_agua) / (rho_h0 / rho_H2_ini - 1)

    def dinamica(t, y):
        v, h, m_H2, T_H2, Q_con, Q_gen, m_agua = y[:7]
        B_v = max(y[9], 0.0) if cond is not None else 0.0      # vapor retenido como gas (modo pulmón)
        m_H2   = max(m_H2, 0.1)
        m_agua = max(m_agua, 0.0)

        drops_hechos  = modo['drops']
        m_carga_actual = p.m_carga - drops_hechos * (toneladas_por_drop * 1000.0)
        # Lastre de agua (O2 atmosferico absorbido) participa en la inercia y en el
        # equilibrio de flotabilidad — igual que t2_gain_scheduled.py, segun §2.3.
        M_eff        = p.m_vacio + m_carga_actual + m_agua
        M_total_nave = M_eff + m_H2

        rho_aire = p.rho_aire_std * np.exp(-(h - p.h0) / p.H_escala)
        V_H2     = (m_H2 * p.R_H2 * T_H2) / p.P_atmos
        V_eq     = (M_eff * p.R_H2 * p.T1) / \
                   (p.P_atmos * (rho_h0 / rho_H2_ini - 1))
        V_gas    = V_H2 + (condensador.volumen_pulmon(B_v) if B_v > 0.0 else 0.0)
        M_I      = inercia.inercial(M_total_nave + B_v, rho_aire, V_casco, V_gas, k_nave)

        # Controlador PD con ganancia adaptativa: Kd con la inercia real (zeta de diseño)
        zeta    = p.zeta
        # (autor, 2026-10-09) ANTES: Kd_auto = 2 * zeta * np.sqrt(p.Kp * M_I / (rho_aire * p.g)) con zeta = 1.35, que
        # con la rigidez real (rho_aire - rho_H2) g daba un amortiguamiento de 1.302; ahora un solo zeta, el real.
        Kd_auto = 2 * zeta * np.sqrt(p.Kp * M_I / ((rho_aire - rho_H2_ini) * p.g))
        V_cmd   = V_eq + p.Kp * (p.h0 - h) - Kd_auto * v
        V_error = V_cmd - V_H2

        # Gain scheduling: Fase 1 (emergencia) → Fase 2 (trim); el cambio lo hace el evento de abajo
        K_mass_eff = p.K_mass_trim if modo['fase2'] else p.K_mass_control
        dm = K_mass_eff * V_error
        if drops_hechos == 0:
            dm = 0.0

        # Pulso proactivo — usa dm_pres[k] específico de cada drop
        if modo['pre'] is not None:
            dm = -dm_pres[modo['pre']]

        dm = np.clip(dm, -dm_opts[0], dm_opts[0])
        if not inyeccion:
            dm = min(dm, 0.0)   # sin reinyectar H2 desde los tanques: solo se quema

        F_emp = rho_aire * p.g * V_H2
        if B_v > 0.0:                                          # el vapor como gas ocupa volumen
            F_emp += rho_aire * p.g * condensador.volumen_pulmon(B_v)
        F_arr = 0.5 * p.Cd * DIRIGIBLE['Area'] * rho_aire * v**2 * np.sign(v)
        F_hel = y[i_hel] if helices is not None else 0.0
        dv    = (F_emp - (M_total_nave + B_v) * p.g - F_arr - F_hel) / M_I
        if helices is not None:
            ant = helices.get('anticipar')
            F_rampa = (max(0.0, ant['rampa'] * (t - (t_drops[0] - ant['t_previa'] - ant.get('adelanto', helices['tau']))))
                       if ant is not None else None)
            if drops_hechos >= 1:
                wn = helices['wn']
                F_est = F_emp - (M_total_nave + B_v) * p.g
                if F_rampa is not None:
                    F_est = min(F_est, F_rampa)
                F_cmd = F_est + M_I * (wn**2 * (h - p.h0) + 2 * helices['zeta'] * wn * v)
                F_cmd = min(max(F_cmd, -helices['F_max'] if helices.get('bidireccional', True) else 0.0),
                            helices['F_max'])
            else:
                F_cmd = min(F_rampa, helices['F_max']) if F_rampa is not None else 0.0
            dF = (F_cmd - F_hel) / helices['tau']
        dQgen  = abs(dm) * p.Q_comb_H2 * p.rend_gen if dm < -0.05 else 0.0
        extra = [dF] if helices is not None else []
        if cond is None:                                       # artículo: agua a bordo al instante
            d_agua = -(1 + p.ratio_O2_H2) * min(dm, 0.0) if retener else 0.0   # 1 kg H2 + 8 kg O2 -> 9 kg H2O
            return [dv, v, dm, (p.T1 - T_H2)/5.0, 0.0, dQgen, d_agua] + extra
        dP_c, d_agua, d_vent, dB_v = condensador.derivadas(max(-dm, 0.0), y[7], B_v, cond, t)
        return [dv, v, dm, (p.T1 - T_H2)/5.0, 0.0, dQgen, d_agua, dP_c, d_vent, dB_v] + extra

    y0  = [0.0, p.h0, m_H2_ini, p.T1, 0.0, 0.0, 0.0] + ([] if cond is None else
                                                       [condensador.potencia_inicial(cond), 0.0, 0.0])
    if helices is not None:
        y0 = y0 + [0.0]
    # Tolerancia CONVERGIDA (VER-02): rtol 1e-6 y atol POR ESTADO = rtol x escala del estado. Antes era
    # la de scipy (rtol 1e-3, atol 1e-6 único), con 1.2 % de error numérico en la excursión.
    if atol is None:
        atol = rtol * np.array(ESCALAS_ESTADO + (ESCALAS_COND if cond is not None else ())
                               + ((1e3,) if helices is not None else ()))

    def ev_trim(t, y):
        """Fase 1 -> fase 2: la masa de gas llega a la de equilibrio con la carga que queda."""
        return max(y[2], 0.1) - m_H2_equilibrio(max(y[6], 0.0), modo['drops'])
    ev_trim.terminal, ev_trim.direction = True, -1

    bordes = list(t_pre_starts) + list(t_drops)
    if cond is not None and cond.get('prearranque') and cond.get('t_arranque') is not None:
        bordes.append(cond['t_arranque'])
    bordes = sorted({b for b in bordes if 0.0 < b < t_final} | {t_final})
    t0, y_act = 0.0, np.array(y0, dtype=float)
    ts, ys = [np.array([t0])], [y_act.reshape(-1, 1)]
    fase2_at_drop = 0
    for tb in bordes:
        while t0 < tb:
            drops = sum(1 for td in t_drops if td <= t0)
            if drops != modo['drops']:
                if modo['fase2'] and drops > fase2_at_drop:
                    modo['fase2'] = False
                modo['drops'] = drops
                if verbose:
                    print(f"  [t={t0:6.1f}s] DROP #{drops}  carga restante: "
                          f"{(p.m_carga - drops * m_drop_kg) / 1000:.1f} t")
            pre = next((k for k in range(n_entregas) if t_pre_starts[k] <= t0 < t_drops[k]), None)
            if verbose and pre is not None and pre != modo['pre']:
                print(f"  [t={t0:6.2f}s] PRECOND #{pre+1} activado  dm={-dm_pres[pre]:.1f} kg/s  "
                      f"t_anticipo={t_anticipos[pre]:.3f}s  → drop en t={t_drops[pre]:.1f}s")
            modo['pre'] = pre
            if not modo['fase2'] and drops > 0 and ev_trim(t0, y_act) <= 0.0:
                modo['fase2'], fase2_at_drop = True, drops
            eventos = [ev_trim] if (not modo['fase2'] and drops > 0) else None
            tramo = solve_ivp(dinamica, [t0, tb], y_act, method=metodo, max_step=max_step, rtol=rtol,
                              atol=atol, events=eventos)
            if tramo.status < 0:
                raise RuntimeError(f'solve_ivp falló en [{t0:.3f}, {tb:.3f}] s: {tramo.message}')
            ts.append(tramo.t[1:])
            ys.append(tramo.y[:, 1:])
            t0, y_act = float(tramo.t[-1]), tramo.y[:, -1].copy()
            if tramo.status == 1:                              # evento: pasa a trim
                modo['fase2'], fase2_at_drop = True, drops
    sol = SimpleNamespace(t=np.concatenate(ts), y=np.concatenate(ys, axis=1))

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
                           dm_pre=dm_pres[0], t_anticipo=t_anticipos[0],
                           dm_opt=dm_opts[0], alpha_star=alpha_stars[0],
                           tons_p=toneladas_por_drop)
        else:
            _graficar_multi(t_arr, h_arr, v_arr, a_arr, dm_arr, m_H2_arr,
                            Q_gen_arr/3.6e6, p, n_entregas=n_entregas,
                            intervalo=intervalo, tons_p=toneladas_por_drop,
                            dm_pres=dm_pres, t_anticipos=t_anticipos,
                            t_pre_starts=t_pre_starts)

    dh_peak = float(np.max(np.abs(h_arr[t_arr >= intervalo] - p.h0)))

    # Métricas de toda la maniobra (precondicionamiento + liberación + asentamiento)
    dh_arr = h_arr - p.h0
    post   = t_arr >= t_drops[0]
    v_crit = 0.10 if p.m_vacio < 10000 else 0.05           # criterio de estabilización (§2.2)
    fuera  = post & ((np.abs(dh_arr) > h_max) | (np.abs(v_arr) >= v_crit))
    if not np.any(fuera):
        t_estab = 0.0
    else:
        ult = np.where(fuera)[0][-1]
        t_estab = float(t_arr[ult + 1] - t_drops[0]) if ult + 1 < len(t_arr) else float('nan')

    hel = {}
    if helices is not None:
        F_arr_h = sol.y[i_hel]
        rho_t = p.rho_aire_std * np.exp(-(h_arr - p.h0) / p.H_escala)
        P_hel = np.abs(F_arr_h) ** 1.5 / (helices['FM'] * np.sqrt(2.0 * rho_t * helices['A']))
        activa = np.abs(F_arr_h) > 0.01 * m_drop_kg * p.g
        t_hel = float(t_arr[np.where(activa)[0][-1]] - t_drops[0]) if np.any(activa) else 0.0
        hel = {'F_hel_max_N': float(np.max(F_arr_h)), 'F_hel_min_N': float(np.min(F_arr_h)),
               'P_hel_max_W': float(np.max(P_hel)), 'E_hel_J': float(np.trapezoid(P_hel, t_arr)),
               't_hel_s': t_hel, 'satura': bool(np.max(np.abs(F_arr_h)) >= helices['F_max'] * (1 - 1e-6))}
    return {**hel, 'Energia_kWh': Q_gen_arr[-1]/3.6e6,
            'h_max_m':     h_max, 'D_lim_m': D_lim,
            'dh_arriba_m': float(np.max(dh_arr)), 'dh_abajo_m': float(np.min(dh_arr)),
            'error_max_m': float(np.max(np.abs(dh_arr))),
            'v_max_ms':    float(np.max(np.abs(v_arr))),
            'a_max_ms2':   float(np.max(np.abs(a_arr))),
            't_estab_s':   t_estab,
            't_anticipo_teorico': t_ant_teoricos[0],
            'dm_pre':      dm_pres[0],
            't_anticipo':  t_anticipos[0],
            'alpha_star':  alpha_stars[0],
            'dm_opt':      dm_opts[0],
            'dm_opts':     dm_opts,
            'alpha_stars': alpha_stars,
            't_anticipos': t_anticipos,
            'M_total':     M_ini,
            'M_inercial':  I_eq,
            'dh_peak_m':   dh_peak,
            'm_H2_neto':   float(m_H2_arr[0] - m_H2_arr[-1]),
            'm_H2_bruto':  float(m_H2_arr[0] - np.min(m_H2_arr)),
            'm_agua_fin':  float(m_agua_arr[-1]),
            'm_agua_venteada': float(sol.y[8][-1]) if cond is not None else 0.0,
            'm_vapor_pulmon_fin': float(sol.y[9][-1]) if cond is not None else 0.0,
            **({'serie': {'t': t_arr, 'h': h_arr, 'v': v_arr, 'm_H2': m_H2_arr, 'agua': m_agua_arr,
                          'P_c': sol.y[7] if cond is not None else None,
                          'vent': sol.y[8] if cond is not None else None,
                          'pulmon': sol.y[9] if cond is not None else None,
                          'F_hel': sol.y[i_hel] if helices is not None else None,
                          't_drops': t_drops}} if devolver_serie else {})}


# =============================================================================
# Gráfica
# =============================================================================

def _graficar_base(t, h, v, a, dm, m_h2, E, p, t_lib,
                   dm_pre, t_anticipo, dm_opt, alpha_star, tons_p=1):

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
    ax_alt.set_xlim(0, t[-1])
    ax_alt.axhline(p.h0, color=COL_ALT, ls=':', alpha=0.5)
    ax_alt.axvline(t_lib, color=COL_LIB, ls='--', lw=1.5, label=f'Release ({t_lib} s)')

    subtitle = (
        f'$\\dot{{m}}_{{\\mathrm{{pre}}}}={dm_pre:.1f}$ kg/s'
        f'  |  $t_{{\\mathrm{{anticip}}}}={t_anticipo:.2f}$ s'
        f'  |  $\\alpha^*={alpha_star:.3f}$'
        f'  |  $\\dot{{m}}_{{\\mathrm{{opt}}}}={dm_opt:.2f}$ kg/s'
        f'  [Appendix A — per-drop]'
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
    nombre    = DIRIGIBLE['nombre'].replace(' ', '')
    tons_str  = f"{int(tons_p*1000)}kg" if tons_p < 1 else f"{int(tons_p)}t"
    path      = rutas.figura(f"minimax_single{tons_str}_{nombre}")
    plt.savefig(path, dpi=200, bbox_inches='tight')
    print(f"  Figura guardada: {path}")
    rutas.mostrar()


def _graficar_multi(t, h, v, a, dm, m_h2, E, p, n_entregas, intervalo,
                    tons_p, dm_pres, t_anticipos, t_pre_starts=None):
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
    flows_str = "  |  ".join(
        f"D{k+1}: $\\dot{{m}}$={dm_pres[k]:.2f} kg/s  $t_a$={t_anticipos[k]:.2f}s"
        for k in range(n_entregas)
    )
    _ = flows_str  # (sin título dentro de la figura)

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
        for tps in t_pre_starts:
            axes[2].axvline(tps, color='royalblue', ls=':', lw=1.4, alpha=0.8)

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
    nombre   = DIRIGIBLE['nombre'].replace(' ', '')
    tons_str = f"{int(tons_p*1000)}kg" if tons_p < 1 else f"{int(tons_p)}t"
    path     = rutas.figura(f"minimax_multi{n_entregas}x{tons_str}_{nombre}")
    plt.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Figura guardada: {path}")
    rutas.mostrar()


# =============================================================================
# Entrada principal
# =============================================================================

if __name__ == '__main__':
    plt.switch_backend('Agg')

    # ── Tabla: single 1-tonne release para las 4 plataformas (Tabla 4) ───────
    print("\n" + "="*78)
    print("  QUEMA LIMPIA — single 1 t  (tab:multiplatform_validation)")
    print("="*78)
    print(f"  {'Platform':<22}{'M_total':>10}{'alpha*':>9}{'m*_H2':>9}"
          f"{'t_ant':>8}{'E_gen':>8}{'dh_peak':>9}{'H2_neto':>9}")
    tabla = {}
    for nom, d in p_base.PLATAFORMAS.items():
        usar_nave(nom, 1000.0)
        # la figura single del caso base la genera drop_flujo_optimo.py (minimax_base_<nave>)
        r = correr(n_entregas=1, toneladas_por_drop=1.0, graficar=False, verbose=False)
        tabla[nom] = {'corto': d['corto'], **{k: r[k] for k in (
            'M_total', 'alpha_star', 'dm_opt', 't_anticipo', 'Energia_kWh', 'dh_peak_m', 'm_H2_neto')}}
        print(f"  {nom:<22}{r['M_total']:>10.0f}{r['alpha_star']:>9.4f}"
              f"{r['dm_opt']:>9.2f}{r['t_anticipo']:>8.2f}"
              f"{r['Energia_kWh']:>8.0f}{r['dh_peak_m']:>9.2f}{r['m_H2_neto']:>9.1f}")

    # ── Multi-drop (energias, §3.1.3): 4 entregas; Zeppelin NT de 250 kg (su capacidad es 1 t) ──
    print("\n" + "="*78)
    print("  MULTI-DROP  (energias acumuladas, §3.1.3)")
    print("="*78)
    multi = {}
    for nom, d in p_base.PLATAFORMAS.items():
        n, tpd = 4, (0.25 if nom == 'Zeppelin NT' else 1.0)
        usar_nave(nom, 1000.0)
        graf = (nom == p_base.nombre_dirigible)   # figura del paper: el caso base
        r = correr(n_entregas=n, toneladas_por_drop=tpd, graficar=graf, verbose=False)
        multi[nom] = {'corto': d['corto'], 'n_entregas': n, 't_por_entrega': tpd,
                      'E_cum_kWh': r['Energia_kWh'], 'm_H2_neto': r['m_H2_neto'],
                      'dh_peak_m': r['dh_peak_m'], 'dm_opts': [float(x) for x in r['dm_opts']]}
        tag = f"{nom} ({n}x{int(tpd*1000)}kg)"
        print(f"  {tag:<26} E_cum={r['Energia_kWh']:>7.0f} kWh | "
              f"H2_neto={r['m_H2_neto']:>6.1f} kg | dh_peak={r['dh_peak_m']:.2f} m | "
              f"m*(drops)={[round(float(x), 2) for x in r['dm_opts']]}")

    rutas.guardar_datos('minimax_plataformas', {'D_lim_m': D_LIM, 'h_max_m': H_MAX_PREC,
                                                'single_1t': tabla, 'multidrop': multi})
