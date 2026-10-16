# -*- coding: utf-8 -*-
# [EN] Sweep of the reactive controller (Section 4, Table 3).
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# t3_barrido_reactivo.py   (reemplaza a t3_barrido.py)
#
# Barrido paramétrico de ṁ_max para las 4 plataformas: caso REACTIVO PURO.
# El controlador PD de dos fases responde SOLO tras el drop (sin pulso proactivo)
# y la combustión reduce la flotabilidad SIN absorción de O₂ atmosférico.
# Este es el caso que acota la fórmula cerrada eq:mmin_final (§2.5).
#
# Salida:
#   - t3_desplazamiento.png : figura ÚNICA (fusiona las antiguas Fig 4 y Fig 5).
#       |Δh| pico vs ṁ_H2  +  franjas tecnológicas  +  marcas de t_peak.
#   - consola: datos de la Tabla 3 (ṁ_min, ṁ_fine, |Δh|, t_estab)  y  el par
#       t_peak / t_estab en cada punto marcado (para la mención en el texto).
# =============================================================================

import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter, NullFormatter
import plot_style  # noqa  (estilo del paper)
import parametros as p
import inercia
import rutas

# ── Plataformas ──────────────────────────────────────────────────────────────
# Color: paleta Okabe-Ito (CVD-safe, validada; ΔE_CVD 11, ΔE_normal 15.6).
# Mismo esquema en apendice_dlim_flujo.py (Fig 3). Las 2 curvas que se apiñan
# (PF3 y FW, M_total parecidos) llevan el par más contrastante: vermillón / azul.
_ORDEN_DIBUJO = ['ZR-3 Los Angeles', 'Flying Whales LCA60T', 'Zeppelin NT', 'LZ 129 Hindenburg']  # orden de la leyenda
DIRIGIBLES = [{'nombre': n, 'Area': p.PLATAFORMAS[n]['Area'], 'm_vacio': p.PLATAFORMAS[n]['m_vacio'],
               'V_casco': p.PLATAFORMAS[n]['V_casco'], 'k': p.PLATAFORMAS[n]['k'],
               'm_carga': p.m_carga, 'color': p.PLATAFORMAS[n]['color'], 'corto': p.PLATAFORMAS[n]['corto']}
              for n in _ORDEN_DIBUJO]

# Grid de evaluación para ṁ_min/ṁ_fine (Tabla 3): "el conjunto evaluado".
# Espaciado unitario hasta 35 kg/s (cubre ṁ_min y ṁ_fine de las 3 grandes:
# 5/8/9 y 18/27/32); más grueso arriba. El Zeppelin agrega EXTRA con espaciado
# unitario 36–50 (su ṁ_min ≈ 37) y grid fino 124–142 (su ṁ_fine ≈ 138).
TASAS_BASE       = list(range(1, 36)) + [38, 40, 45, 50, 60, 75, 100]
TASAS_ZEPP_EXTRA = list(range(36, 51)) + [55, 65, 80, 100, 120,
                    124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 145, 150, 160]

# Grid denso, solo para dibujar curvas lisas (NO interviene en la Tabla 3).
GRID_CURVA = np.logspace(0, np.log10(175), 60)

MARCA_ALTO  = 3      # ṁ del punto "tramo alto"
MARCA_FINAL = 160    # ṁ del punto "final"

FRANJAS_TECH = [
    (0.5,  5.0,   'seagreen',   0.13, 'H$_2$ combustion engine  ($\\lesssim$5 kg/s)'),
    (5.0,  50.0,  'darkorange', 0.13, 'Industrial system  (5–50 kg/s)'),
    (50.0, 300.0, 'firebrick',  0.09, 'Theoretical control limit  ($\\gtrsim$50 kg/s)'),
]

CRITERIO_H_M      = 1.0     # umbral de "se desestabilizó" (destino: |Δh| < 1 m)
CRITERIO_V_PESADO = 0.05    # m/s — umbral de velocidad residual (plataformas grandes)
CRITERIO_V_LIGERO = 0.10    # m/s — Zeppelin NT

# =============================================================================
# Predicción analítica cerrada (§2.5, eq:mmin_final):
#   |Δh|_pico = (2/3)·(ρ_H2/ρ_air)²·m_load³·g / (M_I·ṁ²)   ∝ ṁ⁻²,  M_I = (1+k)·ρ_air·V_casco (inercia.py)
# Asume quema CONSTANTE a ṁ (actuador ideal). El barrido con PD saturado la
# sigue mientras el controlador satura (ṁ ≲ ṁ_min) y se despega arriba.
# =============================================================================

def dh_analitico(diri, dm):
    ratio = (p.P_atmos / (p.R_H2 * p.T1)) / p.rho_aire_std          # ρ_H2/ρ_air
    M_I = inercia.inercial_casco(diri['V_casco'], diri['k'])         # inercia en equilibrio
    C = (2.0 / 3.0) * ratio**2 * diri['m_carga']**3 * p.g / M_I
    return C / np.asarray(dm, float)**2


# =============================================================================
# Motor de física — reactivo puro (solo pérdida de flotabilidad, sin O₂)
# Controlador PD de dos fases (emergencia K_m=60 → trim K_m=0.5), saturado a ṁ_max.
# =============================================================================

def evaluar_caso(diri, tasa_max):
    m_est       = diri['m_vacio']
    Area        = diri['Area']
    m_carga_ini = diri['m_carga']
    es_ligero   = m_est < 10000

    rho_h0     = p.rho_aire_std   # densidad del aire en la altitud de referencia h0
    rho_H2_ini = p.P_atmos / (p.R_H2 * p.T1)
    m_H2_ini   = (m_est + m_carga_ini) / (rho_h0 / rho_H2_ini - 1)
    V_casco    = diri['V_casco']                           # volumen total del casco (dato; inercia.py)

    t_drop, t_final = 30.0, 250.0

    # Fase 1 (emergencia, K_mass grande): respuesta rápida post-drop.
    # Fase 2 (trim, K_mass=0.5): asentamiento suave al llegar al nuevo equilibrio.
    fase2 = [False]

    def dinamica(t, y):
        v, h, m_H2, T_H2, m_agua = y
        m_H2   = max(m_H2,   0.1)
        m_agua = max(m_agua, 0.0)

        ha_soltado     = t >= t_drop
        m_carga_actual = 0.0 if ha_soltado else m_carga_ini
        M_total_nave   = m_est + m_carga_actual + m_H2 + m_agua

        zeta     = p.zeta
        rho_aire = p.rho_aire_std * np.exp(-(h - p.h0) / p.H_escala)

        M_eff   = m_est + m_carga_actual + m_agua
        V_H2    = (m_H2 * p.R_H2 * T_H2) / p.P_atmos
        M_I     = inercia.inercial(M_total_nave, rho_aire, V_casco, V_H2, diri['k'])   # inercia real
        # (autor, 2026-10-09) ANTES: Kd_auto = 2 * zeta * np.sqrt(p.Kp * M_I / (rho_aire * p.g)) con zeta = 1.35, que
        # con la rigidez real (rho_aire - rho_H2) g daba un amortiguamiento de 1.302; ahora un solo zeta, el real.
        Kd_auto = 2 * zeta * np.sqrt(p.Kp * M_I / ((rho_aire - rho_H2_ini) * p.g))

        V_eq = (M_eff * p.R_H2 * p.T1) / (p.P_atmos * (rho_h0 / rho_H2_ini - 1))

        V_cmd   = V_eq + p.Kp * (p.h0 - h) - Kd_auto * v
        V_error = V_cmd - V_H2

        if ha_soltado and not fase2[0]:
            m_H2_eq_nuevo = M_eff / (rho_h0 / rho_H2_ini - 1)
            if m_H2 <= m_H2_eq_nuevo:
                fase2[0] = True

        K_mass_eff = p.K_mass_trim if fase2[0] else p.K_mass_control
        dm = np.clip(K_mass_eff * V_error, -tasa_max, tasa_max)

        d_agua = 0.0   # SIN absorción de O₂
        F_emp = rho_aire * p.g * V_H2
        F_arr = 0.5 * p.Cd * Area * rho_aire * v**2 * np.sign(v)
        dv    = (F_emp - M_total_nave * p.g - F_arr) / M_I
        dT    = (p.T1 - T_H2) / 5.0
        return [dv, v, dm, dT, d_agua]

    y0  = [0.0, p.h0, m_H2_ini, p.T1, 0.0]
    sol = solve_ivp(dinamica, [0, t_final], y0, method='RK45', max_step=0.2)

    t, v, h = sol.t, sol.y[0], sol.y[1]
    m = t >= t_drop
    t_post, h_post, v_post = t[m] - t_drop, h[m], v[m]
    dh_post = h_post - p.h0

    h_desp = float(np.max(np.abs(dh_post)))
    t_peak = float(t_post[int(np.argmax(np.abs(dh_post)))])

    crit_v = CRITERIO_V_LIGERO if es_ligero else CRITERIO_V_PESADO
    se_desestab = (np.max(np.abs(dh_post)) >= CRITERIO_H_M) or (np.max(np.abs(v_post)) >= crit_v)
    if not se_desestab:
        t_estab = 0.0
    else:
        mov = np.abs(v_post) >= crit_v
        if not np.any(mov):
            t_estab = 0.0
        else:
            li = np.where(mov)[0][-1]
            t_estab = float(t_post[li + 1]) if li + 1 < len(t_post) else 250.0

    return round(h_desp, 3), round(t_peak, 1), round(t_estab, 1)


# =============================================================================
# Barrido
# =============================================================================

def barrer():
    for d in DIRIGIBLES:
        scan = sorted(set(TASAS_BASE + (TASAS_ZEPP_EXTRA if d['nombre'] == 'Zeppelin NT' else [])))

        rows = {}
        print(f"  {d['nombre']}")
        for tm in scan:
            hd, tp, te = evaluar_caso(d, tm)
            rows[tm] = (hd, tp, te)
            print(f"    ṁ={tm:>4} kg/s  |Δh|={hd:>7.3f} m  t_peak={tp:>6.1f} s  t_estab={te:>6.1f} s")

        d['rows'] = rows
        # ṁ_min / ṁ_fine: primer valor del CONJUNTO EVALUADO que cumple el umbral
        d['mmin']  = next((tm for tm in scan if rows[tm][0] <= p.h_max_seguridad), None)
        d['mfine'] = next((tm for tm in scan if rows[tm][0] <= p.h_max_precision), None)

        # puntos marcados que no están en el scan → evaluar aparte
        for tm in (MARCA_ALTO, MARCA_FINAL):
            if tm not in rows:
                rows[tm] = evaluar_caso(d, tm)

        # curva lisa: grid denso, solo |Δh|
        d['curva'] = np.array([evaluar_caso(d, float(tm))[0] for tm in GRID_CURVA])
    return DIRIGIBLES


def imprimir_tablas(dirs):
    print("\n" + "=" * 78)
    print("  DATOS PARA LA TABLA 3  (tab:minimos_actualizada)")
    print("=" * 78)
    print(f"  {'Plataforma':<21}{'ṁ_min':>7}{'|Δh|@min':>10}{'ṁ_fine':>8}{'|Δh|@fine':>11}{'t_estab@fine':>14}")
    for d in dirs:
        mm, mp = d['mmin'], d['mfine']
        print(f"  {d['nombre']:<21}{mm:>7}{d['rows'][mm][0]:>10.2f}{mp:>8}"
              f"{d['rows'][mp][0]:>11.3f}{d['rows'][mp][2]:>14.1f}")

    # NOTA: esto es |Δh|_sim vs eq:peak EVALUADOS AL MISMO ṁ_min (métrica de |Δh|).
    # El "6-19% / 24%" del manuscrito es OTRA cosa: ṁ_min_sim vs ṁ_min_analítico
    # (métrica de tasa). Ej. Zeppelin: sim 37 vs analítico 29.8 → +24% (tasa),
    # mientras que en |Δh| a ṁ=37 da +53%. Son comparaciones distintas.
    print("\n  Desviación sim vs analítica |Δh| ∝ ṁ⁻² (a ṁ_min fijo; NO es el % del paper):")
    for d in dirs:
        mm = d['mmin']
        an_min = float(dh_analitico(d, mm))
        print(f"    {d['nombre']:<21} @ṁ_min={mm:>3}: sim {d['rows'][mm][0]:.2f} m  vs  "
              f"analítico {an_min:.2f} m  → +{100*(d['rows'][mm][0]/an_min - 1):.0f}%")

    print("\n" + "=" * 78)
    print("  PARES  t_peak / t_estab  EN CADA PUNTO MARCADO  (para el texto §3.1.2)")
    print("=" * 78)
    for d in dirs:
        mm, mp, r = d['mmin'], d['mfine'], d['rows']
        print(f"\n  {d['nombre']}")
        print(f"    tramo alto  (ṁ={MARCA_ALTO:>3})  : t_peak={r[MARCA_ALTO][1]:>6.1f}  t_estab={r[MARCA_ALTO][2]:>6.1f}")
        print(f"    ṁ_min       (ṁ={mm:>3})  : t_peak={r[mm][1]:>6.1f}  t_estab={r[mm][2]:>6.1f}")
        print(f"    ṁ_fine      (ṁ={mp:>3})  : t_peak={r[mp][1]:>6.1f}  t_estab={r[mp][2]:>6.1f}")
        print(f"    final       (ṁ={MARCA_FINAL:>3})  : t_peak={r[MARCA_FINAL][1]:>6.1f}  t_estab={r[MARCA_FINAL][2]:>6.1f}")


# =============================================================================
# Figura única
# =============================================================================

def _suave(dh):
    """|Δh| pico es monótona decreciente en ṁ; limpia ruido numérico + suaviza."""
    a = np.minimum.accumulate(np.asarray(dh, float))
    return np.exp(savgol_filter(np.log(a), window_length=11, polyorder=2))


def graficar(dirs):
    fig, ax = plt.subplots(figsize=(9.0, 5.9))

    for x0, x1, col, alpha, _ in FRANJAS_TECH:
        ax.axvspan(x0, x1, color=col, alpha=alpha, zorder=0, lw=0)

    ax.axhline(5.0, ls='--', c='0.3', lw=1.1, zorder=1)
    ax.axhline(1.0, ls=':',  c='0.3', lw=1.1, zorder=1)
    ax.text(1.12, 5.7,  'safety 5 m',    ha='left', va='bottom', fontsize=8, color='0.3')
    ax.text(1.12, 1.13, 'precision 1 m', ha='left', va='bottom', fontsize=8, color='0.3')

    dm_c = np.array(GRID_CURVA, float)
    for d in dirs:
        dh_c = _suave(d['curva'])
        ax.plot(dm_c, dh_c, '-', color=d['color'], lw=2.3, label=d['nombre'], zorder=3)

        puntos = [MARCA_ALTO, d['mmin'], d['mfine']]
        if d['nombre'] != 'Zeppelin NT':
            puntos.append(MARCA_FINAL)                # Zeppelin la omite: pega con ṁ_fine

        # Colocación de etiquetas — las 3 curvas grandes se apiñan (PF3 > FW > Hind):
        #   ZR-3       → etiqueta ARRIBA (offset amplio)
        #   Hindenburg → etiqueta ABAJO (final: a la izquierda del punto)
        #   FW (medio) → ṁ_min/ṁ_fine debajo del punto; extremos (alto/final) a la izquierda
        #   Zeppelin (aislada) → ARRIBA centrada
        for tm in puntos:
            y  = float(np.exp(np.interp(np.log(tm), np.log(dm_c), np.log(dh_c))))
            tp = d['rows'][tm][1]
            etiqueta = True
            if d['nombre'] == 'ZR-3 Los Angeles':
                va, oy = 'bottom', 12
                ha, ox = ('right', -9) if tm == MARCA_FINAL else ('left', 9)
            elif d['nombre'] == 'LZ 129 Hindenburg':
                if tm == MARCA_FINAL:                        # a la izq., centrado en el punto
                    va, oy, ha, ox = 'center', 0, 'right', -9
                else:
                    va, oy, ha, ox = 'top', -12, 'left', 9
            elif d['nombre'] == 'Flying Whales LCA60T':
                if tm in (MARCA_ALTO, MARCA_FINAL):          # extremos: a la izq. del punto
                    va, oy, ha, ox = 'center', 0, 'right', -8
                else:                                        # umbrales: debajo del punto
                    va, oy, ha, ox = 'top', -13, 'left', 7
            else:  # Zeppelin NT
                va, oy = 'bottom', 11
                ha, ox = 'center', 0
            ax.plot([tm], [y], 'o', ms=8, color=d['color'], mec='white', mew=1.3, zorder=5)
            if etiqueta:
                ax.annotate(f"{tp:.1f} s", (tm, y), textcoords='offset points',
                            xytext=(ox, oy), ha=ha, va=va, fontsize=11, color='black', zorder=6)

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(1, 178); ax.set_ylim(0.07, 340)
    ax.set_xticks([1, 2, 3, 5, 10, 20, 50, 100, 150])
    ax.set_yticks([0.1, 0.2, 0.5, 1, 2, 3, 5, 10, 20, 50, 100, 200])
    g = FuncFormatter(lambda v, _pos: f'{v:g}')
    ax.xaxis.set_major_formatter(g); ax.yaxis.set_major_formatter(g)
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'H$_2$ combustion rate  $\dot{m}_{H_2}$  [kg/s]', fontsize=10)
    ax.set_ylabel(r'Max. vertical displacement  $|\Delta h|$  [m]', fontsize=10)
    ax.grid(True, which='major', alpha=0.22)

    leg = ax.legend(fontsize=8.5, loc='lower left', framealpha=0.95)
    ax.add_artist(leg)
    parches = [mpatches.Patch(color=c, alpha=a + 0.18, label=lbl)
               for _, _, c, a, lbl in FRANJAS_TECH]
    ax.legend(handles=parches, fontsize=7.6, loc='upper right',
              framealpha=0.95, title='Actuator technology')

    plt.tight_layout()
    fpath = rutas.figura('barrido_reactivo')
    fig.savefig(fpath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Figura guardada: {fpath}")


# =============================================================================
if __name__ == '__main__':
    print("=" * 78)
    print("  BARRIDO REACTIVO PURO — sin O₂, sin precondicionamiento")
    print("=" * 78 + "\n")
    dirs = barrer()
    imprimir_tablas(dirs)
    graficar(dirs)
    rutas.guardar_datos('barrido_reactivo', {d['nombre']: {
        'corto': d['corto'], 'dm_min_kgps': d['mmin'], 'dm_fine_kgps': d['mfine'],
        'dh_en_min_m': d['rows'][d['mmin']][0], 'dh_en_fine_m': d['rows'][d['mfine']][0],
        't_estab_fine_s': d['rows'][d['mfine']][2],
        'dh_analitico_en_min_m': float(dh_analitico(d, d['mmin'])),
        'puntos': {str(tm): {'dh_m': r[0], 't_peak_s': r[1], 't_estab_s': r[2]} for tm, r in sorted(d['rows'].items())},
    } for d in dirs})
