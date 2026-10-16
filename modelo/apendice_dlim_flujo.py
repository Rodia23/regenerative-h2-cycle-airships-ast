# -*- coding: utf-8 -*-
# [EN] Clean-burn rate versus the allowed sink D_lim for the four airships (Section 4).
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
# =============================================================================
# apendice_dlim_flujo.py
# D_lim vs flujo instalado de la quema limpia (minimax.limpio) para los 4 dirigibles.
# La quema limpia no depende de h_max: ṁ* = K_pre·α*^(3/2) ∝ D_lim^(-1/2), α* = 1/(1+s).
# La nave termina (1+s)·D_lim abajo, así que D_lim ≤ h_max/(1+s) mantiene la nave en el umbral.
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import plot_style
import parametros as p
import inercia
import minimax as mm
import rutas

h_max      = p.h_max_precision      # umbral de PRECISIÓN (|Δh|, §2.4): solo para la cota con reinyección
m_cargo    = float(p.m_carga)
D_lim_oper = p.D_lim_oper           # punto de operación de las simulaciones

# ── Dirigibles ───────────────────────────────────────────────────────────────
# Color: paleta Okabe-Ito (CVD-safe, validada). Mismo esquema que
# t3_barrido_reactivo.py (Fig 4): ZR-3=vermillón, FW=azul, Hind=verde, Zepp=oro.
DIRIGIBLES = [{'nombre': n, 'color': d['color'], 'corto': d['corto']} for n, d in p.PLATAFORMAS.items()]
for d in DIRIGIBLES:
    eq = mm.equilibrio(d['nombre'], m_cargo)
    d['M']   = eq['M_total']
    d['M_I'] = eq['M_I']                               # (1 + k)·rho·V_casco, con masa añadida (inercia.py)

# Grilla D_lim: escala log de 0.005 m hasta h_max = 1 m
D_arr = np.logspace(np.log10(5e-3), np.log10(h_max), 800)

# ── Figura ───────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8.0, 5.2))
for d in DIRIGIBLES:
    dm_arr = np.array([mm.limpio(d['M_I'], m_cargo, D)['dm_star'] for D in D_arr])
    ax.plot(D_arr, dm_arr, color=d['color'], lw=2.3, label=d['nombre'], zorder=3)

# ── D_lim con el que la bajada total (1+s)·D_lim llega al umbral de precisión ──
ax.axvline(h_max / 2.0, color='black', lw=1.5, ls='--', alpha=0.5, zorder=1)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(5e-3, h_max)
ax.set_xticks([0.01, 0.05, 0.1, 0.5, 1.0])
ax.set_xticklabels(['0.01', '0.05', '0.1', '0.5', '1.0'])
ax.set_xlabel(r'Max. preconditioning depth  $D_\mathrm{lim}$  [m]', fontsize=10)
ax.set_ylabel(r'Installed burn rate  $\dot{m}_{H_2}^{*}$  [kg/s]', fontsize=10)
# (sin título dentro de la figura: va en el caption del manuscrito)
ax.grid(True, which='both', alpha=0.25)
ax.tick_params(labelsize=9)
ax.legend(fontsize=8.5, loc='upper right', framealpha=0.93)

plt.tight_layout()
fname = rutas.figura('dlim_flujo')
plt.savefig(fname, dpi=180, bbox_inches='tight')
rutas.mostrar()
print(f"Guardado: {fname}")

# ── Tabla consola ────────────────────────────────────────────────────────────
print(f"\n  Quema limpia, D_lim = {D_lim_oper} m (ley general: α* = 1/2; exacta: α* = 1/(1+s))")
print(f"  {'Dirigible':<24} {'α*':>7} {'ṁ* [kg/s]':>10} {'ley gral.':>10} {'bajada total [m]':>17} "
      f"{'cota reiny. 1 m':>16} {'sobrequema':>11}")
for d in DIRIGIBLES:
    r = mm.resolver(d['nombre'], m_cargo, h_max, D_lim_oper)
    d.update({'alpha': r['alpha'], 'dm_star': r['dm_star'], 'dm_star_aprox': r['dm_star_aprox'],
              't_anticip': r['t_anticip'], 't_reactiva': r['t_reactiva'],
              'hundimiento_total': r['hundimiento_total'], 'dm_cota': r['dm_cota'],
              'sobrequema_cota': r['sobrequema_cota']})
    print(f"  {d['nombre']:<24} {d['alpha']:>7.4f} {d['dm_star']:>10.3f} {d['dm_star_aprox']:>10.3f} "
          f"{d['hundimiento_total']:>17.3f} {d['dm_cota']:>16.3f} {100 * d['sobrequema_cota']:>10.0f}%")

rutas.guardar_datos('dlim_flujo', {
    'D_lim_oper_m': D_lim_oper, 'h_max_m': h_max, 'm_carga_kg': m_cargo,
    'alpha_star': 0.5,                              # ley general (s = 1)
    'plataformas': {d['nombre']: {'corto': d['corto'], 'M_total_kg': d['M'], 'M_inercial_kg': d['M_I'],
                                  'alpha': d['alpha'], 'dm_star_kgps': d['dm_star'],
                                  'dm_star_aprox_kgps': d['dm_star_aprox'], 't_anticip_s': d['t_anticip'],
                                  't_reactiva_s': d['t_reactiva'], 'hundimiento_total_m': d['hundimiento_total'],
                                  'dm_cota_reinyeccion_kgps': d['dm_cota'],
                                  'sobrequema_cota': d['sobrequema_cota']}
                    for d in DIRIGIBLES},
})
