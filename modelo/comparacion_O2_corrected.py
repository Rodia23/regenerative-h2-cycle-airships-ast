# -*- coding: utf-8 -*-
# [EN] Atmospheric versus stored oxygen: hydrogen per tonne and regeneration (Sections 3.3 and 6).
"""
comparacion_O2_corrected.py
============================
Comparison of O2 source (on-board tanks vs. atmospheric air) for descent-phase
H2 consumption after a 1-tonne payload release.

PHYSICAL MODEL — key distinction between the two scenarios:

  Scenario A — O2 from on-board tanks (closed system):
    When H2 burns with on-board O2:
      - H2 leaves ballonets       →  volume decreases → buoyancy ↓
      - O2 leaves storage tank    →  becomes H2O water ballast (stays on board)
      - H2O retained on board     →  mass conserved, ΔM_vehicle = 0
    Net downward force per kg H2 burned = (ρ_air / ρ_H2) × g
    → factor_A = ρ_air / ρ_H2

  Scenario B — O2 from atmosphere (open system, this work):
    When H2 burns with atmospheric O2:
      - H2 leaves ballonets       →  volume decreases → buoyancy ↓
      - 8 kg O2 absorbed from air →  becomes H2O water ballast → ΔM_vehicle = +8·ṁ
    Net downward force per kg H2 burned = (ρ_air / ρ_H2 + 8) × g  =  L_eff × g
    → factor_B = ρ_air / ρ_H2 + 8  =  L_eff

NOTE — what the original comparacion_O2.py computed (incorrectly for Scenario A):
    factor_A_wrong = (ρ_air - ρ_H2) / ρ_H2 = ρ_air/ρ_H2 - 1
    This corresponds to VENTING H2 to atmosphere (both buoyancy and H2 weight are
    lost), not to combustion where the H2 mass is conserved as on-board water.
"""

import numpy as np
import matplotlib.pyplot as plt
import plot_style
import parametros as p
import rutas

# =============================================================================
# CONDITIONS — ISA sea level (parametros.py; coastal operation, Lima)
# =============================================================================
P_STP        = p.P_atmos        # Pa    — ISA sea-level pressure
T_STP        = p.T1             # K     — ISA sea-level temperature (15 °C)
rho_air_stp  = p.rho_aire_std   # kg/m³ — ISA sea-level air density
rho_H2_stp   = P_STP / (p.R_H2 * T_STP)   # ≈ 0.0853 kg/m³

RATIO_O2_H2  = p.ratio_O2_H2   # stoichiometric O2/H2 mass ratio

# =============================================================================
# EFFECTIVE FORCE FACTORS  (net downward force per kg H2 burned / g)
# =============================================================================
factor_A = rho_air_stp / rho_H2_stp           # O2 from tanks: buoyancy only
factor_B = rho_air_stp / rho_H2_stp + RATIO_O2_H2   # air: buoyancy + absorbed O2 mass

# =============================================================================
# MINIMUM H2 REQUIRED  (to cancel m_load × g)
# =============================================================================
def m_H2_required(m_load, factor):
    return m_load / factor

m_load = 1000.0   # kg payload

mA = m_H2_required(m_load, factor_A)
mB = m_H2_required(m_load, factor_B)

reduction_pct = (mA - mB) / mA * 100.0
savings_kg    = mA - mB

sep = "=" * 50
print(sep)
print("O2 SOURCE COMPARISON - corrected model")
print(sep)
print(f"  rho_air / rho_H2    = {rho_air_stp/rho_H2_stp:.4f}")
print(f"  L_eff  (air)        = {factor_B:.4f}")
print()
print(f"  Scenario A (O2 tanks):  factor = {factor_A:.4f}")
print(f"    m_H2 required = {mA:.2f} kg")
print()
print(f"  Scenario B (air):       factor = {factor_B:.4f}")
print(f"    m_H2 required = {mB:.2f} kg")
print()
print(f"  Reduction (A->B): {reduction_pct:.1f}%  (saves {savings_kg:.1f} kg H2 per delivery)")
print()
agua_A, agua_B = (1 + RATIO_O2_H2) * mA, (1 + RATIO_O2_H2) * mB
Egen_A = mA * p.Q_comb_H2 * p.rend_gen / 3.6e6
Egen_B = mB * p.Q_comb_H2 * p.rend_gen / 3.6e6
print(f"  Water produced (A): {agua_A:.1f} kg/delivery")
print(f"  Water produced (B): {agua_B:.1f} kg/delivery")
print()
print(f"  E_gen (A): {Egen_A:.0f} kWh/delivery")
print(f"  E_gen (B): {Egen_B:.0f} kWh/delivery")
print()

G_ef = rutas.leer_datos('irradiancia_casco')['estaciones']['anual']['efectiva']   # sobre el casco, anual
print(f"  Solar regeneration (Lima, {G_ef:.2f} kWh/m2/day on the hull, annual; target = Scenario B):")
print(f"  {'Platform':<20} {'reg(kg/d)':>10} {'days':>8} {'hours':>8}")
regen = {}
for name, d in p.PLATAFORMAS.items():
    reg  = G_ef * p.eta_panel * p.A_paneles_fraccion * d['A_total'] / p.E_H2_electrolisis
    days = mB / reg
    regen[name] = {'corto': d['corto'], 'A_total_m2': d['A_total'], 'regen_kg_dia': reg,
                   'dias': days, 'horas': days * 24}
    print(f"  {name:<20} {reg:>10.0f} {days:>8.2f} {days*24:>8.1f}")

print(sep)

# =============================================================================
# FIGURE
# =============================================================================
FIGSIZE   = (7.5, 5.0)
DPI       = 300
FS_LABEL  = 11
FS_LEGEND = 10
FS_TITLE  = 13
LW_MAIN   = 2.5
LW_DASH   = 1.5
COLOR_A   = '#1f77b4'   # blue  — O2 tanks
COLOR_B   = '#d62728'   # red   — air

n_pts    = 300
x_max    = mA * 1.5
m_H2_arr = np.linspace(0, x_max, n_pts)

deficit_A = factor_A * m_H2_arr
deficit_B = factor_B * m_H2_arr

fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
ax.set_facecolor('#fcfcfc')
ax.grid(True, linestyle='-', alpha=0.3, color='gray')

ax.plot(m_H2_arr, deficit_A, color=COLOR_A, lw=LW_MAIN,
        label=f'Scenario A — O$_2$ from tanks   ($m_{{H_2}}$ = {mA:.1f} kg)')
ax.plot(m_H2_arr, deficit_B, color=COLOR_B, lw=LW_MAIN,
        label=f'Scenario B — O$_2$ from air      ($m_{{H_2}}$ = {mB:.1f} kg)')

ax.axhline(m_load, color='black', ls='--', lw=LW_DASH,
           label=f'Payload $m_{{\\mathrm{{load}}}}$ = {m_load:.0f} kg')

ax.axvline(mA, color=COLOR_A, ls=':', lw=1.8)
ax.axvline(mB, color=COLOR_B, ls=':', lw=1.8)

# double-headed arrow showing savings
arrow_y = m_load * 1.12
ax.annotate('', xy=(mA, arrow_y), xytext=(mB, arrow_y),
            arrowprops=dict(arrowstyle='<->', color='green', lw=2.0))
ax.text((mA + mB) / 2, arrow_y + m_load * 0.04,
        f'{savings_kg:.1f} kg saved ({reduction_pct:.1f}%)',
        color='green', fontsize=FS_LEGEND, ha='center', va='bottom', fontweight='bold')

ax.set_xlabel('H$_2$ burned [kg]', fontsize=FS_LABEL, labelpad=10)
ax.set_ylabel('Net buoyancy deficit [kg equiv.]', fontsize=FS_LABEL, labelpad=10)

ax.set_xlim(0, x_max)
ax.set_ylim(0, m_load * 1.8)
ax.legend(fontsize=FS_LEGEND, loc='upper left', framealpha=0.95, edgecolor='gray')

plt.tight_layout()

output_path = rutas.figura('deficit_O2')
plt.savefig(output_path, bbox_inches='tight')
print(f"\nFigure saved: {output_path}")
rutas.mostrar()

rutas.guardar_datos('fuente_O2', {
    'rho_aire_sobre_rho_H2': factor_A, 'L_eff': factor_B, 'm_carga_kg': m_load,
    'm_H2_tanques_kg': mA, 'm_H2_aire_kg': mB, 'reduccion_pct': reduction_pct, 'ahorro_kg': savings_kg,
    'agua_tanques_kg': agua_A, 'agua_aire_kg': agua_B, 'E_gen_tanques_kWh': Egen_A, 'E_gen_aire_kWh': Egen_B,
    'regeneracion_solar': regen,
})
