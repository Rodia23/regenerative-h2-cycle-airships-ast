# [EN] Net heaviness per kg of hydrogen and minimum condensation power (supporting step, Section 3.3).
import numpy as np
import matplotlib.pyplot as plt
import plot_style
from matplotlib.ticker import ScalarFormatter, FuncFormatter
from t2_gain_scheduled import _Params
import rutas

# 1. Official simulation parameters (todos desde parametros.py)
p = _Params()
L_v = p.L_v_agua / 1e3          # kJ/kg (Water latent heat)
t_objetivo = p.t_cond_objetivo  # Target threshold [min]

# ── VISUALIZATION PARAMETERS ──────────────────────────────────────────────────
FIGSIZE   = (7.5, 4.2)
DPI       = 150
LW_CURVA  = 2.2
LW_REF    = 1.2
FS_LABEL  = 11
FS_TITLE  = 12
FS_ANNOT  = 10
COL_CURVA = '#2c7bb6'
COL_ZONA  = 'green'
COL_REF   = 'black'      # Cambiado a negro
COL_GRID  = 'gray'
LOC_LEG   = 'upper right'
S_SCATTER = 50
# ──────────────────────────────────────────────────────────────────────────────

# 2. Water mass for critical condensation window (with O₂ absorption)
# Using L_eff = ρ_air/ρ_H2 + 8: each kg H₂ burned with atmospheric O₂
# creates L_eff × g of net downward force (buoyancy loss + O₂ ballast).
# Descent-phase H₂ requirement: m_H2 = m_load / L_eff  (O₂-corrected, §2.5)
# ≈ 44.7 kg para m_load = 1 t. Es lo que se quema en la ventana crítica de
# condensación; el resto del ciclo (agua residual) condensa pasivamente en ruta.
rho_aire = p.rho_aire_std   # densidad del aire en la altitud de referencia h0
rho_h2   = p.P_atmos / (p.R_H2 * p.T1)
L_eff    = rho_aire / rho_h2 + p.ratio_O2_H2   # effective specific lift factor
m_h2_quemada = p.m_carga / L_eff          # descent-phase H₂ (O₂-corrected)
m_h2o = m_h2_quemada * (1 + p.ratio_O2_H2)  # water produced in critical window (1+8 = 9 kg/kg)

# Minimum power threshold calculation
Q_min_kW = (m_h2o * L_v) / (t_objetivo * 60)

# 3. Data generation
Q_range = np.logspace(1.5, 4, 500) 
t_range = (m_h2o * L_v) / (Q_range * 60)

# 4. Plot rendering
plt.figure(figsize=FIGSIZE, dpi=DPI)
ax = plt.gca()

# Background and grid configuration
ax.set_facecolor('#fcfcfc')
plt.grid(True, which="both", ls="-", alpha=0.1, color=COL_GRID)

# Plot main curve and operational zone
plt.plot(Q_range, t_range, label='Condensation', color=COL_CURVA, lw=LW_CURVA, zorder=3)
plt.axhspan(1, t_objetivo, color=COL_ZONA, alpha=0.08, label='Operational zone (t < 10 min)') # Ajustado el fondo a 1 por el log

# Plot threshold markers (Negro)
plt.axvline(Q_min_kW, color=COL_REF, ls='--', lw=LW_REF, alpha=0.6)
plt.axhline(t_objetivo, color=COL_REF, ls='--', lw=LW_REF, alpha=0.6)
plt.scatter([Q_min_kW], [t_objetivo], color=COL_REF, s=S_SCATTER, edgecolors='white', zorder=5)

# Threshold annotation (Negro)
# Threshold annotation (Negro con subíndices matemáticos)
ax.annotate(f'Critical threshold:\n$Q_{{min}} \\approx {Q_min_kW:.0f}$ kW\n$t = 10$ min',
            xy=(Q_min_kW, t_objetivo), xytext=(Q_min_kW * 1.5, t_objetivo * 2.5),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color=COL_REF),
            fontsize=FS_ANNOT, fontweight='bold', color=COL_REF,
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec=COL_REF, alpha=0.9))

# Axis formatting (ESCALA LOGARÍTMICA EN AMBOS EJES)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1, 1000) # El log no puede empezar en 0, así que empieza en 1
plt.xlim(30, 10000)

# ── FORZAR NÚMEROS REALES EN LUGAR DE NOTACIÓN CIENTÍFICA ──
# Commas only for numbers with 5+ digits (MDPI style)
def _fmt(x, _):
    return f'{int(x):,}' if x >= 10000 else f'{int(x)}'
ax.xaxis.set_major_formatter(FuncFormatter(_fmt))
ax.yaxis.set_major_formatter(ScalarFormatter())

plt.xlabel('Condensation power $Q_{cond}$ [kW]', fontsize=FS_LABEL, labelpad=10)
plt.ylabel('Condensation time [min]', fontsize=FS_LABEL, labelpad=10)
# (sin título dentro de la figura: va en el caption del manuscrito)
# Agregar leyenda y ajustar los márgenes
plt.legend(frameon=True, loc=LOC_LEG, facecolor='white', framealpha=1, fontsize=9)
plt.tight_layout()

# Aquí es donde ocurre la magia de guardado
plt.savefig(rutas.figura('condensacion'), dpi=DPI, bbox_inches='tight')
rutas.mostrar()

rutas.guardar_datos('condensacion', {
    'L_eff': L_eff, 'm_carga_kg': p.m_carga, 'm_H2_descenso_kg': m_h2_quemada,
    'm_H2O_descenso_kg': m_h2o, 't_objetivo_min': t_objetivo, 'Q_cond_min_kW': Q_min_kW,
})
print(f"  L_eff = {L_eff:.4f}  m_H2 = {m_h2_quemada:.2f} kg  m_H2O = {m_h2o:.1f} kg  "
      f"Q_cond,min = {Q_min_kW:.0f} kW  (t = {t_objetivo:.0f} min)")