# -*- coding: utf-8 -*-
# [EN] Operation at altitude on the design hot day (Section 9).
# =============================================================================
# exp15_altura.py — EXP-15: ¿conviene hacer la maniobra más arriba, donde el aire es más frío?
# (Pregunta del autor, 2026-09-25; comentario de la discusión del paper 1.)
#
# Perfil: el día caluroso del caso de diseño (T_aire_exterior a nivel del mar) con el gradiente
# térmico de la ISA y la presión hidrostática consistente con ese perfil:
#     T(h) = T0 - G h,   P(h) = P0 (T/T0)^(g/(R G)),   rho(h)/rho0 = (P/P0)(T0/T)
# Misma regla de diseño que en tierra: el escape sale del condensador (T_salida_cond - T_aire_exterior)
# K por encima del aire. Con eso:
#   - phi_max(T_salida, P) y el H2 por tonelada (la razón rho_aire/rho_H2 no cambia con la altura:
#     los dos gases están a la misma T y P);
#   - la masa de rechazo por kW con las dos bandas de EXP-14 (el radiador de la banda fría ve un aire
#     más frío; el kW/kg por K de FlyZero no cambia; NO se modela que el aire menos denso rinda menos);
#   - la sustentación relativa de la misma nave (densidad del aire).
# Los perfiles reales difieren del estándar (p. ej. la inversión térmica de la costa de Lima): cada
# sitio con sus datos. Esto es una cota de tendencia, no un caso de diseño.
# =============================================================================

import math
import parametros as p
import rutas
from exp05_condensacion_tiempo_real import phi_max
import exp14_calor_por_temperatura as e14

ALTURAS_M = (0.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0)
H_REF_M = 1000.0                                  # altura de comparación que cita el texto


def perfil(h):
    """Temperatura [°C], presión [Pa] y densidad relativa del aire a la altura h [m]."""
    T0 = p.T_aire_exterior + 273.15
    T = T0 - p.gradiente_ISA * h
    P = p.P_atmos * (T / T0) ** (p.g / (p.R_aire * p.gradiente_ISA))
    return T - 273.15, P, (P / p.P_atmos) * (T0 / T)


def fila(h, k_fria, k_rho):
    T_aire, P, rho_rel = perfil(h)
    acerc = p.T_salida_cond - p.T_aire_exterior
    phi = phi_max(T_aire + acerc, P)
    B_t = 1000.0 / (k_rho + (1.0 + p.ratio_O2_H2) * phi)
    b = e14.banda_fria(P, phi)
    qt = e14.q_total(0.0, phi, b['T_fin_C'], T_in_C=T_aire)
    q_cal = qt - b['Q_banda']
    kgkw = []
    for i in range(2):
        m_f = e14.masa_banda_fria(b, k_fria[i], T_aire0=T_aire)
        m_c = 1.0 / (e14.KWKG_CALIENTE[i] * 1e3)
        kgkw.append((q_cal * m_c + b['Q_banda'] * m_f) / qt * 1e3)
    return {'h_m': h, 'T_aire_C': T_aire, 'P_Pa': P, 'rho_rel': rho_rel, 'phi': phi, 'H2_por_t_kg': B_t,
            'T_rocio_C': b['T_rocio_C'], 'T_fin_C': b['T_fin_C'], 'q_MJkg': qt / 1e6,
            'q_fria_MJkg': b['Q_banda'] / 1e6, 'kg_por_kW': kgkw}


def main():
    r14 = rutas.leer_datos('exp14_calor_por_temperatura')
    k_fria = tuple(r14['k_fria_kWkg_por_K'])
    k_rho = p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1)) - 1.0      # la misma que EXP-14 y EXP-06
    filas = [fila(h, k_fria, k_rho) for h in ALTURAS_M]
    f0 = filas[0]
    ref = next(f for f in filas if f['h_m'] == H_REF_M)
    cambio = {'H2_por_t': ref['H2_por_t_kg'] / f0['H2_por_t_kg'] - 1.0,
              'kg_por_kW': [a / b - 1.0 for a, b in zip(ref['kg_por_kW'], f0['kg_por_kW'])],
              'sustentacion': ref['rho_rel'] - 1.0}
    rutas.guardar_datos('exp15_altura', {'H_ref_m': H_REF_M, 'filas': filas, 'cambio_ref': cambio})
    print('h [m]  T aire  P [atm]  rho/rho0  phi    H2/t [kg]  kg/kW (2026, 2035)')
    for f in filas:
        print(f"{f['h_m']:5.0f}  {f['T_aire_C']:6.1f}  {f['P_Pa'] / p.P_atmos:6.3f}  {f['rho_rel']:7.3f}  {f['phi']:.3f}"
              f"  {f['H2_por_t_kg']:8.2f}   {f['kg_por_kW'][0]:.4f}, {f['kg_por_kW'][1]:.4f}")
    print(f"a {H_REF_M:g} m: H2 por t {100 * cambio['H2_por_t']:+.1f} %, masa de rechazo "
          f"{100 * cambio['kg_por_kW'][0]:+.0f} / {100 * cambio['kg_por_kW'][1]:+.0f} %, sustentación "
          f"{100 * cambio['sustentacion']:+.0f} %")


if __name__ == '__main__':
    main()
