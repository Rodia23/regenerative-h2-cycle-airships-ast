# -*- coding: utf-8 -*-
# [EN] Comparison with venting and with gas compression, hydrogen and helium (Section 7, Table 6).
# =============================================================================
# exp12_helio.py — EXP-12: comparación con helio (punto 3 del checklist del paper 1;
# Paso 4 del plan, sección 3.4), a nivel de CAJA NEGRA, igual que el H2.
#
# Mismo mecanismo de extracción: sacar gas de sustentación de la envolvente y dejar
# que el balonet admita aire. Lo que cambia es qué pasa con el gas:
#   H2 quemado, agua retenida   L = rho_a/rho_H2 - 1 + r phi  (r = M_H2O/M_H2)
#   H2 venteado                 L = rho_a/rho_H2 - 1
#   He comprimido a bordo       L = rho_a/rho_He   (COSH: el He queda en tanques)
#   He venteado                 L = rho_a/rho_He - 1
# (L = kg de pesadez ganados por kg de gas procesado.)
#
# Requisitos de cada caja negra, sin elegir equipo:
#   H2: calor a rechazar q (balance de EXP-05) -> unidad de rechazo de calor.
#   He: trabajo mínimo de compresión isotérmica reversible hasta la presión de
#       almacenamiento, W = (h2 - h1) - T (s2 - s1) con propiedades de gas real de
#       NIST; ese trabajo también sale como calor -> unidad de compresión + rechazo.
# Se evalúan: (1) por kg de gas y por tonelada liberada; (2) el caso base del
# artículo (PF, 1 t, h_max = 1 m, D_lim = 0.10 m) con la misma forma cerrada minimax
# (minimax.resolver con rho_gas y L_gas); (3) liberación por escalones con la misma
# potencia instalada (tiempo por tonelada = energía por tonelada / potencia).
# La masa de los tanques de He y del condensador del H2 NO se evalúa (depende del
# equipo): queda como brecha de cada ruta.
# Salidas: resultados/datos/exp12_helio.json, resultados/tablas/exp12_helio.tex
# =============================================================================

import io
import math
import sys
import parametros as p
import rutas
import minimax as mm
from exp05_condensacion_tiempo_real import q_rech, phi_max

# NIST Chemistry WebBook (SRD 69; He: Ortiz-Vega et al.), isoterma de 15 °C, gas real; consulta 2026-09-24.
# P [MPa] -> (densidad [kg/m3], entalpía [kJ/kg], entropía [kJ/kg K])
NIST_HE_15C = {0.101325: (0.1691952, 1501.784, 27.77087),
               5.0: (8.15230, 1517.71, 19.6769),
               20.0: (30.4595, 1565.77, 16.8116),
               30.0: (43.8120, 1597.28, 15.9785)}
# Auditoría 2 (2026-10-06, M1): H2 comprimido a bordo (COSH con el propio H2 de sustentación). Isoterma de 15 °C de
# normal-hidrógeno con la ecuación de estado de Leachman et al. (2009), la de NIST WebBook SRD 69 (referencia [50] del
# artículo), evaluada con CoolProp 8.0.0 (NIST no accesible desde el contenedor). Sus densidades a 1.5, 35 y 70 MPa
# reproducen las de parametros.NIST_RHO_H2_15C (consulta a NIST del 2026-09-29): lo comprueba verificacion.py.
# P [MPa] -> (densidad [kg/m3], entalpía [kJ/kg], entropía [kJ/kg K])
NIST_H2_15C = {0.101325: (0.08520494, 3788.9325, 52.888566),
               1.5: (1.2508613, 3794.5972, 41.758621),
               35.0: (23.994753, 3980.5547, 28.545058),
               70.0: (40.172161, 4231.4588, 25.604522)}
P_H2_COMP_MPA = 35.0                          # los tanques del ciclo (parametros.P_tanque_bar)
P_ALMACEN_MPA = (5.0, 20.0, 30.0)
POTENCIAS_MW = (10.0, 100.0)                  # potencia instalada de la caja negra (como EXP-06)
PLATAFORMA = p.nombre_dirigible               # caso base del artículo


def w_iso_h2(P_MPa):
    """Trabajo mínimo [J/kg] para comprimir H2 isotérmica y reversiblemente de 1 atm a P (gas real, Leachman)."""
    T = p.T1
    r1, h1, s1 = NIST_H2_15C[0.101325]
    r2, h2, s2 = NIST_H2_15C[P_MPa]
    return ((h2 - h1) - T * (s2 - s1)) * 1e3


def w_iso_he(P_MPa):
    """Trabajo mínimo [J/kg] para comprimir He isotérmica y reversiblemente de 1 atm a P (gas real, NIST)."""
    T = p.T1
    r1, h1, s1 = NIST_HE_15C[0.101325]
    r2, h2, s2 = NIST_HE_15C[P_MPa]
    return ((h2 - h1) - T * (s2 - s1)) * 1e3


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    rho_a = p.rho_aire_std
    rho_h2 = p.P_atmos / (p.R_H2 * p.T1)                        # gas ideal, como el resto del modelo
    rho_he = p.P_atmos / (p.R_univ / p.M_He * p.T1)
    r = 1.0 + p.ratio_O2_H2
    phi_1atm = phi_max(p.T_salida_cond)
    phi_obj = p.phi_objetivo

    # ── (1) Por kg de gas y por tonelada liberada ────────────────────────────
    W = {P: w_iso_he(P) for P in P_ALMACEN_MPA}
    W_ideal = {P: p.R_univ / p.M_He * p.T1 * math.log(P * 1e6 / p.P_atmos) for P in P_ALMACEN_MPA}
    rutas_gas = [
        {'ruta': 'H2 quemado, agua retenida (phi = 0.99)', 'gas': 'H2', 'L': rho_a / rho_h2 - 1 + r * phi_obj,
         'rho': rho_h2, 'q_Jkg': q_rech(0.0, phi_obj), 'perdida_frac': 1.0 - phi_obj},
        {'ruta': f'H2 quemado, agua retenida a 1 atm (phi = {phi_1atm:.2f})', 'gas': 'H2',
         'L': rho_a / rho_h2 - 1 + r * phi_1atm, 'rho': rho_h2, 'q_Jkg': q_rech(0.0, phi_1atm),
         'perdida_frac': 1.0 - phi_1atm},
        {'ruta': 'H2 venteado', 'gas': 'H2', 'L': rho_a / rho_h2 - 1, 'rho': rho_h2, 'q_Jkg': 0.0, 'perdida_frac': 1.0},
        {'ruta': 'He comprimido a bordo (COSH)', 'gas': 'He', 'L': rho_a / rho_he, 'rho': rho_he,
         'q_Jkg': W[20.0], 'perdida_frac': 0.0},
        {'ruta': 'He venteado', 'gas': 'He', 'L': rho_a / rho_he - 1, 'rho': rho_he, 'q_Jkg': 0.0, 'perdida_frac': 1.0},
    ]
    for g in rutas_gas:
        g['kg_por_t'] = 1000.0 / g['L']
        g['m3_por_t'] = g['kg_por_t'] / g['rho']                 # volumen de succión a 1 atm
        g['calor_GJ_por_t'] = g['kg_por_t'] * g['q_Jkg'] / 1e9
        g['perdido_kg_por_t'] = g['kg_por_t'] * g['perdida_frac']
    he_t = 1000.0 / (rho_a / rho_he)
    almacen = [{'P_MPa': P, 'W_iso_MJkg': W[P] / 1e6, 'W_iso_ideal_MJkg': W_ideal[P] / 1e6,
                'calor_GJ_por_t': he_t * W[P] / 1e9, 'm3_almacen_por_t': he_t / NIST_HE_15C[P][0]}
               for P in P_ALMACEN_MPA]
    # revisión 2026-10-06 (comparación en dos etapas): masa del tanque de COSH, SOLO el compuesto, con el escalado de
    # pared delgada del Tipo IV del DOE (compuesto proporcional a P·V, parametros.tanque_a_presion): cota inferior
    d24 = p.TANQUE_DOE24
    kg_por_PV = d24['compuesto_kg'] / (d24['P_bar'] * 1e5 * d24['V_L'] / 1e3)       # [kg/J]
    for a in almacen:
        a['tanque_compuesto_t_por_t'] = kg_por_PV * a['P_MPa'] * 1e6 * a['m3_almacen_por_t'] / 1e3
    for g in rutas_gas:                                             # agua retenida por tonelada (las que queman)
        quema = g['gas'] == 'H2' and g['q_Jkg'] > 0.0                  # r phi B de agua, densidad 1000 kg/m3
        g['agua_m3_por_t'] = g['kg_por_t'] * (1.0 + p.ratio_O2_H2) * (1.0 - g['perdida_frac']) / 1000.0 if quema else 0.0
    # revisión 2026-10-06 (corrección del autor): para devolver la sustentación, las dos rutas sacan gas de sus tanques
    # a bordo. Cada kg sacado devuelve rho_aire/rho_gas, porque su masa ya está a bordo. El ciclo saca los B kg que
    # quemó: devuelve B rho_a/rho_H2 al instante, y el resto, B (r phi - 1) (el O2 tomado del aire, en el agua, menos el
    # H2 perdido), vuelve cuando la electrólisis parte el agua: el H2 rellena los tanques y el O2 vuelve al aire. COSH descomprime todo su gas. Compuesto de esos
    # tanques con el mismo escalado P·V: el ciclo a P_tanque_bar (35 MPa, NIST), COSH a 20 MPa. Venteos: sin tanques.
    rho_h2_t = p.NIST_RHO_H2_15C[p.P_tanque_bar / 10.0]
    alm_ref = next(a for a in almacen if a['P_MPa'] == 20.0)
    for g in rutas_gas:
        if g['gas'] == 'H2' and g['q_Jkg'] > 0.0:
            g['reserva_kg'] = g['kg_por_t']
            g['reserva_P_MPa'] = p.P_tanque_bar / 10.0
            g['reserva_compuesto_t'] = kg_por_PV * p.P_tanque_bar * 1e5 * g['kg_por_t'] / rho_h2_t / 1e3
        elif g['gas'] == 'He' and g['q_Jkg'] > 0.0:
            g['reserva_kg'] = g['kg_por_t']
            g['reserva_P_MPa'] = 20.0
            g['reserva_compuesto_t'] = alm_ref['tanque_compuesto_t_por_t']
        else:
            g['reserva_kg'] = None
            continue
        g['al_instante_kg'] = g['reserva_kg'] * rho_a / g['rho']
        g['resto_kg'] = max(0.0, 1000.0 - g['al_instante_kg'])
    # devolver la tonelada entera de golpe con H2: mismos moles que el He (gas ideal: n = 1000 / M_aire)
    todo_h2 = 1000.0 / (rho_a / rho_h2)
    reserva_todo = {'H2_kg': todo_h2, 'H2_compuesto_t': kg_por_PV * p.P_tanque_bar * 1e5 * todo_h2 / rho_h2_t / 1e3,
                    'H2_kmol': todo_h2 / p.M_H2 / 1e3, 'He_kmol': he_t / p.M_He / 1e3, 'aire_kmol': 1000.0 / p.M_aire / 1e3,
                    'rho_H2_tanque': rho_h2_t, 'kg_por_PV': kg_por_PV}

    # ── (2) Caso base del artículo con la misma forma cerrada ────────────────
    d = p.PLATAFORMAS[PLATAFORMA]
    base_h2 = mm.resolver(d, p.m_carga, p.h_max_precision, p.D_lim_oper)
    base_he = mm.resolver(d, p.m_carga, p.h_max_precision, p.D_lim_oper,
                          rho_gas=rho_he, L_gas=rho_a / rho_he)
    q99 = q_rech(0.0, phi_obj)
    base = {
        'H2': {'dm_kgps': base_h2['dm_star'], 'alpha': base_h2['alpha'], 'M_total_kg': base_h2['M_total'],
               'succion_m3ps': base_h2['dm_star'] / rho_h2, 'potencia_MW': base_h2['dm_star'] * q99 / 1e6,
               'que': 'calor a rechazar (condensar phi = 0.99)'},
        'He': {'dm_kgps': base_he['dm_star'], 'alpha': base_he['alpha'], 'M_total_kg': base_he['M_total'],
               'succion_m3ps': base_he['dm_star'] / rho_he,
               'potencia_MW': {P: base_he['dm_star'] * W[P] / 1e6 for P in P_ALMACEN_MPA},
               'que': 'trabajo mínimo de compresión (sale como calor)'},
    }

    # ── (3) Escalones con la misma potencia instalada ────────────────────────
    h2_t = 1000.0 / rutas_gas[0]['L']
    escalones = [{'P_MW': P, 'H2_min_por_t': h2_t * q99 / (P * 1e6) / 60.0,
                  'He_min_por_t': {Pa: he_t * W[Pa] / (P * 1e6) / 60.0 for Pa in P_ALMACEN_MPA}}
                 for P in POTENCIAS_MW]

    # Auditoría 2 (M1): H2 comprimido a bordo, a la presión de los tanques del ciclo, con la misma forma cerrada
    assert abs(P_H2_COMP_MPA - p.P_tanque_bar / 10.0) < 1e-9
    L_h2c = rho_a / rho_h2
    W_h2c = w_iso_h2(P_H2_COMP_MPA)
    base_h2c = mm.resolver(d, p.m_carga, p.h_max_precision, p.D_lim_oper, rho_gas=rho_h2, L_gas=L_h2c)
    kg_h2c = 1000.0 / L_h2c
    h2_comprimido = {'P_MPa': P_H2_COMP_MPA, 'L': L_h2c, 'kg_por_t': kg_h2c, 'W_iso_MJkg': W_h2c / 1e6,
                     'W_iso_ideal_MJkg': p.R_H2 * p.T1 * math.log(P_H2_COMP_MPA * 1e6 / p.P_atmos) / 1e6,
                     'calor_GJ_por_t': kg_h2c * W_h2c / 1e9, 'dm_kgps': base_h2c['dm_star'],
                     'potencia_MW': base_h2c['dm_star'] * W_h2c / 1e6,
                     'm3_almacen_por_t': kg_h2c / NIST_H2_15C[P_H2_COMP_MPA][0],
                     'reserva_compuesto_t': kg_por_PV * P_H2_COMP_MPA * 1e6 * kg_h2c / NIST_H2_15C[P_H2_COMP_MPA][0] / 1e3,
                     'al_instante_kg': kg_h2c * rho_a / rho_h2, 'perdido_kg_por_t': 0.0}
    # Auditoría 2 (M3): el resto de la sustentación vuelve al ritmo del electrolizador de la tabla 6, dimensionado para
    # regenerar en t_regeneracion_h el H2 de una tonelada a 1 atm (balance_masa: la ruta 1 de aquí)
    g1 = rutas_gas[1]
    ritmo_el = g1['kg_por_t'] * (1.0 - g1['perdida_frac']) / p.t_regeneracion_h          # [kg/h]
    for g in rutas_gas[:2]:
        g['t_resto_h'] = g['kg_por_t'] * (1.0 - g['perdida_frac']) / ritmo_el

    salida = {'rho_He': rho_he, 'rho_H2': rho_h2, 'rho_He_sobre_rho_H2': rho_he / rho_h2,
              'rho_He_NIST': NIST_HE_15C[0.101325][0], 'rutas': rutas_gas, 'almacen': almacen,
              'base': base, 'plataforma_base': PLATAFORMA, 'escalones': escalones,
              'factor_flujo_base': base_he['dm_star'] / base_h2['dm_star'],
              'factor_flujo_L': rutas_gas[0]['L'] / (rho_a / rho_he),
              'factor_calor_20MPa': rutas_gas[0]['calor_GJ_por_t'] / (he_t * W[20.0] / 1e9),
              'P_ref_MPa': 20.0, 'reserva_todo': reserva_todo, 'h2_comprimido': h2_comprimido,
              'electrolizador_kgph': ritmo_el}
    rutas.guardar_datos('exp12_helio', salida)
    tabla(rutas_gas)

    print(f'rho_He/rho_H2 = {rho_he / rho_h2:.3f} (gas ideal; NIST real: '
          f'{NIST_HE_15C[0.101325][0] / 0.08520506:.3f})\n')
    print(f'{"ruta":<46}{"L":>7}{"kg/t":>8}{"m3/t":>8}{"calor GJ/t":>12}{"perdido kg/t":>14}')
    for g in rutas_gas:
        print(f'{g["ruta"]:<46}{g["L"]:7.2f}{g["kg_por_t"]:8.1f}{g["m3_por_t"]:8.0f}'
              f'{g["calor_GJ_por_t"]:12.2f}{g["perdido_kg_por_t"]:14.1f}')
    print('\nHe comprimido: trabajo mínimo isotérmico (NIST, gas real) y almacenamiento por tonelada')
    for a in almacen:
        print(f'  {a["P_MPa"]:5.0f} MPa: W = {a["W_iso_MJkg"]:.2f} MJ/kg (ideal {a["W_iso_ideal_MJkg"]:.2f}), '
              f'{a["calor_GJ_por_t"]:.2f} GJ/t, {a["m3_almacen_por_t"]:.2f} m3 de He comprimido por t')
    print(f'\nCaso base {PLATAFORMA} 1 t, h_max = {p.h_max_precision:g} m, D_lim = {p.D_lim_oper:g} m (misma forma cerrada):')
    b2, be = base['H2'], base['He']
    print(f'  H2: {b2["dm_kgps"]:.2f} kg/s (alpha {b2["alpha"]:.3f}), succión {b2["succion_m3ps"]:.0f} m3/s, '
          f'{b2["potencia_MW"]:.0f} MW de calor')
    print(f'  He: {be["dm_kgps"]:.2f} kg/s (alpha {be["alpha"]:.3f}), succión {be["succion_m3ps"]:.0f} m3/s, '
          + ', '.join(f'{v:.0f} MW a {k:g} MPa' for k, v in be['potencia_MW'].items()))
    print(f'  factor de flujo He/H2: {salida["factor_flujo_base"]:.2f} (regla de tres con L: {salida["factor_flujo_L"]:.2f})')
    print('\nEscalones con la misma potencia instalada (minutos por tonelada):')
    for e in escalones:
        print(f'  {e["P_MW"]:5.0f} MW: H2 {e["H2_min_por_t"]:.1f} min/t; He '
              + ', '.join(f'{v:.2f} min/t a {k:g} MPa' for k, v in e['He_min_por_t'].items()))


def tabla(rutas_gas):
    L = [r'\begin{tabular}{l r r r r r}', r'\toprule',
         r'Ruta & $L$ & Gas & Succión & Calor & Perdido \\',
         r' & [kg/kg] & [kg/t] & [m$^3$/t] & [GJ/t] & [kg/t] \\', r'\midrule']
    for g in rutas_gas:
        nom = (g['ruta'].replace('H2', r'H$_2$').replace('He ', r'He ').replace('phi', r'$\phi$'))
        L.append(f"{nom} & {g['L']:.2f} & {g['kg_por_t']:.1f} & {g['m3_por_t']:.0f} & "
                 f"{g['calor_GJ_por_t']:.2f} & {g['perdido_kg_por_t']:.1f}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('exp12_helio'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/exp12_helio.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
