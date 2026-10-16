# -*- coding: utf-8 -*-
# [EN] Air admission into the hull (Section 5.1).
# =============================================================================
# req_aire_casco.py — requisito de la admisión de aire al casco (paper 1).
# El modelo supone que, cuando las celdas de gas se encogen al quemar H2, el aire
# entra al casco al instante (volumen del casco constante; Ec. 2 cuenta ese aire en
# la inercia). Aquí ese supuesto se vuelve un requisito de caja negra: el sistema de
# aire del casco debe admitir V = ṁ/rho_H2 (y expulsarlo cuando el controlador
# devuelve H2 desde los tanques, a lo sumo a la misma tasa instalada).
#   - cota ideal: ṁ* de la liberación de 1 t con los umbrales de 1 y 5 m;
#   - escalones que caben con la potencia de rechazo instalada (EXP-06, 10 MW, 5 m);
#   - comparación: el aire de combustión que toma la unidad de conversión.
# No hay datos publicados de sopladores de balonet: se da el caudal, no un equipo.
# =============================================================================

import parametros as p
import rutas
import minimax as mm
from exp05_condensacion_tiempo_real import AIRE_ESTEQ

P_ESCALON_MW = 10.0


def main():
    rho_h = p.P_atmos / (p.R_H2 * p.T1)
    pasos = {x['corto']: x['paso_max_kg'] for x in rutas.leer_datos('exp06_requisitos_minimos')['pasos']
             if x['h_m'] == p.h_max_seguridad and x['P_instalada_MW'] == P_ESCALON_MW}
    filas = []
    for n, d in p.PLATAFORMAS.items():
        r1 = mm.resolver(n, 1000.0, p.h_max_precision, p.D_lim_oper)
        r5 = mm.resolver(n, 1000.0, p.h_max_seguridad, p.D_lim_oper)
        m_esc = pasos[d['corto']]
        re = mm.resolver(n, m_esc, p.h_max_seguridad, p.D_lim_oper)
        V_casco = d['V_casco']
        filas.append({'plataforma': n, 'corto': d['corto'],
                      'mdot_1m': r1['dm_star'], 'aire_1m_m3s': r1['dm_star'] / rho_h,
                      'mdot_5m': r5['dm_star'], 'aire_5m_m3s': r5['dm_star'] / rho_h,
                      'escalon_kg': m_esc, 'mdot_escalon': re['dm_star'], 'aire_escalon_m3s': re['dm_star'] / rho_h,
                      'aire_comb_1m_m3s': r1['dm_star'] * AIRE_ESTEQ / p.rho_aire_std,
                      'V_casco_m3': V_casco, 'frac_casco_por_s_1m': r1['dm_star'] / rho_h / V_casco})
    B_t = 1000.0 / mm.L_eff
    rutas.guardar_datos('req_aire_casco', {'rho_H2': rho_h, 'P_escalon_MW': P_ESCALON_MW,
                                           'volumen_por_t_m3': B_t / rho_h, 'filas': filas,
                                           'verif_sim': verificar_con_simulacion(rho_h)})
    print(f'aire que entra al casco por tonelada liberada: {B_t / rho_h:.0f} m3')
    for f in filas:
        print(f"  {f['corto']:<5} 1 m: {f['aire_1m_m3s']:6.1f} m3/s  5 m: {f['aire_5m_m3s']:6.1f}  "
              f"escalón {f['escalon_kg']:.0f} kg: {f['aire_escalon_m3s']:5.2f} m3/s  "
              f"(aire de combustión {f['aire_comb_1m_m3s']:.0f} m3/s; {100*f['frac_casco_por_s_1m']:.3f} % del casco por s)")


def verificar_con_simulacion(rho_h):
    """Caso base simulado (1 t, 1 m, condensación instantánea): máximo |dV_gas/dt| frente a ṁ*/rho_H2."""
    import numpy as np
    import drop_flujo_optimo_perdrop as sim
    sim.usar_nave(p.nombre_dirigible, 1000.0)
    r = sim.correr(n_entregas=1, toneladas_por_drop=1.0, graficar=False, verbose=False, devolver_serie=True)
    s = r['serie']
    V = np.asarray(s['m_H2']) / rho_h
    dVdt = np.gradient(V, np.asarray(s['t']))
    return {'max_dVdt_sim_m3s': float(np.max(np.abs(dVdt))), 'dm_opt_sim': r['dm_opt'],
            'aire_formula_m3s': r['dm_opt'] / rho_h}


if __name__ == '__main__':
    main()
