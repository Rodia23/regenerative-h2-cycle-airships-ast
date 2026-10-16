# -*- coding: utf-8 -*-
# [EN] Mass budget of the black boxes: per kg/s of hydrogen and per tonne released (Section 5.1, Table S7).
# =============================================================================
# balance_masa.py — balance de masa por caja negra (paper 1), con una clase de
# tecnología representativa por caja y sin elegir equipo. Dos referencias:
#   - base de la TABLA 6 del paper (la escribe req_masa_rechazo.py, por kg/s de H2, desde 2026-09-29): el flujo
#     IDEAL del caso base (parametros.nombre_dirigible, 1 t de una vez, ṁ* de minimax.resolver, phi objetivo):
#     la masa que pediría la cota ideal (autor, 2026-09-27);
#   - paquete de P_ref_masa (10 MW de rechazo), que usan las derivaciones y EXP-16.
#   conversión:   motor alternativo de H2 a lambda 1 (Rotax 916 iS como proxy de aviación a DEUTZ
#                 TCG 7.8 H2), para la potencia de eje eta PCI ṁ (hallazgo I2; antes turbina CT7)
#   rechazo:      sistema de gestión térmica completo, FlyZero p. 4 (caja negra; antes dos bandas de EXP-14,
#                 reemplazadas por los hallazgos I5/S1 de la auditoría, 2026-09-27)
#   buffer:       batería a nivel paquete, dimensionada por potencia (EXP-11)
#   regeneración: electrolizador PEM que recupera en 24 h el H2 de una tonnelada
#   tanques:      Tipo IV a 350 bar para TODO el H2 de una tonelada (EXP-11, retardo 5 s): el regenerado,
#                 que se guarda comprimido para no cambiar la flotabilidad, más la reposición (autor, 2026-09-27;
#                 antes solo la reposición)
#   paneles:      opcionales, dos tecnologías con su eficiencia (a-Si y CIGS); días por tonelada en
#                 invierno, verano y anual con la irradiancia sobre el casco (irradiancia_casco.py)
# Además: qué fracción de la electricidad de la electrólisis cubre la propia quema.
# =============================================================================

import parametros as p
import rutas
import minimax
from exp05_condensacion_tiempo_real import q_rech, phi_max, rho_liq_sat
from exp02_quema_centrada import L_ef

# Las cajas que van por tonelada liberada (fijas); el resto escala con el flujo. El tanque de agua (D, 2026-10-13)
# guarda el agua retenida de la tonelada: en una suelta única va con ella; en el reparto, con todo lo entregado
FIJAS = ('regeneracion', 'tanques_por_t', 'tanque_agua', 'compresor')


def tanque_agua():
    """Tanque de agua de la clase MIL-T-53029C (parametros.TANQUE_AGUA): tela de las dos caras de la almohada vacía
    por kg de agua a su capacidad, con el agua a la salida del condensador; y el agua retenida por tonelada liberada
    con la phi objetivo (la del ciclo, 14 atm), que es lo que guarda."""
    ta = p.TANQUE_AGUA
    area = ta['caras'] * (ta['lado_vacio_ft'] * p.FT) ** 2                         # [m2]
    agua_cap = ta['capacidad_gal'] * p.GAL * rho_liq_sat(p.T_salida_cond)       # [kg] agua a su capacidad
    kg_por_kg = tuple(area * w * p.OZ / p.YD2 / agua_cap for w in ta['tela_ozyd2'])
    agua_t = (1.0 + p.ratio_O2_H2) * p.phi_objetivo * 1000.0 / L_ef(p.phi_objetivo)   # [kg] por tonelada liberada
    return {'area_m2': area, 'agua_cap_kg': agua_cap, 'kg_por_kg_agua': kg_por_kg, 'agua_por_t_kg': agua_t,
            'kg_por_t': tuple(k * agua_t for k in kg_por_kg)}



def f_caida_presion():
    """Fracción de ventilación por caída de presión (R-03, revisión 2026-10-01), una por cada dp de
    parametros.dp_rechazo_Pa: f = dp / (eta_vent rho_aire cp dT_aire), con el aire exterior."""
    rho_ext = p.P_atmos / (p.R_aire * (p.T_aire_exterior + 273.15))
    return tuple(dp / (p.eta_ventilador * rho_ext * p.cp_aire * p.dT_aire_enfr) for dp in p.dp_rechazo_Pa)

def masa_buffer(P_buf, t_buf, kWkg):
    """Batería a nivel paquete que da P_buf [W] durante t_buf [s]: manda la potencia o la energía (regla de EXP-11)."""
    return max(P_buf * t_buf / 3600.0 / p.Whkg_bateria, P_buf / 1e3 / kWkg)


def main():
    P = p.P_ref_masa
    phi = phi_max(p.T_salida_cond)
    eta = p.rend_gen
    q = q_rech(eta, phi)
    mdot = P / q                                        # quema que produce P de calor a rechazar [kg/s]
    P_eje = eta * p.Q_comb_H2 * mdot
    kwkg_turb = [w / 1e3 / m for w, m in p.TURBOEJE_TCDS]          # solo para la nota de la turbina
    kwkg_motor = [w / 1e3 / m for w, m in p.MOTORES_ALTERNATIVOS]
    kgkw_rech = [1.0 / k for k in p.kWkg_TMS_sistema]     # [kg/kW] 2026, 2035
    e11 = rutas.leer_datos('exp11_balance_energia')
    buf = [b['m_bat_kg'] for b in e11['buffer'] if b['P_instalada_MW'] == P / 1e6]
    k_rho = p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1)) - 1.0
    B_t = 1000.0 / (k_rho + (1.0 + p.ratio_O2_H2) * phi)       # H2 por tonelada (1 atm, 40 °C)
    recup_t = B_t * phi                                         # H2 que vuelve por tonelada (agua retenida / r)
    P_pem = recup_t * p.E_H2_electrolisis * 1e3 / 24.0          # [W] para recuperar una tonelada en 24 h
    # H2 que los tanques deben poder guardar por tonelada = todo el quemado (regenerado + reposición)
    todo_t = [g['H2_por_t_kg'] for g in e11['regeneracion'] if g['caso'] == 'ventear_tau5s' and g['corto'] != 'ZNT']
    tanque_por_kg = (1.0 - p.wt_COPV) / p.wt_COPV

    cajas = [
        ('conversion', P_eje / 1e3 / max(kwkg_motor), P_eje / 1e3 / min(kwkg_motor)),
        ('rechazo', P / 1e3 * min(kgkw_rech), P / 1e3 * max(kgkw_rech)),
        ('buffer', min(buf), max(buf)),
        ('regeneracion', P_pem / 1e3 / max(p.kWkg_PEM), P_pem / 1e3 / min(p.kWkg_PEM)),
        ('tanques_por_t', min(todo_t) * tanque_por_kg, max(todo_t) * tanque_por_kg),
        # compresor de todo el H2 regenerado hacia los tanques, al ritmo del electrolizador (arquitectura (b)
        # del autor, EXP-16): equipo eléctrico de catálogo del tamaño pedido, cota superior (fuentes.HPAC)
        ('compresor', p.COMPRESOR_REF['masa_kg'][0], p.COMPRESOR_REF['masa_kg'][1]),
        # el agua retenida de la tonelada, en un tanque de tela colapsable (fuentes.MILT53029)
        ('tanque_agua',) + tanque_agua()['kg_por_t'],
    ]
    total = (sum(c[1] for c in cajas), sum(c[2] for c in cajas))

    # ── Flujo IDEAL del caso base (tabla 6): 1 t de una vez en el caso base ──
    pf = p.NAVE_BASE
    rr = minimax.resolver(pf, 1000.0, p.h_max_precision, p.D_lim_oper)
    phi_i = p.phi_objetivo
    mdot_i = rr['dm_star']
    q_i = q_rech(eta, phi_i)
    P_rej_i = mdot_i * q_i                                  # [W] calor a rechazar con la conversión a eta
    P_eje_i = eta * p.Q_comb_H2 * mdot_i
    t_buf_i = p.n_tau_arranque * p.tau_condensador + rr['t_anticip'] + rr['t_reactiva']
    buf_i = [masa_buffer(fr * P_rej_i, t_buf_i, kw) for fr, kw in zip(p.perdida_rad, p.kWkg_bateria)]
    cajas_i = [
        ('conversion', P_eje_i / 1e3 / max(kwkg_motor), P_eje_i / 1e3 / min(kwkg_motor)),
        ('rechazo', P_rej_i / 1e3 * min(kgkw_rech), P_rej_i / 1e3 * max(kgkw_rech)),
        ('buffer', min(buf_i), max(buf_i)),
    ] + [c for c in cajas if c[0] in FIJAS]   # van por tonelada
    total_i = (sum(c[1] for c in cajas_i), sum(c[2] for c in cajas_i))
    turbina_i = (P_eje_i / 1e3 / max(kwkg_turb), P_eje_i / 1e3 / min(kwkg_turb))
    # el motor es una ELECCIÓN (el balance solo pide eta >= eta_min): su eje sobra durante la quema frente a los
    # ventiladores. Arreglo alternativo: quemador (eta = 0, sin trabajo) y el buffer mueve los ventiladores; la
    # masa del quemador no se cuenta (auditoría R2-02, autor 2026-09-27)
    eje_sobra = (P_eje_i - max(p.perdida_rad) * P_rej_i, P_eje_i - min(p.perdida_rad) * P_rej_i)
    P_rej_q = mdot_i * q_rech(0.0, phi_i)
    buf_q = [masa_buffer(fr * P_rej_q, t_buf_i, kw) for fr, kw in zip(p.perdida_rad, p.kWkg_bateria)]
    cajas_q = [
        ('rechazo', P_rej_q / 1e3 * min(kgkw_rech), P_rej_q / 1e3 * max(kgkw_rech)),
        ('buffer', min(buf_q), max(buf_q)),
    ] + [c for c in cajas if c[0] in FIJAS]
    total_q = (sum(c[1] for c in cajas_q), sum(c[2] for c in cajas_q))
    quemador = {'P_rej_W': P_rej_q, 'cajas': [{'caja': c[0], 'min_kg': c[1], 'max_kg': c[2]} for c in cajas_q],
                'total_kg': total_q, 'veces_carga': (total_q[0] / pf['carga_max'], total_q[1] / pf['carga_max'])}
    ideal = {'plataforma': p.nombre_dirigible, 'm_kg': 1000.0, 'phi': phi_i, 'mdot_kgps': mdot_i, 'q_Jkg': q_i,
             'P_rej_W': P_rej_i, 'P_rej_eta0_W': mdot_i * q_rech(0.0, phi_i), 'P_eje_W': P_eje_i, 't_buf_s': t_buf_i,
             'cajas': [{'caja': c[0], 'min_kg': c[1], 'max_kg': c[2]} for c in cajas_i], 'total_kg': total_i,
             'veces_carga': (total_i[0] / pf['carga_max'], total_i[1] / pf['carga_max']),
             'turbina_kg': turbina_i, 'kWkg_motor': kwkg_motor, 'eje_sobra_W': eje_sobra, 'quemador': quemador}

    # Flujo que cabe en la carga útil (autor, 2026-09-28: la brecha de masa se da como flujo alcanzable, no en
    # toneladas). Regeneración, tanques y compresor van por tonelada (fijos); el resto escala con el flujo, así que
    # con toda la carga útil para el equipo cabe la fracción x = (carga - fijos) / (lo que escala) del flujo ideal.
    # Extremos: quemador con el extremo liviano de cada caja; motor con el extremo pesado.
    cabe = {}
    for nom, cajas_x, lado in (('quemador_liviano', cajas_q, 1), ('motor_pesado', cajas_i, 2)):
        fijo = sum(c[lado] for c in cajas_x if c[0] in FIJAS)
        escala = sum(c[lado] for c in cajas_x if c[0] not in FIJAS)
        # (2026-10-14, objeción 1.8) la tonelada que se suelta va a bordo hasta la suelta: sale del presupuesto
        x = (pf['carga_max'] - fijo - ideal['m_kg']) / escala
        cabe[nom] = {'fijo_kg': fijo, 'escala_kg': escala, 'fraccion': x, 'mdot_kgps': x * mdot_i,
                     'veces_menos': 1.0 / x, 'escalon_kg': 1000.0 * x ** (2.0 / 3.0)}   # escalón: ley m^{3/2}
    ideal['cabe_en_carga'] = cabe

    # electricidad de la electrólisis: cuánto cubre la propia quema (1 atm). Revisión 2026-10-01 (R-03): los
    # ventiladores con la f por caída de presión; la f de radiadores de avión (EXP-11) queda como cota alta
    aa = [a for a in e11['autoabasto'] if a['caso'] == '1atm']
    E_el_t = recup_t * p.E_H2_electrolisis                       # kWh por tonelada
    sobra_av = [a['E_gen_kWh_t'] - a['E_pide_kWh_t'] for a in aa]
    cubre_avion = [s / E_el_t for s in sobra_av]
    Q_t = aa[0]['E_pide_kWh_t'] / aa[0]['f']                     # calor por tonelada a 1 atm (sin compresión)
    sobra = [aa[0]['E_gen_kWh_t'] - f * Q_t for f in f_caida_presion()]
    cubre = [s / E_el_t for s in sobra]

    # irradiancia sobre los paneles del casco, medida en Lima (SENAMHI) x factor de forma, por estación
    G_ef = {e: x['efectiva'] for e, x in rutas.leer_datos('irradiancia_casco')['estaciones'].items()}
    plataformas = []
    for n, d in p.PLATAFORMAS.items():
        A_pv = p.A_paneles_fraccion * d['A_total']
        pv = {}
        for tec, (ef, kgm2) in p.PV_TEC.items():
            pv[tec] = {'masa_kg': A_pv * kgm2,
                       'dias_por_t': {e: recup_t / (g * ef * A_pv / p.E_H2_electrolisis) for e, g in G_ef.items()}}
        plataformas.append({'plataforma': n, 'corto': d['corto'], 'carga_max_kg': d['carga_max'],
                            'frac_carga': (total[0] / d['carga_max'], total[1] / d['carga_max']),
                            'A_pv_m2': A_pv, 'pv': pv})

    rutas.guardar_datos('balance_masa', {
        'P_ref_W': P, 'eta': eta, 'phi': phi, 'q_Jkg': q, 'mdot_kgps': mdot, 'P_eje_W': P_eje,
        'kWkg_turbina': kwkg_turb, 'kgkW_rechazo': kgkw_rech, 'P_pem_W': P_pem, 'tanque_por_kg_H2': tanque_por_kg,
        'recuperado_por_t_kg': recup_t, 'E_el_por_t_kWh': E_el_t, 'sobra_por_t_kWh': sobra, 'cubre': cubre,
        'cubre_avion': cubre_avion,
        'cajas': [{'caja': c[0], 'min_kg': c[1], 'max_kg': c[2]} for c in cajas], 'fijas': FIJAS,
        'tanque_agua': tanque_agua(),
        'total_kg': total, 'plataformas': plataformas, 'kWkg_motor': kwkg_motor, 'ideal': ideal})
    print(f'IDEAL {ideal["plataforma"]} 1 t: quema {mdot_i:.3f} kg/s, eje {P_eje_i/1e6:.0f} MW, rechazo {P_rej_i/1e6:.0f} MW '
          f'({ideal["P_rej_eta0_W"]/1e6:.0f} MW sin exportar trabajo), buffer {t_buf_i:.1f} s')
    for c in cajas_i:
        print(f'  {c[0]:<14}{c[1]/1e3:8.1f} -{c[2]/1e3:8.1f} t')
    print(f'  TOTAL        {total_i[0]/1e3:8.1f} -{total_i[1]/1e3:8.1f} t = {ideal["veces_carga"][0]:.0f}-'
          f'{ideal["veces_carga"][1]:.0f} veces la carga útil; turbina CT7: {turbina_i[0]/1e3:.0f}-{turbina_i[1]/1e3:.0f} t')
    print(f'Referencia {P/1e6:g} MW: quema {mdot:.3f} kg/s, eje {P_eje/1e6:.2f} MW, PEM {P_pem/1e3:.0f} kW')
    for c in cajas:
        print(f'  {c[0]:<14}{c[1]:8.0f} -{c[2]:8.0f} kg')
    print(f'  TOTAL        {total[0]:8.0f} -{total[1]:8.0f} kg')
    print(f'  la quema cubre {100*min(cubre):.0f}-{100*max(cubre):.0f} % de la electrólisis ({E_el_t:.0f} kWh/t)')
    for x in plataformas:
        print(f"  {x['corto']:<5} {100*x['frac_carga'][0]:5.0f}-{100*x['frac_carga'][1]:4.0f} % carga; PV "
              + '; '.join(f"{t}: {v['masa_kg']/1e3:.1f} t, {v['dias_por_t']['invierno']:.2f} (invierno)-"
                          f"{v['dias_por_t']['verano']:.2f} (verano) d/t" for t, v in x['pv'].items()))


if __name__ == '__main__':
    main()
