# -*- coding: utf-8 -*-
# [EN] Net mass and energy balance per tonne released, with the compression of the regenerated hydrogen (Section 6).
# =============================================================================
# exp16_balance_neto.py — EXP-16: balance NETO por tonelada liberada, de materia y de energía.
# (Pedido del autor, 2026-09-27: juntar en una sola cuenta lo que EXP-11, EXP-14 y balance_masa
# calculan por separado, y agregar la COMPRESIÓN del H2 regenerado hacia los tanques.)
#
# Arquitectura (b) del autor: electrolizador -> compresor -> tanques (acumulador) -> válvula ->
# celdas de gas. Los tanques salen llenos de tierra; TODO el H2 regenerado se comprime a los tanques
# (aclarado por el autor 2026-09-27: guardado comprimido no cambia la flotabilidad; nada va directo a
# las celdas), y el compresor trabaja al ritmo del electrolizador durante t_regeneracion_h (antes se
# había leído "la mitad a los tanques": corregido). Energía de compresión: trabajo isotérmico de 20 a
# 350 bar / eficiencia típica (DOE Record 9013) y, como extremo pesimista, el promedio medido en
# estaciones (3.1 kWh/kg).
#
# Casos por tonelada liberada:
#   A  cota ideal: 1 t de golpe, condensación instantánea y total (phi = 1);
#   B  1 atm de equilibrio: phi_max a la salida del condensador (40 °C), sin retardo;
#   C  escalones con 10 MW instalados y tau = 5 s (EXP-10), en cada nave grande: los ventiladores
#      andan a plena potencia desde n_tau·tau antes de cada quema hasta su fin (cota del buffer, EXP-11).
# Materia: H2 quemado + aire = agua retenida + vapor venteado + inertes (N2 + Ar + O2 sobrante = 0,
# combustión estequiométrica); regeneración: agua -> H2 + O2; reposición = H2 quemado - regenerado.
# Energía (kWh por tonelada): química quemada -> trabajo (eta) -> ventiladores -> excedente;
# electrólisis y compresión; ELECTRICIDAD EXTERNA NETA = electrólisis + compresión - excedente (motor), o
# electrólisis + compresión + ventiladores (quemador: la batería se recarga desde afuera; hallazgo M7);
# días con paneles (invierno, verano) con esa neta, y energía en tierra del H2 de reposición.
# Masa del compresor: COTA SUPERIOR con un compresor eléctrico de catálogo del tamaño pedido (Hydro-Pac
# C06-05, 260 kg; contraste: Bauer H15.2-20, industrial y más grande).
# =============================================================================

import parametros as p
import rutas
from exp05_condensacion_tiempo_real import phi_max, q_rech, AIRE_ESTEQ
from balance_masa import f_caida_presion

KWH = 3.6e6
P_ESCALON_MW = 10.0
# Revisión 2026-10-01 (R-03): los ventiladores de la unidad de rechazo con la f por caída de presión (en vuelo
# estacionario no hay arrastre de aire de impacto); la f de radiadores de avión (FlyZero) queda como cota alta
F_VENT = f_caida_presion()
F_AVION = p.perdida_rad


def e_comp():
    """Electricidad de compresión [kWh por kg de H2 comprimido]: (típica DOE, medida promedio DOE)."""
    return p.W_iso_20_350 / p.eta_compresion_H2, p.e_comp_medido


def materia(B, R):
    """Balance de materia por tonelada [kg]: B = H2 quemado, R = H2 regenerado (agua retenida / r)."""
    r = 1.0 + p.ratio_O2_H2                      # agua por kg de H2 (8.937)
    aire = B * AIRE_ESTEQ
    return {'H2_quemado': B, 'aire': aire, 'agua_formada': r * B, 'agua_retenida': r * R,
            'vapor_venteado': r * (B - R), 'inertes': aire - p.ratio_O2_H2 * B,
            'H2_regenerado': R, 'O2_liberado': p.ratio_O2_H2 * R, 'H2_reposicion': B - R}


def energia(B, R, E_vent_kWh, caso, extra=None, E_vent_q_kWh=None):
    """Cuenta de energía por tonelada [kWh]: devuelve un dict con cada término y la neta externa, con motor
    (su excedente descuenta) y con quemador (eta = 0: los ventiladores van con la batería, que se recarga desde
    afuera; hallazgo M7, autor 2026-09-29)."""
    E_quim = B * p.Q_comb_H2 / KWH
    W_gen = p.rend_gen * E_quim
    E_el = R * p.E_H2_electrolisis
    ec = e_comp()
    E_comp = [R * x for x in ec]                       # se comprime todo el H2 regenerado
    sobra = [W_gen - v for v in E_vent_kWh]
    neta = [E_el + c - s for c in E_comp for s in sobra]
    E_vent_q = list(E_vent_kWh if E_vent_q_kWh is None else E_vent_q_kWh)
    neta_q = [E_el + c + v for c in E_comp for v in E_vent_q]
    d = {'caso': caso, 'E_quimica': E_quim, 'W_generado': W_gen, 'E_ventiladores': list(E_vent_kWh),
         'excedente': sobra, 'E_electrolisis': E_el, 'E_compresion': E_comp,
         'E_externa_neta': (min(neta), max(neta)),
         'E_ventiladores_quemador': E_vent_q, 'E_externa_neta_quemador': (min(neta_q), max(neta_q)),
         'E_tierra_reposicion': (B - R) * (p.E_H2_electrolisis + ec[0])}
    if extra:
        d.update(extra)
    return d


def main():
    k_rho = p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1)) - 1.0
    r = 1.0 + p.ratio_O2_H2
    casos = []
    # A: cota ideal, phi = 1, sin retardo
    B = 1000.0 / (k_rho + r)
    Q = B * q_rech(p.rend_gen, 1.0)
    Q0 = B * q_rech(0.0, 1.0)                           # con quemador no se exporta trabajo: más calor
    casos.append({'id': 'A', 'nombre': 'cota ideal (1 t de golpe, phi = 1)', 'corto': None, 'B': B, 'R': B,
                  'E_vent': [f * Q / KWH for f in F_VENT], 'E_vent_q': [f * Q0 / KWH for f in F_VENT],
                  'E_vent_av': [f * Q / KWH for f in F_AVION], 'E_vent_q_av': [f * Q0 / KWH for f in F_AVION]})
    # B: 1 atm de equilibrio
    phi = phi_max(p.T_salida_cond)
    B = 1000.0 / (k_rho + r * phi)
    Q = B * q_rech(p.rend_gen, phi)
    Q0 = B * q_rech(0.0, phi)
    casos.append({'id': 'B', 'nombre': '1 atm de equilibrio (phi_max, sin retardo)', 'corto': None, 'B': B,
                  'R': phi * B, 'E_vent': [f * Q / KWH for f in F_VENT],
                  'E_vent_q': [f * Q0 / KWH for f in F_VENT],
                  'E_vent_av': [f * Q / KWH for f in F_AVION], 'E_vent_q_av': [f * Q0 / KWH for f in F_AVION]})
    # C: escalones con 10 MW y tau = 5 s (EXP-10), naves grandes
    for x in rutas.leer_datos('exp10_escalones_condensador')['casos']:
        if x['caso'] != 'ventear_tau5s' or x['P_instalada_MW'] != P_ESCALON_MW or x['corto'] == 'ZNT':
            continue
        n_esc = 1000.0 / x['escalon_kg']
        t_vent = p.n_tau_arranque * p.tau_condensador + x['dur_quema_s']
        E_vent = [f * P_ESCALON_MW * 1e6 * t_vent * n_esc / KWH for f in F_VENT]
        E_vent_av = [f * P_ESCALON_MW * 1e6 * t_vent * n_esc / KWH for f in F_AVION]
        casos.append({'id': 'C', 'nombre': f"escalones de {x['escalon_kg']:.0f} kg con 10 MW, tau 5 s",
                      'corto': x['corto'], 'B': x['H2_por_t_kg'], 'R': min(x['frac_retenida'], 1.0) * x['H2_por_t_kg'],
                      'E_vent': E_vent, 'E_vent_av': E_vent_av, 'n_escalones': n_esc})

    # paneles: energía diaria por nave grande y tecnología (irradiancia sobre el casco, balance_masa)
    ir = rutas.leer_datos('irradiancia_casco')['estaciones']
    A_pv = {x['corto']: x['A_pv_m2'] for x in rutas.leer_datos('balance_masa')['plataformas']}
    filas = []
    for c in casos:
        # escalones (C): los ventiladores andan a la potencia instalada, igual con motor o con quemador
        e = energia(c['B'], c['R'], c['E_vent'], c['id'],
                    {'nombre': c['nombre'], 'corto': c['corto'], 'n_escalones': c.get('n_escalones')},
                    E_vent_q_kWh=c.get('E_vent_q'))
        e['materia'] = materia(c['B'], c['R'])
        e_av = energia(c['B'], c['R'], c['E_vent_av'], c['id'], E_vent_q_kWh=c.get('E_vent_q_av'))
        e['avion'] = {k: e_av[k] for k in ('E_ventiladores', 'excedente', 'E_externa_neta', 'E_ventiladores_quemador',
                                            'E_externa_neta_quemador')}
        dias, dias_q = {}, {}
        for corto in ([c['corto']] if c['corto'] else [d_['corto'] for d_ in p.PLATAFORMAS.values()
                                                       if d_['corto'] != 'ZNT']):
            for tec, (ef, _) in p.PV_TEC.items():
                for est in ('invierno', 'verano'):
                    diario = ir[est]['efectiva'] * ef * A_pv[corto]
                    dias[f'{corto}|{tec}|{est}'] = [n / diario for n in e['E_externa_neta']]
                    dias_q[f'{corto}|{tec}|{est}'] = [n / diario for n in e['E_externa_neta_quemador']]
        e['dias_paneles'] = dias                     # con motor
        e['dias_paneles_quemador'] = dias_q          # con quemador (la batería de los ventiladores también)
        filas.append(e)

    # compresor: flujo y potencia pedidos (regenerar 1 t por día) vs el equipo de referencia
    R_B = next(f for f in filas if f['caso'] == 'B')['materia']['H2_regenerado']
    flujo = R_B / p.t_regeneracion_h                   # al ritmo del electrolizador: no se acumula H2 sin comprimir
    # con ~20 bar de entrada (la presión de partida de la compresión) vale el caudal a 300 psig, el menor listado
    compresor = {'flujo_kgph': flujo, 'P_kW': [flujo * x for x in e_comp()],
                 'ref': p.COMPRESOR_REF, 'cota_valida': flujo <= min(p.COMPRESOR_REF['flujo_kgph']),
                 'masa_cota_sup_kg': p.COMPRESOR_REF['masa_kg']}
    bm = rutas.leer_datos('balance_masa')
    rutas.guardar_datos('exp16_balance_neto', {
        'e_comp_kWhkg': e_comp(), 'filas': filas,
        'compresor': compresor, 'masa_cajas_t': [x / 1e3 for x in bm['total_kg']],
        'masa_con_compresor_t': [x / 1e3 for x in bm['total_kg']]})   # balance_masa ya incluye el compresor
    tabla(filas)
    print(f"compresión: {e_comp()[0]:.2f} (típica DOE) - {e_comp()[1]:.2f} (medida) kWh/kg; "
          f"todo el H2 regenerado va a los tanques")
    for f in filas:
        m = f['materia']
        print(f"[{f['caso']}{' ' + f['corto'] if f['corto'] else ''}] {f['nombre']}")
        print(f"   materia [kg/t]: H2 {m['H2_quemado']:.1f}, aire {m['aire']:.0f}, agua retenida {m['agua_retenida']:.0f}, "
              f"vapor venteado {m['vapor_venteado']:.0f}, regenerado {m['H2_regenerado']:.1f}, reposición {m['H2_reposicion']:.1f}")
        print(f"   energía [kWh/t]: química {f['E_quimica']:.0f}, trabajo {f['W_generado']:.0f}, ventiladores "
              f"{f['E_ventiladores'][1]:.0f}-{f['E_ventiladores'][0]:.0f}, electrólisis {f['E_electrolisis']:.0f}, "
              f"compresión {f['E_compresion'][0]:.0f}-{f['E_compresion'][1]:.0f} -> EXTERNA NETA "
              f"{f['E_externa_neta'][0]:.0f}-{f['E_externa_neta'][1]:.0f}; en tierra (reposición) {f['E_tierra_reposicion']:.0f}")
    print(f"compresor pedido: {compresor['flujo_kgph']:.2f} kg/h, {compresor['P_kW'][0]:.1f}-{compresor['P_kW'][1]:.1f} kW; "
          f"cota superior de masa {p.COMPRESOR_REF['masa_kg'][0]:.0f}-{p.COMPRESOR_REF['masa_kg'][1]:.0f} kg "
          f"(equipo eléctrico de {p.COMPRESOR_REF['flujo_kgph'][0]}-{p.COMPRESOR_REF['flujo_kgph'][1]} kg/h, "
          f"{p.COMPRESOR_REF['motor_kW']} kW)")


def tabla(filas):
    """Tabla en inglés (para integrar al paper cuando el autor decida): cuenta de energía por tonelada."""
    def rg(a, b, d=0):
        a, b = sorted((a, b))
        return f'{a:.{d}f}' if f'{a:.{d}f}' == f'{b:.{d}f}' else f'{a:.{d}f}--{b:.{d}f}'
    L = [r'\begin{tabular}{@{} l r r r r r r r @{}}', r'\toprule',
         r'Case & H\textsubscript{2} burned & Work & Fans & Electrolysis & Compression & Net external & Make-up H\textsubscript{2} \\',
         r' & (kg/t) & (kWh/t) & (kWh/t) & (kWh/t) & (kWh/t) & (kWh/t) & (kg/t) \\', r'\midrule']
    nombres = {'A': 'Ideal bound', 'B': '1 atm, equilibrium'}
    for f in filas:
        n = nombres.get(f['caso'], f"Steps, 10 MW ({f['corto']})")
        L.append(f"{n} & {f['materia']['H2_quemado']:.1f} & {f['W_generado']:.0f} & {rg(*f['E_ventiladores'])} & "
                 f"{f['E_electrolisis']:.0f} & {rg(*f['E_compresion'])} & {rg(*f['E_externa_neta'])} & "
                 f"{f['materia']['H2_reposicion']:.1f}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('art_balance_neto'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/exp16_balance_neto.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
