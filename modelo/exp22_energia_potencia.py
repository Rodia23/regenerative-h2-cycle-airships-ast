# [EN] Energy and power of the release and of the regeneration; delivery law (Sections 5.3 and 6, Figures 6-8).
# =============================================================================
# exp22_energia_potencia.py — EXP-22: energía y potencia de las dos etapas del ciclo, por tonelada liberada (autor,
# 2026-10-12). En correr_todo.py desde 2026-10-13 (después de EXP-21): EXP-21/22 es la fuente de los números de energía
# del paper (sección 6, "Energy and regeneration", y la ley de carga contra tiempo de entrega de la 5.3); sus macros del
# paper van en numeros.tex (exportar_numeros.py) y sus chequeos en verificacion.py. La bitácora
# (bitacora/EXP22_energia_potencia.tex) usa además sus propias macros (\valVeintidos...) y tablas en
# resultados/exploracion/.
#
# Etapa 1, la suelta (combustión, condensación y freno); etapa 2, la regeneración (electrólisis y compresión del H2).
# Dos miradas: la ENERGÍA por tonelada (lo que cuesta; casi no depende del ritmo) y la POTENCIA (lo que pesa el equipo;
# es energía / tiempo, así que depende del tiempo que se elija para cada etapa).
#   A. Ciclo completo por tonelada, los tres arreglos (NHD 14 atm, HD 14 atm, HD 1 atm; balance de EXP-21).
#   B. Suelta con freno (LCA60T, 1 t, 14 atm) contra el tiempo de recuperación: potencia de rechazo, potencia neta del
#      motor y potencia de las hélices al soltar (límite de POTENCIA); energía neta del motor y energía de las hélices,
#      que crece con el tiempo porque sostienen la nave toda la recuperación (límite de ENERGÍA).
#   C. Tamaño de la suelta (0.5-20 t, LCA60T) con 10 MW y con 10 min por tonelada: energía de las hélices por tonelada.
#   D. Regeneración contra su duración: potencia del electrolizador y su masa (PEM, kWkg_PEM), por tonelada.
#   E. La entrega de referencia de la sección 5.3 (EXP-20): el LCA60T reparte su carga de negocio en sueltas de 1 t
#      recuperadas en 10 min, una tras otra, y regenera en la base. Por suelta y por la entrega: rechazo, H2 quemado,
#      agua, hélices, neto a la planta, regeneración y electricidad de tierra; el H2 de los gases en el tiempo (cada
#      suelta deja su agua como lastre: equilibrio con la carga a bordo, minimax.equilibrio). Los tramos entre puntos
#      de entrega y el regreso a la base no se modelan; en la base regenera el PEM del equipo (1 t en
#      t_regeneracion_h), o uno de tierra en el mismo tiempo para toda la entrega.
#   F. Potencia de rechazo por tonelada (P/m = Q/t) contra la suelta, las cuatro naves: sin freno crece como m^0.5;
#      con freno, lo más que cabe en la carga útil, (C - (1 + b) m) / (a m). El cruce es la suelta más grande que se
#      compensa sin freno con un equipo que cabe en la carga útil.
# Fuentes: exp21_balance_energia (balance por tonelada), freno.reparto (hélices: P = F^1.5 / (FM raiz(2 rho A)), con la
# fuerza que el freno toma tras lo que absorbe la quema, y E = P t / 2.5 con la fuerza en rampa a cero),
# exp20_reparto_sueltas (carga de negocio y equipo), fig_presupuesto_masa (a, b, Q), minimax, parametros.
#   Después (2026-10-12/13): la contabilidad con signo de la suelta (hélices y tether) y de la regeneración; H, la
#   carga por viaje contra el tiempo de las sueltas seguidas, N = C' T / (T + a Q); I, el módulo de regeneración
#   independiente de la nave, su nomograma con solo paneles y el costo por tonelada recuperada, k_t m_H2 + a / T.
# Salidas: resultados/figuras/exp22_ciclo_energia.png, exp22_suelta_tiempo.png, exp22_potencia_tiempo.png,
#          exp22_potencia_por_tonelada.png, exp22_entrega_tiempo.png, exp22_suelta_signos_helices.png,
#          exp22_suelta_signos_amarre.png, exp22_regen_signos.png, exp22_carga_tiempo.png, exp22_modulo_regen.png,
#          exp22_regen_nomograma.png; resultados/datos/exp22_*.json;
#          resultados/tablas/exp22_energia_potencia.txt; resultados/exploracion/exp22_numeros.tex y exp22_tabla_*.tex
# =============================================================================
import os
import sys

import numpy as np

import parametros as p
import rutas

KWH = 3.6e6
NAVE = 'Flying Whales LCA60T'
S_REF, T_REF = 1.0, 600.0        # [t], [s] la entrega de referencia de la sección 5.3 (macro RepartoFWDiez)


SUELTAS_HEL = (1.0, 2.0, 10.0, 20.0)   # [t] límites de las hélices por tamaño de suelta, 14 atm (sección 5.3)


def _por(n, x):
    return [n * v for v in x] if isinstance(x, list) else n * x


def _cruce(B, f_diff):
    """El t_s donde f_diff cambia de signo a lo largo de B (interpolación logarítmica); None si no cruza."""
    v = [f_diff(b) for b in B]
    for i in range(len(v) - 1):
        if v[i] * v[i + 1] <= 0:
            a = v[i] / (v[i] - v[i + 1])
            return B[i]['t_s'] * (B[i + 1]['t_s'] / B[i]['t_s']) ** a
    return None


def entrega(d21, neto, Q, cand):
    """E. La entrega de referencia: el LCA60T reparte su carga de negocio en sueltas de S_REF recuperadas en T_REF cada
    una, una tras otra, y regenera en la base. Por suelta y por la entrega completa, y el H2 de los gases tras cada
    suelta (equilibrio con la carga a bordo: la suelta baja S_REF y deja su agua como lastre)."""
    import freno
    import minimax as mm
    r20 = rutas.leer_datos('exp20_reparto_sueltas')
    f1 = next(f for f in r20['naves']['FW']['sueltas'] if abs(f['s_t'] - S_REF) < 1e-9)
    N = f1['N_t'][f'{T_REF:g}']
    n = int(N / S_REF + 1e-9)                                                    # sueltas enteras
    # (D, 2026-10-13) más el tanque de agua, que guarda el agua de toda la carga entregada (w t por tonelada)
    equipo = r20['comp_t'] + S_REF * (r20['b_sin_comp_t_por_t'] + r20['aQ_s'] / T_REF) + r20['w_t_por_t'] * N  # [t]
    P = Q * S_REF / T_REF                                                        # [W] rechazo instalado
    r = freno.reparto(NAVE, S_REF * 1e3, P, cand)
    fl = d21['filas'][1]                                                         # HD, 14 atm
    m, k = fl['masa_kg'], fl['kWh']
    hel = sorted(e / KWH for e in r['E_hel_J'])
    paso = {'t_s': r['t_rec_s'], 'P_rech_W': P, 'H2_kg': S_REF * m['H2_quemado'], 'agua_kg': S_REF * m['agua_retenida'],
            'calor_kWh': sorted(S_REF * x for x in k['calor_rechazado']), 'motor_kWh': [S_REF * x for x in neto],
            'P_hel_W': sorted(r['P_hel_W']), 'hel_kWh': hel,
            'neto_kWh': [S_REF * neto[0] - hel[1], S_REF * neto[1] - hel[0]],          # [pesimista, optimista]
            'electrolisis_kWh': S_REF * k['electrolisis'], 'compresion_kWh': [S_REF * x for x in k['compresion_H2']],
            'H2_regen_kg': S_REF * m['H2_regenerado'], 'H2_repos_kg': S_REF * m['H2_reposicion'],
            'reposicion_kWh': S_REF * k['reposicion_tierra']}
    paso['regen_kWh'] = [paso['electrolisis_kWh'] + c for c in paso['compresion_kWh']]
    paso['externa_kWh'] = [paso['regen_kWh'][0] - paso['neto_kWh'][1], paso['regen_kWh'][1] - paso['neto_kWh'][0]]
    total = {c: _por(n, v) for c, v in paso.items() if c not in ('P_rech_W', 'P_hel_W')}
    P_pem = paso['electrolisis_kWh'] / p.t_regeneracion_h                         # [kW] PEM del equipo: S_REF en un día
    carga0 = n * S_REF + equipo                                                   # [t] a bordo al salir
    cargas = [carga0 - j * (S_REF - paso['agua_kg'] / 1e3) for j in range(n + 1)]
    gas = [mm.equilibrio(NAVE, c * 1e3)['m_gas'] for c in cargas]
    return {'nave': NAVE, 's_t': S_REF, 't_paso_s': T_REF, 'N_t': N, 'n': n, 'equipo_t': equipo, 'carga0_t': carga0,
            'carga_util_t': p.PLATAFORMAS[NAVE]['carga_max'] / 1e3, 'paso': paso, 'total': total, 'P_pem_kW': P_pem,
            't_regen_pem_h': n * S_REF * p.t_regeneracion_h,
            'P_base_MW': [x / p.t_regeneracion_h / 1e3 for x in total['regen_kWh']],   # toda la entrega en un día
            'cargas_t': cargas, 'gas_kg': gas}


def solar(k14, netos):
    """G. El resto de la regeneración de una tonelada (lo que no cubre el neto de la suelta, si se hubiera almacenado)
    en horas de sol pleno (1 kW/m2) y en días de Lima, con los paneles del casco del LCA60T de la sección 6
    (balance_masa: área A_paneles_fraccion del casco; irradiancia efectiva sobre el casco por estación,
    irradiancia_casco). netos: {'hel': [pes, opt], 'amarre': [pes, opt]} [kWh]."""
    bm = rutas.leer_datos('balance_masa')
    A = next(x for x in bm['plataformas'] if x['plataforma'] == NAVE)['A_pv_m2']
    G = {e: x['efectiva'] for e, x in rutas.leer_datos('irradiancia_casco')['estaciones'].items()}   # [kWh/m2/día]
    E_reg = [k14['electrolisis'] + max(k14['compresion_H2']), k14['electrolisis'] + min(k14['compresion_H2'])]
    resto = {c: [E_reg[i] - netos[c][i] for i in (0, 1)] for c in netos}                       # [kWh] pes, opt
    tec = {}
    for t_, (ef, _) in p.PV_TEC.items():
        P = A * ef                                                   # [kW] a 1 kW/m2 (sol pleno)
        tec[t_] = {'ef': ef, 'P_pleno_kW': P,
                   'horas': {c: [r / P for r in resto[c]] for c in resto},
                   'dias': {c: {e: [r / (P * g) for r in resto[c]] for e, g in G.items()} for c in resto}}
    return {'A_m2': A, 'G_kWh_m2_dia': G, 'E_reg_kWh': E_reg, 'resto_kWh': resto, 'tec': tec}


def por_tonelada():
    """F. Potencia de rechazo por tonelada liberada, P/m = Q/t, contra la suelta, las cuatro naves. Sin freno, la del
    tiempo natural de la suelta (q dm*/m, crece como m^0.5); con freno, lo más que cabe en la carga útil,
    (C - (1 + b) m) / (a m). El cruce: la suelta más grande que se compensa sin freno con un equipo que cabe."""
    from scipy.optimize import brentq
    import minimax as mm
    import fig_presupuesto_masa as fpm
    import fig_esferoides as fe
    from exp05_condensacion_tiempo_real import q_rech
    dap = fpm.datos_apilado()
    a, b = dap['t_por_MW'], dap['fijas_t']
    bp, cu = dap['fijas_sin_comp_t'], dap['compresor']             # (2026-10-13) el compresor en unidades de catálogo
    q = q_rech(0.0, p.phi_objetivo)
    naves = {}
    for c in fe.ORDEN:
        nom = next(n for n, v in p.PLATAFORMAS.items() if v['corto'] == c)
        C = p.PLATAFORMAS[nom]['carga_max'] / 1e3

        def sin(m):                                                               # [MW/t]
            return q * mm.resolver(nom, m * 1e3, p.h_max_precision, p.D_lim_oper)['dm_star'] / 1e6 / m

        def cabe(m):                                                              # [MW/t]
            return (C - fpm.compresor_t(m, cu) - (1.0 + bp) * m) / (a * m)

        ms = np.unique(np.append(np.geomspace(0.002, C, 120), min(1.0, C)))
        mx = brentq(lambda m: sin(m) - cabe(m), 1e-4, min(fpm.suelta_unica_max(C, bp, 0.0, cu), 1.0) * 0.999)
        naves[c] = {'nombre': fe.CORTO[c], 'carga_t': C, 'color': p.PLATAFORMAS[nom]['color'], 'm_t': ms.tolist(),
                    'P_sin_MW_t': [sin(x) for x in ms], 'P_1t_MW_t': sin(min(1.0, C)), 'cruce_t': mx,
                    'cruce_MW_t': cabe(mx)}
    return {'a_t_por_MW': a, 'b_t_por_t': b, 'Q_J_por_t': dap['calor_t_J'], 'q_Jkg': q, 'naves': naves}


def datos():
    import exp21_balance_energia as e21
    import freno
    d21 = e21.datos()
    cand = {c['clave']: c for c in rutas.leer_datos('req_masa_rechazo')['candidatos']}
    fs = d21['fs']
    hd14 = d21['filas'][1]
    k = hd14['kWh']
    # neto del motor sin las hélices: trabajo - ventiladores - compresor del escape [kWh por tonelada]
    _n = [k['trabajo_motor'][i] - f * k['calor_rechazado'][i] - k['compresor_escape']
          for i in range(len(d21['etas'])) for f in fs]
    neto = [min(_n), max(_n)]                                       # [pesimista, optimista]
    Q = k['calor_eta0'] * KWH                                     # [J por tonelada]
    # B. barrido del tiempo de recuperación (por la potencia instalada), 1 t
    P_grid = np.geomspace(0.05e6, 400e6, 90)
    B = []
    for P in P_grid:
        r = freno.reparto(NAVE, 1000.0, float(P), cand)
        B.append({'P_W': float(P), 't_s': r['t_rec_s'], 'm_e_kg': r['m_e_kg'], 'P_hel_W': r['P_hel_W'],
                  'E_hel_J': r['E_hel_J']})
    # límites: potencia (hélices al soltar <= neto del motor / t) y energía (hélices <= neto del motor)
    def cruce(f_diff):
        v = [f_diff(b) for b in B]
        for i in range(len(v) - 1):
            if v[i] * v[i + 1] <= 0:
                a = v[i] / (v[i] - v[i + 1])
                return B[i]['t_s'] * (B[i + 1]['t_s'] / B[i]['t_s']) ** a
        return None
    lim = {'potencia_pes': cruce(lambda b: neto[0] * KWH / b['t_s'] - max(b['P_hel_W'])),
           'potencia_opt': cruce(lambda b: neto[1] * KWH / b['t_s'] - min(b['P_hel_W'])),
           'energia_pes': cruce(lambda b: neto[0] * KWH - max(b['E_hel_J'])),
           'energia_opt': cruce(lambda b: neto[1] * KWH - min(b['E_hel_J']))}
    # (autor, 2026-10-13) lo mismo a 1 atm (la figura de potencias del paper): sin compresor del escape; las hélices,
    # las del barrido a 14 atm para el mismo tiempo de recuperación (ver figura_suelta_signos)
    k1 = d21['filas'][2]['kWh']
    _n1 = [k1['trabajo_motor'][i] - f * k1['calor_rechazado'][i] for i in range(len(d21['etas'])) for f in fs]
    neto1 = [min(_n1), max(_n1)]
    lim1 = {'potencia_pes': cruce(lambda b: neto1[0] * KWH / b['t_s'] - max(b['P_hel_W'])),
            'potencia_opt': cruce(lambda b: neto1[1] * KWH / b['t_s'] - min(b['P_hel_W'])),
            'energia_pes': cruce(lambda b: neto1[0] * KWH - max(b['E_hel_J'])),
            'energia_opt': cruce(lambda b: neto1[1] * KWH - min(b['E_hel_J']))}
    # (2026-10-13) el mismo límite de potencia por tamaño de suelta, a 14 atm, para la sección 5.3 (antes salía de EXP-20,
    # con la suelta entera en las hélices): el motor, al ritmo del rechazo que recupera la suelta s en t, entrega
    # s neto / t; las hélices sujetan (s - m_e) g al soltar
    lim_s = {}
    for s_t in SUELTAS_HEL:
        Bs = []
        for t_s in np.geomspace(30.0, 400.0 * 3600.0, 160):
            r = freno.reparto(NAVE, s_t * 1000.0, s_t * Q / t_s, cand)
            Bs.append({'t_s': r['t_rec_s'], 'P_hel_W': r['P_hel_W']})
        lim_s[f'{s_t:g}'] = {
            'potencia_pes': _cruce(Bs, lambda b: s_t * neto[0] * KWH / b['t_s'] - max(b['P_hel_W'])),
            'potencia_opt': _cruce(Bs, lambda b: s_t * neto[1] * KWH / b['t_s'] - min(b['P_hel_W']))}
    # tabla de B en tiempos redondos: el rechazo que recupera la tonelada en t es P = Q / t
    B_tabla = []
    for t_s in (60.0, 600.0, 1800.0, 3600.0, 10800.0, 36000.0):
        r = freno.reparto(NAVE, 1000.0, Q / t_s, cand)
        B_tabla.append({'t_s': r['t_rec_s'], 'P_W': Q / t_s, 'm_e_kg': r['m_e_kg'], 'P_hel_W': sorted(r['P_hel_W']),
                        'E_hel_kWh': sorted(e / KWH for e in r['E_hel_J']), 'P_neto_W': [x * KWH / t_s for x in neto]})
    # C. tamaño de la suelta
    C = []
    for m_t in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
        r10 = freno.reparto(NAVE, m_t * 1000.0, 10e6, cand)
        rT = freno.reparto(NAVE, m_t * 1000.0, m_t * Q / 600.0, cand)       # 10 min por tonelada
        C.append({'m_t': m_t, 't10_s': r10['t_rec_s'], 'E10_kWh_t': sorted(e / KWH / m_t for e in r10['E_hel_J']),
                  'P10_hel_W': r10['P_hel_W'], 'tT_s': rT['t_rec_s'], 'ET_kWh_t': sorted(e / KWH / m_t for e in rT['E_hel_J']),
                  'PT_hel_W': rT['P_hel_W']})
    # D. regeneración: potencia del electrolizador y masa del PEM por tonelada
    E_el = k['electrolisis']                                        # [kWh por tonelada]
    D = [{'t_h': th, 'P_el_kW': E_el / th, 'masa_PEM_t': [E_el / th / max(p.kWkg_PEM) / 1e3,
                                                          E_el / th / min(p.kWkg_PEM) / 1e3]}
         for th in (6.0, 12.0, 24.0, 72.0, 168.0)]
    E_ = entrega(d21, neto, Q, cand)
    G_ = solar(k, {'hel': E_['paso']['neto_kWh'], 'amarre': neto})
    I_ = modulo_regen(d21)
    return {'d21': d21, 'neto_kWh': neto, 'Q_J': Q, 'B': B, 'B_tabla': B_tabla, 'limites_s': lim, 'C': C, 'D': D,
            'neto_1atm_kWh': neto1, 'limites_1atm_s': lim1, 'limites_sueltas_s': lim_s,
            'E_el_kWh': E_el,
            'E_comp_kWh': k['compresion_H2'], 'E': E_, 'F': por_tonelada(), 'G': G_, 'H': carga_tiempo(), 'I': I_,
            'I2': (J_ := presupuesto_regen(I_, G_, d21['filas'][1]['masa_kg']['H2_regenerado'])),
            'I3': capacidad_inversa(J_, I_), 'I2U': (U_ := regen_universal(J_, I_, G_)), 'I4': costo_por_tonelada(J_, U_)}


def chequeos(d):
    import exp21_balance_energia as e21
    out = []

    def ch(n, a, b, tol):
        out.append((n, a, b, abs(a - b) <= tol * max(abs(b), 1e-12)))

    ch('el balance de la suelta cierra (EXP-21)', sum(ok for *_, ok in e21.chequeos(d['d21'])), len(e21.chequeos(d['d21'])),
       0)
    for b in d['B'][::15]:
        ch(f"P x t = calor por tonelada (P = {b['P_W'] / 1e6:.3g} MW)", b['P_W'] * b['t_s'], d['Q_J'], 1e-9)
        ch(f"hélices: E = P t / 2.5 (P = {b['P_W'] / 1e6:.3g} MW)", max(b['E_hel_J']), max(b['P_hel_W']) * b['t_s'] / 2.5,
           1e-12)
    b10 = min(d['B'], key=lambda b: abs(b['P_W'] - 1e7))
    fr = next(x for x in rutas.leer_datos('freno')['filas'] if x['corto'] == 'FW' and x['m_kg'] == 1000.0
              and abs(x['P_W'] - 1e7) < 1.0)
    import freno
    cand = {c['clave']: c for c in rutas.leer_datos('req_masa_rechazo')['candidatos']}
    r = freno.reparto(NAVE, 1000.0, 1e7, cand)
    ch('hélices con 10 MW = la tabla del freno', max(r['E_hel_J']), max(fr['E_hel_J']), 1e-12)
    k14 = d['d21']['filas'][1]['kWh']
    ch('neto sin hélices (pesimista) - hélices con 10 MW = neto de EXP-21 (pesimista)', d['neto_kWh'][0] - k14['helices'][1],
       k14['neto_planta'][0], 1e-9)
    ch('neto sin hélices (optimista) - hélices con 10 MW = neto de EXP-21 (optimista)', d['neto_kWh'][1] - k14['helices'][0],
       k14['neto_planta'][1], 1e-9)
    for x in d['D']:
        ch(f"regeneración en {x['t_h']:g} h: potencia x tiempo = electrólisis por tonelada", x['P_el_kW'] * x['t_h'],
           d['E_el_kWh'], 1e-12)
    ch('regeneración en 24 h: masa del PEM = la caja fija del modelo (clase liviana, φ 0.99)',
       d['D'][2]['masa_PEM_t'][0] * 1e3,
       rutas.leer_datos('req_masa_rechazo')['requisitos']['fijas_kg']['regeneracion'][0] * d['E_el_kWh'] /
       rutas.leer_datos('balance_masa')['E_el_por_t_kWh'], 1e-6)
    # E. la entrega de referencia
    E, ps = d['E'], d['E']['paso']
    ch('entrega: con el rechazo Q s / t la suelta se recupera en t (freno.reparto)', ps['t_s'], T_REF, 1e-9)
    ch('entrega: equipo + carga de negocio = carga útil del LCA60T (ec. del reparto)', E['equipo_t'] + E['N_t'],
       E['carga_util_t'], 1e-12)
    ch('entrega: sueltas enteras, lo que sobra de la carga de negocio es menos de una suelta',
       float(0.0 <= E['N_t'] - E['n'] * S_REF < S_REF), 1.0, 0)
    g = E['gas_kg']
    ch('entrega: el H2 de los gases baja lo quemado en cada suelta (equilibrio, el agua queda como lastre)',
       max(abs(g[j] - g[j + 1] - ps['H2_kg']) for j in range(E['n'])) + ps['H2_kg'], ps['H2_kg'], 1e-9)
    ch('entrega: el PEM del equipo (S_REF en t_regeneracion_h) x su tiempo = la electrólisis de la entrega',
       E['P_pem_kW'] * E['t_regen_pem_h'], E['total']['electrolisis_kWh'], 1e-12)
    ch('entrega: H2 regenerado + reposición = H2 quemado', E['total']['H2_regen_kg'] + E['total']['H2_repos_kg'],
       E['total']['H2_kg'], 1e-9)
    ch('entrega: neto a la planta (pesimista) = k_min Q s de EXP-20 - hélices', ps['neto_kWh'][0],
       rutas.leer_datos('exp20_reparto_sueltas')['k'][0] * d['Q_J'] * S_REF / KWH - ps['hel_kWh'][1], 1e-6)
    # G. el resto de la regeneración con sol: el mismo método que los días por tonelada de la sección 6 (balance_masa)
    Gs = d['G']
    bm = rutas.leer_datos('balance_masa')
    pv_fw = next(x for x in bm['plataformas'] if x['plataforma'] == NAVE)['pv']
    for t_, x in Gs['tec'].items():
        ch(f'solar, {t_}: con la electrólisis a 1 atm, días por tonelada = balance_masa (anual)',
           bm['E_el_por_t_kWh'] / (x['P_pleno_kW'] * Gs['G_kWh_m2_dia']['anual']), pv_fw[t_]['dias_por_t']['anual'], 1e-9)
    # H. carga por viaje contra el tiempo de las sueltas seguidas
    H = d['H']
    for pt in H['puntos']:
        Cp = H['curvas'][pt['s_t']]['C_prima_t'] if pt['s_t'] in H['curvas'] else None
        if Cp is not None:
            ch(f"carga-tiempo, {pt['s_t']:g} t cada {pt['t_s'] / 60:.0f} min: N = C' T / ((1 + w) T + a Q)", pt['N_t'],
               Cp * pt['T_h'] / ((1.0 + H['w']) * pt['T_h'] + H['aQ_h']), 1e-9)
    ch('carga-tiempo: el techo de la suelta única = el de EXP-20 (LCA60T)', H['techo_t'],
       rutas.leer_datos('exp20_reparto_sueltas')['naves']['FW']['unica_max_t'], 1e-12)
    ch('carga-tiempo: el techo de la suelta única = el de la carta de diseño (EXP-18, LCA60T)', H['techo_t'],
       H['techo_paper_t'], 1e-12)
    # (2026-10-13) los límites de las hélices por tamaño de suelta: con 1 t, los del barrido B (otra malla)
    for k_ in ('potencia_pes', 'potencia_opt'):
        ch(f'hélices, 1 t a 14 atm: límite por tamaño de suelta = el del barrido B ({k_})',
           d['limites_sueltas_s']['1'][k_], d['limites_s'][k_], 1e-2)
    # I. el módulo de regeneración: las mismas reglas que las cajas fijas del modelo
    I = d['I']
    bm = rutas.leer_datos('balance_masa')
    ch('módulo: tanque por kg de H2 a 350 bar = el de balance_masa', I['tanques']['350']['kg_por_kg'],
       bm['tanque_por_kg_H2'], 1e-12)
    ch('módulo: masa del PEM por kg/h (clase pesada) = la caja del modelo por kW (balance_masa)',
       I['pem_kg_por_kgph'][1], p.E_H2_electrolisis / min(p.kWkg_PEM), 1e-12)
    ch('módulo: potencia por kg/h x 24 h = la regeneración de una tonelada por kg de H2 (pesimista)',
       I['kW_por_kgph'][1] * d['d21']['filas'][1]['masa_kg']['H2_regenerado'],
       d['E']['paso']['regen_kWh'][1], 1e-9)
    # I2. el módulo alimentado solo con paneles
    J = d['I2']
    m1 = J['tabla']['a-Si']['1|1']
    ch('presupuesto: con 1 kg/h a plena potencia, el H2 de un día = las horas de sol pleno', m1['H2_kg'], J['G_h_dia'],
       1e-12)
    ch('presupuesto: la capacidad del casco = los paneles a sol pleno / energía por kg (a-Si)',
       J['cap_casco_kgph']['a-Si'], d['G']['tec']['a-Si']['P_pleno_kW'] / J['e_kg'], 1e-12)
    m = J['tabla']['CIGS']['12|30']
    ch('presupuesto: total = paneles + PEM + compresor + tanques (CIGS, 12 kg/h, 30 días)', m['total_t'][1],
       m['paneles_t'] + m['pem_t'][1] + m['comp_t'] + m['tanques_t'], 1e-12)
    # I2 universal: con el compresor lleno (c múltiplo exacto de una unidad) el punto real cae sobre 1 + u
    import math
    U = d['I2U']
    c4 = 4 * d['I']['cap_comp_kgph']
    for t_, (ef, kgm2) in p.PV_TEC.items():
        K = U['tec'][t_]['K_kg_por_kgph']
        G, dd = d['I2']['G_h_dia'], 7.0
        M = (c4 * d['I2']['e_kg'] / ef * kgm2 + c4 * U['k_pem'] + math.ceil(c4 / d['I']['cap_comp_kgph'] - 1e-9)
             * d['I']['comp_kg_unidad'] + U['tanque_por_kg'] * c4 * G * dd)
        ch(f'universal, {t_}: con 4 compresores llenos, M / (c K) = 1 + u', M / (c4 * K),
           1.0 + G * dd / U['tec'][t_]['x_star_h'], 1e-9)
    # I4. costo por tonelada recuperada = la masa del módulo (compresor continuo) / las toneladas que recupera
    C4 = d['I4']
    G4, T4 = d['I2']['G_h_dia'], 7.0
    m4 = d['I2']['m_h2_t']
    ef, kgm2 = p.PV_TEC['a-Si']
    c4 = m4 / (G4 * T4)                                                   # [kg/h] para recuperar 1 t en T4 días
    M4 = (c4 * (d['I2']['e_kg'] / ef * kgm2 + p.E_H2_electrolisis / min(p.kWkg_PEM) + d['I2U']['k_comp'])
          + d['I2']['tanque_por_kg'] * m4) / 1e3
    ch('costo por tonelada (a-Si, 1 semana, PEM pesado) = masa del módulo que recupera 1 t en ese tiempo', M4,
       C4['fuentes']['a-Si']['costo']['7'][1], 1e-12)
    ch('costo por tonelada: en T de equilibrio cuesta 1 t por t (red, PEM pesado)',
       C4['piso_t_por_t'] + C4['fuentes']['red']['a'][1] / C4['fuentes']['red']['T_equilibrio_d'][1], 1.0, 1e-12)
    # F. potencia por tonelada
    F = d['F']
    ch('por tonelada: el calor por tonelada de la carta de diseño = el de EXP-21', F['Q_J_por_t'], d['Q_J'], 1e-9)
    ap = {nv['corto']: nv for nv in rutas.leer_datos('fig_presupuesto_masa')['apilado']['naves']}
    for c in ('Hind', 'FW', 'ZR'):
        ch(f"por tonelada, {c}: sin freno, 1 t = el rechazo de la barra NHD del desglose", F['naves'][c]['P_1t_MW_t'],
           ap[c]['sin']['P_MW'], 1e-6)
        ch(f"por tonelada, {c}: la suelta más grande sin freno es menor que 1 t", float(F['naves'][c]['cruce_t'] < 1.0),
           1.0, 0)
    return out


def figuras(d):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401
    from matplotlib.ticker import FixedLocator, NullFormatter
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    # ── A. ciclo completo por tonelada: suelta y regeneración, el caso conservador (eta 0.30, ventiladores 0.20) ──
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    C = {'calor': '#E6BC30', 'vent': '#3CB371', 'comp': '#9b8fd4', 'hel': '#CC79A7', 'neto': '#2b2926',
         'quim': '#4169E1', 'perd': '#e8a09a', 'compH2': '#1EA2E1', 'ext': '#ffffff'}
    eti = {'calor': 'heat rejected', 'vent': 'fans', 'comp': 'exhaust compressor', 'hel': 'propellers',
           'neto': 'net work from the engine', 'quim': 'chemical energy back in the tanks',
           'perd': 'electrolysis losses (heat)', 'compH2': 'hydrogen compression'}
    f_max, i_eta = max(d['d21']['fs']), 0
    x = 0.0
    xt, xl = [], []
    for fl in d['d21']['filas']:
        k = fl['kWh']
        calor = k['calor_rechazado'][i_eta]
        vent = f_max * calor
        trab = k['trabajo_motor'][i_eta]
        hel = k['helices'][1]
        neto = trab - vent - k['compresor_escape'] - hel
        suelta = [('calor', calor), ('vent', vent), ('comp', k['compresor_escape']), ('hel', hel), ('neto', neto)]
        quim = fl['masa_kg']['H2_regenerado'] * p.Q_comb_H2 / KWH
        regen = [('quim', quim), ('perd', k['electrolisis'] - quim), ('compH2', k['compresion_H2'][1])]
        for j, (barra, rot) in enumerate(((suelta, 'release'), (regen, 'regeneration'))):
            xc = x + j * 0.42
            base = 0.0
            for c, v in barra:
                if v <= 0:
                    continue
                ax.bar(xc, v, 0.38, bottom=base, color=C[c], edgecolor='white', lw=0.7)
                if v > 90:
                    ax.text(xc, base + v / 2, f'{v:.0f}', ha='center', va='center', fontsize=7.5,
                            color='white' if c in ('neto', 'quim', 'comp') else '#111111')
                base += v
            ax.text(xc, base + 30, rot, ha='center', fontsize=7.5, color='#444444')
            if j == 1:
                ext = base - neto
                ax.annotate(f'external electricity\n{ext:.0f} kWh', xy=(xc + 0.19, base), xytext=(4, -2),
                            textcoords='offset points', fontsize=7.5, va='top', color='#111111')
        xt.append(x + 0.21)
        xl.append(fl['nombre'].replace(', 10 MW', ''))
        x += 1.35
    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=8.5)
    ax.tick_params(axis='x', length=0)
    ax.grid(axis='x', visible=False)
    ax.set_ylabel('Energy per tonne released (kWh)')
    ax.set_ylim(0, 2900)
    ax.legend(handles=[Patch(color=C[c], label=eti[c]) for c in ('calor', 'vent', 'comp', 'hel', 'neto', 'quim', 'perd',
                                                                   'compH2')],
              loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, fontsize=7.5, frameon=False)
    fig.tight_layout()
    fig.savefig(rutas.figura('exp22_ciclo_energia'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ── B. suelta con freno contra el tiempo de recuperación: potencia (izq.) y energía (der.) ──
    B = d['B']
    t = np.array([b['t_s'] for b in B]) / 60.0
    neto = d['neto_kWh']
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 3.6))
    for x_ in (ax, bx):
        x_.set_xscale('log')
        x_.set_yscale('log')
    ax.plot(t, [b['P_W'] / 1e6 for b in B], color='#E6BC30', lw=1.8, label='heat rejection')
    ax.fill_between(t, neto[0] * KWH / (t * 60) / 1e6, neto[1] * KWH / (t * 60) / 1e6, color='#2b2926', alpha=0.25,
                    lw=0, label='net power of the engine')
    ax.fill_between(t, [min(b['P_hel_W']) / 1e6 for b in B], [max(b['P_hel_W']) / 1e6 for b in B], color='#CC79A7',
                    alpha=0.55, lw=0, label='propellers at the release')
    ax.set_ylabel('Power for 1 t (MW)')
    ax.set_ylim(0.005, 600)
    bx.axhspan(neto[0], neto[1], color='#2b2926', alpha=0.25, lw=0, label='net work of the engine')
    bx.fill_between(t, [min(b['E_hel_J']) / KWH for b in B], [max(b['E_hel_J']) / KWH for b in B], color='#CC79A7',
                    alpha=0.55, lw=0, label='propellers, whole recovery')
    bx.axhline(d['d21']['filas'][1]['kWh']['calor_eta0'], color='#E6BC30', lw=1.8, label='heat rejected')
    bx.set_ylabel('Energy for 1 t (kWh)')
    bx.set_ylim(1, 1e4)
    for x_ in (ax, bx):
        x_.set_xlim(t.min(), t.max())
        x_.xaxis.set_major_locator(FixedLocator([1, 10, 60, 600, 2880]))
        x_.set_xticklabels(['1 min', '10 min', '1 h', '10 h', '48 h'])
        x_.xaxis.set_minor_formatter(NullFormatter())
        x_.set_xlabel('Recovery time of the tonne')
        x_.grid(True, which='major', alpha=0.25)
        x_.legend(loc='lower left' if x_ is ax else 'upper left', fontsize=7, frameon=False)
    for clave, x_, ls in (('potencia_pes', ax, ':'), ('potencia_opt', ax, '--'), ('energia_pes', bx, ':'),
                          ('energia_opt', bx, '--')):
        v = d['limites_s'][clave]
        if v:
            x_.axvline(v / 60, color='#555555', lw=0.9, ls=ls)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(rutas.figura('exp22_suelta_tiempo'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    figura_potencia_tiempo(d)
    figura_por_tonelada(d)
    figura_entrega(d)
    for dispositivo in ('helices', 'amarre'):   # (autor, 2026-10-12) con hélices y con tether; 14 y 1 atm
        for atm in (14, 1):
            figura_suelta_signos(d, dispositivo, atm)
    figura_potencia_paper(d, 1)                   # (autor, 2026-10-13) para el paper: solo potencias, 1 atm
    figura_regen_signos(d)
    figura_carga_tiempo(d)
    figura_carga_naves(d)                          # (autor, 2026-10-13) para el paper, 4x1 en toneladas
    figura_modulo_regen(d)
    figura_regen_nomograma(d)


# H. Carga por viaje contra el tiempo de las sueltas seguidas (autor, 2026-10-12). De la ecuación del reparto (EXP-20),
# N = C - m_c - s (b' + a Q / t), y del tiempo de la entrega, T = n t con n = N / s sueltas (continuo), sale
#     N = C' T / (T + a Q),   C' = C - m_c - b' s:
# en a Q horas de sueltas seguidas se entrega la mitad de C'. Cada tamaño de suelta empieza en la suelta única (N = s,
# con su recuperación mínima). Las configuraciones que el autor marcó (s, t) van como puntos.
T_MAX = 168.0                    # [h] el eje de la figura H llega a una semana
PUNTOS_H = ((0.5, 300.0), (1.0, 300.0), (1.0, 600.0), (2.0, 600.0), (5.0, 3600.0), (10.0, 3600.0), (10.0, 21600.0),
            (20.0, 21600.0))
SUELTAS_H = (0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 24.0, 30.0)   # (autor, 2026-10-12) con las del máximo (24 y 30 t)
                                                                 # y una de 100 kg, casi sobre la envolvente s -> 0


# I. El módulo de regeneración, independiente de la nave por ahora (autor, 2026-10-12): PEM + compresor de H2 de alta
# presión + tanques, dimensionados a tope de capacidad con energía y agua ideales. Por kg de H2 por hora: potencia
# (e_el + e_comp) y masa (PEM con kWkg_PEM; compresor en unidades de catálogo, COMPRESOR_REF). Los tanques: masa y
# volumen por kg de H2 guardado (tanque_a_presion), a 350 y 700 bar.
C_EJE = np.linspace(0.05, 10.0, 400)            # [kg/h] capacidad a tope


def modulo_regen(d21):
    """I. El módulo de regeneración, solo, independiente de la nave (autor, 2026-10-12): a tope de capacidad, con energía
    y agua ideales. Por kg de H2 por hora: la potencia del PEM (e_el) y del compresor de H2 (e_comp), la masa del PEM
    (kWkg_PEM) y la del compresor (unidades de catálogo, COMPRESOR_REF). Los tanques, aparte: masa y volumen por kg de
    H2 guardado (tanque_a_presion: el sistema Tipo IV del DOE, modular, así que escala lineal con lo guardado), a 350 y
    a 700 bar."""
    import math
    e_el = p.E_H2_electrolisis
    e_comp = sorted(d21['e_comp_kWhkg'])                                 # [típica DOE, medida] [kWh/kg]
    cap = min(p.COMPRESOR_REF['flujo_kgph'])                              # [kg/h] por unidad de catálogo
    m_u = p.COMPRESOR_REF['masa_kg'][0]                                   # [kg] por unidad

    def modulo(c):                                                        # c [kg/h] a tope de capacidad
        P_pem = c * e_el                                                  # [kW]
        n = math.ceil(c / cap - 1e-9)
        return {'c_kgph': c, 'P_pem_kW': P_pem, 'P_comp_kW': [c * x for x in e_comp],
                'pem_kg': [P_pem / max(p.kWkg_PEM), P_pem / min(p.kWkg_PEM)], 'n_comp': n, 'comp_kg': n * m_u}
    eje = [modulo(float(c)) for c in C_EJE]
    tanques = {}
    for P in (350.0, 700.0):
        t = p.tanque_a_presion(P)
        tanques[f'{P:g}'] = {'w': t['w'], 'kg_por_kg': (1.0 - t['w']) / t['w'],
                             'L_por_kg': t['V_L'] / p.TANQUE_DOE24['util_kg'], 'util_kg_modulo': p.TANQUE_DOE24['util_kg']}
    return {'e_el_kWh_kg': e_el, 'e_comp_kWh_kg': e_comp, 'kW_por_kgph': [e_el + x for x in e_comp],
            'pem_kg_por_kgph': [e_el / max(p.kWkg_PEM), e_el / min(p.kWkg_PEM)], 'cap_comp_kgph': cap,
            'comp_kg_unidad': m_u, 'eje': eje, 'tanques': tanques}


def figura_modulo_regen(d):
    """El módulo de regeneración a tope de capacidad (parte I), independiente de la nave: izquierda, la potencia que
    pide (PEM y compresor de H2) contra su capacidad en kg de H2 por hora; centro, su masa (compresor de catálogo en
    escalones y PEM, entre la clase pesada y la liviana); derecha, la masa de los tanques contra el H2 que guardan, a 350
    y a 700 bar."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    I = d['I']
    c = np.array([x['c_kgph'] for x in I['eje']])
    pem = np.array([x['P_pem_kW'] for x in I['eje']])
    cmp_ = np.array([x['P_comp_kW'][1] for x in I['eje']])
    m_comp = np.array([x['comp_kg'] for x in I['eje']]) / 1e3
    m_pem = [np.array([x['pem_kg'][i] for x in I['eje']]) / 1e3 for i in (0, 1)]
    C = {'pem': '#E6BC30', 'comp': '#a8841c', 'tanq': '#d9c27a', 'neto': '#111111'}
    fig, (ap, am, at) = plt.subplots(1, 3, figsize=(7.4, 3.3))
    # potencia
    ap.fill_between(c, 0, pem, color=C['pem'], lw=0)
    ap.fill_between(c, pem, pem + cmp_, color=C['comp'], lw=0)
    ap.plot(c, pem + cmp_, color=C['neto'], lw=1.5)
    ap.set_xlabel('Capacity (kg H$_2$/h)')
    ap.set_ylabel('Power at full capacity (kW)')
    k = I['kW_por_kgph']
    ap.annotate(f"{k[0]:.1f}–{k[1]:.1f} kW\nper kg/h", xy=(c[-1] * 0.05, (pem + cmp_)[-1] * 0.8), fontsize=7.5)
    # masa del módulo
    am.fill_between(c, 0, m_comp, color=C['comp'], lw=0)
    am.fill_between(c, m_comp, m_comp + m_pem[1], color=C['pem'], lw=0)
    am.plot(c, m_comp + m_pem[1], color=C['neto'], lw=1.5)
    am.plot(c, m_comp + m_pem[0], color=C['neto'], lw=1.1, ls=(0, (4, 2)))
    am.set_xlabel('Capacity (kg H$_2$/h)')
    am.set_ylabel('Mass of the module (t)')
    am.annotate(f"compressor: {I['comp_kg_unidad']:.0f} kg\nper {I['cap_comp_kgph']:.2f} kg/h", xy=(c[-1] * 0.05,
                (m_comp + m_pem[1])[-1] * 0.8), fontsize=7.5)
    # tanques
    H2 = np.linspace(0, 100, 50)
    for P, ls in (('350', '-'), ('700', (0, (4, 2)))):
        x = I['tanques'][P]
        at.plot(H2, H2 * x['kg_por_kg'] / 1e3, color=C['neto'], lw=1.5, ls=ls,
                label=f"{P} bar: {x['kg_por_kg']:.1f} kg/kg, {x['L_por_kg']:.0f} L/kg")
    at.set_xlabel('Hydrogen stored (kg)')
    at.set_ylabel('Mass of the tanks (t)')
    at.legend(loc='lower right', frameon=False)
    for ax in (ap, am, at):
        ax.grid(True, alpha=0.25)
        ax.set_xlim(0, None)
        ax.set_ylim(0, None)
    fig.legend(handles=[Patch(color=C['pem'], label='PEM stack'), Patch(color=C['comp'], label='hydrogen compressor'),
                        Line2D([], [], color=C['neto'], lw=1.5, label='in all (dashed: lighter PEM)')],
               loc='lower center', ncol=3, frameon=False)
    plot_style.normalizar_paper(fig)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(rutas.figura('exp22_modulo_regen'), dpi=300, bbox_inches='tight')
    plt.close(fig)


# I2. El módulo alimentado solo con paneles: lo que cuesta en carga útil (autor, 2026-10-12). Para una capacidad c
# [kg de H2 por hora a plena potencia], los paneles se dimensionan para esa potencia pico (sin batería), así que el
# módulo rinde por día como si trabajara G horas a tope (G = horas de sol pleno por día, Lima, promedio anual):
#   masa = c (paneles + PEM + compresor) + tanques vacíos que esperan todo lo regenerado en el intervalo,
#   tanques = tanque_por_kg x c G días;
# en % de la carga útil del LCA60T (solo como referencia de escala). Paneles con la energía pesimista por kg.
CAP_I2 = (1.0, 3.0, 6.0, 12.0)                  # [kg/h]
DIAS_I2 = np.geomspace(0.25, 30.0, 200)         # [días] intervalo de regeneración


def presupuesto_regen(I, Gs, m_h2_t):
    import math
    C = p.PLATAFORMAS[NAVE]['carga_max'] / 1e3
    e_kg = max(I['kW_por_kgph'])                                          # [kWh/kg] pesimista
    G = Gs['G_kWh_m2_dia']['anual']                                       # [h de sol pleno por día]
    tq = I['tanques']['350']['kg_por_kg']

    def masas(c, dias, tec):
        ef, kgm2 = p.PV_TEC[tec]
        pv = c * e_kg / ef * kgm2 / 1e3                                   # [t]
        pem = [c * p.E_H2_electrolisis / k / 1e3 for k in sorted(p.kWkg_PEM, reverse=True)]   # [liviano, pesado]
        comp = math.ceil(c / I['cap_comp_kgph'] - 1e-9) * I['comp_kg_unidad'] / 1e3
        h2 = c * G * dias                                                 # [kg] regenerado en el intervalo
        tanq = tq * h2 / 1e3
        return {'paneles_t': pv, 'pem_t': pem, 'comp_t': comp, 'tanques_t': tanq, 'H2_kg': h2,
                'total_t': [pv + x + comp + tanq for x in pem], 'area_m2': c * e_kg / ef}
    curvas = {tec: {f'{c:g}': [masas(c, float(x), tec) for x in DIAS_I2] for c in CAP_I2} for tec in p.PV_TEC}
    tabla = {tec: {f'{c:g}|{dd:g}': masas(c, dd, tec) for c in (1.0, 5.0, 12.0) for dd in (1.0, 7.0, 30.0)}
             for tec in p.PV_TEC}
    return {'C_t': C, 'e_kg': e_kg, 'G_h_dia': G, 'tanque_por_kg': tq, 'dias': DIAS_I2.tolist(), 'curvas': curvas,
            'm_h2_t': m_h2_t,
            'tabla': tabla, 'cap_casco_kgph': {tec: Gs['A_m2'] * p.PV_TEC[tec][0] / e_kg for tec in p.PV_TEC},
            'tanque_t_por_t': tq * m_h2_t / 1e3}                    # tanque por tonelada liberada regenerada [t/t]


def capacidad_inversa(J, I):
    """I2, lectura inversa (autor, 2026-10-12): con un presupuesto de masa y unos días de viaje, la mayor capacidad
    (kg/h) cuyo módulo cabe, y lo que regenera. Los tanques van vacíos: el H2 sale del agua que ya va a bordo como
    lastre (por kg de H2 se electrolizan 1 + ratio_O2_H2 kg de agua y se ventea el oxígeno), así que no se suma."""
    import math
    G, e_kg, tq = J['G_h_dia'], J['e_kg'], J['tanque_por_kg']

    def masa(c, dias, tec):
        ef, kgm2 = p.PV_TEC[tec]
        pv = c * e_kg / ef * kgm2
        pem = c * p.E_H2_electrolisis / min(p.kWkg_PEM)                     # PEM de la clase pesada
        comp = math.ceil(c / I['cap_comp_kgph'] - 1e-9) * I['comp_kg_unidad']
        h2 = c * G * dias
        return {'total_t': (pv + pem + comp + tq * h2) / 1e3, 'paneles_t': pv / 1e3, 'pem_t': pem / 1e3,
                'comp_t': comp / 1e3, 'tanques_t': tq * h2 / 1e3, 'H2_kg': h2}
    filas = []
    for tec in p.PV_TEC:
        for presup in (2.0, 6.0, 10.0, 20.0):
            for dias in (1.0, 2.0, 7.0):
                cs = np.arange(0.01, 80.0, 0.01)
                ok = [c for c in cs if masa(c, dias, tec)['total_t'] <= presup]
                c = max(ok) if ok else 0.0
                filas.append({'tec': tec, 'presupuesto_t': presup, 'dias': dias, 'c_kgph': c, **masa(c, dias, tec)})
    return {'agua_por_kg': 1.0 + p.ratio_O2_H2, 'O2_por_kg': p.ratio_O2_H2, 'filas': filas}


def costo_por_tonelada(J, U):
    """Cierre de la regeneración (autor, 2026-10-13): la carga útil que cuesta recuperar a bordo una tonelada de
    sustentación (el H2 de una tonelada liberada) en T días:
        costo = k_t m_H2 + a / T,   a = K m_H2 / G   [t por t; t·día por t],
    con K la parte de potencia por kg/h (paneles + PEM + compresor continuo) y G las horas de trabajo a tope por día
    (sol pleno con paneles; 24 h con la red de tierra, sin paneles). El primer término, los tanques, es el piso: no
    depende del tiempo. Cuesta menos que lo que recupera (1 t por t) desde T = a / (1 - k_t m_H2)."""
    m, tq, e_kg = J['m_h2_t'], J['tanque_por_kg'], J['e_kg']
    piso = tq * m / 1e3
    fuentes = {'a-Si': (e_kg / p.PV_TEC['a-Si'][0] * p.PV_TEC['a-Si'][1], J['G_h_dia']),
               'CIGS': (e_kg / p.PV_TEC['CIGS'][0] * p.PV_TEC['CIGS'][1], J['G_h_dia']),
               'red': (0.0, 24.0)}
    out = {}
    for f, (k_pv, G) in fuentes.items():
        a = [(k_pv + p.E_H2_electrolisis / k + U['k_comp']) * m / G / 1e3 for k in sorted(p.kWkg_PEM, reverse=True)]
        out[f] = {'G_h_dia': G, 'a': a, 'T_equilibrio_d': [x / (1.0 - piso) for x in a],
                  'costo': {f'{T:g}': [piso + x / T for x in a] for T in (1.0, 2.0, 7.0, 30.0)}}
    return {'piso_t_por_t': piso, 'fuentes': out}


def regen_universal(J, I, Gs):
    """I2 universal (autor, 2026-10-13): la masa del módulo es proporcional a la capacidad c, así que
        M = c (K + k_t G días),   K = paneles + PEM + compresor por kg/h (compresor continuo: masa / caudal de la unidad),
    y con u = (horas de sol pleno acumuladas, G días) / x*, x* = K / k_t, queda una sola curva para cualquier capacidad,
    panel, estación o lugar:  M / (c K) = 1 + u.  x* es donde los tanques pesan lo mismo que la parte de potencia.
    Los puntos reales llevan el compresor en unidades enteras: se separan de la curva con capacidades chicas."""
    import math
    e_kg, tq = J['e_kg'], J['tanque_por_kg']
    k_comp = I['comp_kg_unidad'] / I['cap_comp_kgph']                       # [kg por kg/h] continuo
    k_pem = p.E_H2_electrolisis / min(p.kWkg_PEM)                            # [kg por kg/h] PEM pesado
    tec = {}
    for t_, (ef, kgm2) in p.PV_TEC.items():
        K = e_kg / ef * kgm2 + k_pem + k_comp                                # [kg por kg/h]
        tec[t_] = {'K_kg_por_kgph': K, 'x_star_h': K / tq, 'x_star_dias_anual': K / tq / J['G_h_dia']}
    puntos = []
    for t_, (ef, kgm2) in p.PV_TEC.items():
        K = tec[t_]['K_kg_por_kgph']
        for est, G in Gs['G_kWh_m2_dia'].items():
            for c in CAP_I2:
                for dias in (0.25, 1.0, 3.0, 7.0, 14.0, 30.0):
                    M = (c * e_kg / ef * kgm2 + c * k_pem + math.ceil(c / I['cap_comp_kgph'] - 1e-9)
                         * I['comp_kg_unidad'] + tq * c * G * dias)
                    puntos.append({'tec': t_, 'estacion': est, 'c_kgph': c, 'dias': dias,
                                   'u': G * dias / tec[t_]['x_star_h'], 'y': M / (c * K)})
    return {'k_comp': k_comp, 'k_pem': k_pem, 'tanque_por_kg': tq, 'tec': tec, 'puntos': puntos}


CAP_N = (1.0, 2.0, 5.0, 10.0, 20.0)             # [kg de H2 por hora de sol pleno] capacidades del nomograma
HSOL_N = (10.0, 30.0, 100.0, 300.0)            # [h de sol pleno] líneas de igual tiempo


def figura_regen_nomograma(d):
    """Nomograma del módulo de regeneración (autor, 2026-10-13), las tres variables a la vista: el H2 recuperado
    (abajo kg; arriba toneladas de suelta), la masa total del módulo y la capacidad (kg de H2 por hora de sol pleno).
        M = c (paneles + PEM) + compresor + k_t H,   H = c h   (h: horas de sol pleno acumuladas).
    Líneas llenas: capacidad fija; discontinuas: horas de sol fijas (con su equivalente en días de Lima). Compresor
    continuo por encima de una unidad (masa / caudal de la unidad), una unidad como mínimo. PEM de la clase pesada."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.ticker import FixedLocator, NullFormatter
    from matplotlib.lines import Line2D
    import fig_esferoides as fe
    J, I, U = d['I2'], d['I'], d['I2U']
    tq, e_kg, G, m_h2_t = J['tanque_por_kg'], J['e_kg'], J['G_h_dia'], J['m_h2_t']
    k_pem = p.E_H2_electrolisis / min(p.kWkg_PEM)

    def masa(c, H, k_pv):                                                    # [t]
        comp = max(I['comp_kg_unidad'], U['k_comp'] * c)
        return (c * (k_pv + k_pem) + comp + tq * H) / 1e3

    cmap = plt.get_cmap('viridis')
    cols = [cmap(x) for x in np.linspace(0.85, 0.1, len(CAP_N))]
    Hx = np.geomspace(3.0, 2500.0, 300)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 4.6), sharey=True)
    cargas = sorted(((v['carga_max'] / 1e3, fe.CORTO[v['corto']]) for v in p.PLATAFORMAS.values()))
    for ax, tec in zip(axs, ('a-Si', 'CIGS')):
        ef, kgm2 = p.PV_TEC[tec]
        k_pv = e_kg / ef * kgm2
        for c, col in zip(CAP_N, cols):
            ax.plot(Hx, [masa(c, H, k_pv) for H in Hx], color=col, lw=1.8, zorder=3)
        for h in HSOL_N:
            # solo las capacidades del nomograma; el rótulo, en el extremo bajo, donde las líneas están separadas
            cs = np.geomspace(min(CAP_N), max(CAP_N), 200)
            ys = np.array([masa(c, c * h, k_pv) for c in cs])
            ax.plot(cs * h, ys, color='#6b675f', lw=0.9, ls=(0, (4, 2)), zorder=2)
            ax.plot(cs[0] * h, ys[0], 'o', ms=2.5, color='#6b675f', zorder=2)
            # (autor, 2026-10-13) los días en las tres estaciones: invierno / año / verano en Lima
            Ge = d['G']['G_kWh_m2_dia']
            # en el panel a-Si las curvas se juntan abajo: el rótulo un poco a la derecha del punto (autor, 2026-10-13)
            dx, ha = (5, 'left') if tec == 'a-Si' else (0, 'center')
            ax.annotate(f"{h:g} h\n{h / Ge['invierno']:.0f}/{h / Ge['anual']:.0f}/{h / Ge['verano']:.0f} d",
                        xy=(cs[0] * h, ys[0]), xytext=(dx, -4), textcoords='offset points', fontsize=7, ha=ha,
                        va='top', color='#55534e')
        for C_, nom in cargas:
            ax.axhline(C_, color='#b8b6b0', lw=0.8, ls=(0, (1, 1.5)), zorder=1)
        if ax is axs[1]:                                                      # las naves, como marcas a la derecha
            der = ax.secondary_yaxis('right')
            der.set_yticks([C_ for C_, _ in cargas])
            der.set_yticklabels([f'{nom} ({C_:.0f} t)' for C_, nom in cargas], fontsize=plot_style.T_ROTULO,
                                color='#8a8882')
            der.tick_params(axis='y', length=0)
        for x_ in (ax,):
            x_.set_xscale('log')
            x_.set_yscale('log')
        ax.set_xlim(Hx[0], Hx[-1])
        ax.xaxis.set_major_locator(FixedLocator([10, 100, 1000]))
        ax.set_xticklabels(['10', '100', '1000'])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel('Hydrogen recovered (kg)')
        ax.set_title(f"{tec} panels, {100 * ef:.1f} % efficient", fontsize=8, loc='left', pad=22)
        ax.grid(True, which='major', alpha=0.25)
        top = ax.secondary_xaxis('top', functions=(lambda x: x / m_h2_t, lambda t: t * m_h2_t))
        top.set_xticks([0.1, 0.3, 1, 3, 10, 30])
        top.set_xticklabels(['0.1', '0.3', '1', '3', '10', '30'])
        top.set_xlabel('Tonnes of release recovered', fontsize=8)
    axs[0].set_ylabel('Mass of the regeneration module (t)')
    axs[0].set_ylim(0.3, 70)
    axs[0].yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20, 50]))
    axs[0].set_yticklabels(['0.5', '1', '2', '5', '10', '20', '50'])
    axs[0].yaxis.set_minor_formatter(NullFormatter())
    fig.legend(handles=[Line2D([], [], color=col, lw=1.8, label=f'{c:g} kg H$_2$ per hour of sun')
                        for c, col in zip(CAP_N, cols)] +
                       [Line2D([], [], color='#6b675f', lw=0.9, ls=(0, (4, 2)), label='equal hours of full sun')],
               loc='lower center', bbox_to_anchor=(0.5, 0.05), ncol=3, frameon=False)
    fig.text(0.5, 0.005, f"Days in Lima, winter / year / summer: {d['G']['G_kWh_m2_dia']['invierno']:.1f} / {G:.1f} / "
                         f"{d['G']['G_kWh_m2_dia']['verano']:.1f} h of full sun per day; the curves do not depend on "
                         "the season.\nPEM of the heavy class; compressor at least one catalog unit.", ha='center',
             va='bottom', fontsize=7.5, color='#333333')
    plot_style.normalizar_paper(fig)                 # todos los rótulos >= T_ROTULO (Elsevier pide >= 7 pt)
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    plot_style.guardar_paper(fig, rutas.figura('exp22_regen_nomograma'))   # figura del paper: PNG y, para el envío, PDF
    plt.close(fig)


def carga_tiempo():
    # (2026-10-13) el compresor de H2 en unidades de catálogo (fig_presupuesto_masa.compresor_t): m_c(s) es una unidad
    # mientras alcance; la envolvente de arriba (s -> 0) lleva una
    import fig_presupuesto_masa as fpm
    r20 = rutas.leer_datos('exp20_reparto_sueltas')
    C = p.PLATAFORMAS[NAVE]['carga_max'] / 1e3
    mc, bp, aQ, cu = r20['comp_t'], r20['b_sin_comp_t_por_t'], r20['aQ_s'], r20['compresor']
    w = r20['w_t_por_t']                                               # (D) tanque de agua, por tonelada entregada
    curvas, no_caben = {}, []
    for s in SUELTAS_H:
        cs = fpm.compresor_t(s, cu)
        if C - cs - s * (1.0 + bp + w) <= 0.0:                          # ni con recuperación infinita (30 t)
            no_caben.append(s)
            continue
        t_min = s * aQ / (C - cs - s * (1.0 + bp + w))                  # [s] la suelta única (N = s)
        t = np.geomspace(t_min * 1.0001, 3.0e6, 400)
        N = (C - cs - s * (bp + aQ / t)) / (1.0 + w)
        T = N * t / s / 3600.0                                          # [h]
        curvas[s] = {'t_min_s': t_min, 'comp_t': cs, 'C_prima_t': C - cs - bp * s, 'T_h': T.tolist(),
                     'N_t': N.tolist()}
    puntos = []
    for s, t in PUNTOS_H:
        N = (C - fpm.compresor_t(s, cu) - s * (bp + aQ / t)) / (1.0 + w)
        puntos.append({'s_t': s, 't_s': t, 'N_t': N, 'T_h': N * t / s / 3600.0, 'ritmo_t_h': s / t * 3600.0})
    # la suelta única (N = s) de cada tamaño, y su techo con recuperación infinita (el mismo de la carta de diseño)
    techo = fpm.suelta_unica_max(C, bp + w, 0.0, cu)
    su = np.linspace(0.05, techo * 0.9995, 400)
    unica = {'s_t': su.tolist(), 'T_h': [x * aQ / (C - fpm.compresor_t(x, cu) - x * (1.0 + bp + w)) / 3600.0 for x in su]}
    # la envolvente superior: sueltas cada vez más chicas, s -> 0 (C' -> C - m_c); la inferior es la suelta única
    Te = np.geomspace(0.05, T_MAX, 400)
    envolvente = {'T_h': Te.tolist(), 'N_t': ((C - mc) * Te / ((1.0 + w) * Te + aQ / 3600.0)).tolist()}
    return {'C_t': C, 'aQ_h': aQ / 3600.0, 'w': w, 'curvas': curvas, 'no_caben_t': no_caben, 'puntos': puntos,
            'unica': unica,
            'techo_t': techo,
            'techo_paper_t': rutas.leer_datos('exp18_por_nave')['diseno']['naves']['FW']['limite_t'],
            'envolvente': envolvente, 'C_menos_mc_t': C - mc}


REDONDOS = (0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75)   # [t] sueltas "bonitas" por debajo de 1 t
# las sueltas de la figura del paper, elegidas por el autor (2026-10-13) entre las que da sueltas_parejas; fijas para que
# no cambien con el modelo (los chequeos avisan si alguna deja de caber)
SUELTAS_PAPER = {'Hind': (2.0, 5.0, 7.0), 'FW': (8.0, 15.0, 23.0), 'ZR': (0.5, 1.0, 2.0), 'ZNT': (0.05, 0.1, 0.2)}


def _redondo(s):
    """La suelta redonda más cercana: toneladas enteras desde 1 t; por debajo, de REDONDOS."""
    return float(round(s)) if s >= 1.0 else min(REDONDOS, key=lambda r: abs(r - s))


def sueltas_parejas(C, bp, aQ_h, T_fin, cu, w=0.0):
    """Tres sueltas que reparten parejo, al final del eje (T_fin), el espacio entre las envolventes: en T_fin la curva
    de la suelta s vale (C - m_c(s) - b' s) f, f = T/(T + aQ); arriba (C - m_c) f, con una unidad, abajo la suelta única
    que se recupera en T_fin. Se toman 1/4, 2/4 y 3/4 de ese tramo y se redondean a una suelta bonita (autor,
    2026-10-13; desde ese día con el compresor en unidades de catálogo)."""
    import fig_presupuesto_masa as fpm
    f = T_fin / ((1.0 + w) * T_fin + aQ_h)
    s_u = fpm.suelta_unica_max(C, bp + w, aQ_h / T_fin, cu)                # la suelta única que se recupera en T_fin
    arriba, abajo = (C - fpm.compresor_t(0.0, cu)) * f, s_u
    malla = np.linspace(1e-4, s_u, 4000)
    Nf = np.array([(C - fpm.compresor_t(x, cu) - bp * x) * f for x in malla])
    return [_redondo(float(malla[np.argmin(abs(Nf - (arriba - k * (arriba - abajo) / 4.0)))])) for k in (1, 2, 3)]


def figura_carga_naves(d):
    """Figura para el paper (autor, 2026-10-13): la carga por viaje contra el tiempo de las sueltas seguidas, una nave
    por fila (4x1), en toneladas, con el compresor de catálogo fijo (como la Fig. 6). Envolventes: arriba, sueltas cada
    vez más chicas, N = (C - m_c) T / (T + aQ), del color de la nave y discontinua; abajo, la suelta única (N = s), en
    negro. Entre ellas, tres sueltas repartidas parejo (sueltas_parejas), rotuladas junto al punto donde nacen.
    (2026-10-13) el compresor de H2 en unidades de catálogo (fig_presupuesto_masa.compresor_t); las sueltas, las que
    eligió el autor (SUELTAS_PAPER)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    from matplotlib.lines import Line2D
    import fig_esferoides as fe
    import fig_presupuesto_masa as fpm
    r20 = rutas.leer_datos('exp20_reparto_sueltas')
    aQ, bp, mc, cu = r20['aQ_s'] / 3600.0, r20['b_sin_comp_t_por_t'], r20['comp_t'], r20['compresor']
    w = r20['w_t_por_t']                                               # (D) tanque de agua
    T = np.geomspace(0.05, T_MAX, 500)
    fig, axs = plt.subplots(4, 1, figsize=(7.4, 9.6), sharex=True)
    for ax, c in zip(axs, ('Hind', 'FW', 'ZR', 'ZNT')):
        nom = next(n for n, v in p.PLATAFORMAS.items() if v['corto'] == c)
        C = p.PLATAFORMAS[nom]['carga_max'] / 1e3
        col = p.PLATAFORMAS[nom]['color']
        ax.plot(T, (C - mc) * T / ((1.0 + w) * T + aQ), color=col, lw=1.2, ls=(0, (4, 2)))
        su = np.linspace(1e-4, fpm.suelta_unica_max(C, bp + w, 0.0, cu) * 0.9995, 1200)
        ax.plot([x * aQ / (C - fpm.compresor_t(x, cu) - x * (1.0 + bp + w)) for x in su], su, color='#111111', lw=1.3)
        for s in SUELTAS_PAPER[c]:
            cs = fpm.compresor_t(s, cu)
            T0 = s * aQ / (C - cs - s * (1.0 + bp + w))
            Tf = T[T >= T0]
            ax.plot(Tf, (C - cs - bp * s) * Tf / ((1.0 + w) * Tf + aQ), color=col, lw=1.6)
            ax.plot(T0, s, 'o', ms=6.5, color=col, mec='white', mew=0.8, zorder=5)
            # abajo a la derecha del punto (bajo la suelta única no pasa ninguna línea); pegado al cero, arriba
            arriba = s < 0.15 * C
            ax.annotate(f'{s:g} t', xy=(T0, s), xytext=(-5, 5) if arriba else (5, -6), textcoords='offset points',
                        fontsize=7.5, ha='right' if arriba else 'left', va='bottom' if arriba else 'top',
                        color='#222222', zorder=6)
        ax.axhline(C, color='#8a8882', lw=0.8, ls=(0, (1, 1.5)))
        ax.set_xscale('log')
        ax.set_xlim(T[0], T[-1])
        ax.set_ylim(0, 1.1 * C)
        ax.set_title(f'{fe.CORTO[c]}, payload {C:.0f} t' if C >= 1.5 else f'{fe.CORTO[c]}, payload {C:.1f} t',
                     fontsize=8, loc='left')
        ax.set_ylabel('Cargo per trip (t)')
        ax.grid(True, alpha=0.25)
    axs[-1].xaxis.set_major_locator(FixedLocator([0.1, 0.25, 0.5, 1, 2, 4, 8, 24, 72, 168]))
    axs[-1].set_xticklabels(['6 min', '15 min', '30 min', '1 h', '2 h', '4 h', '8 h', '1 day', '3 days', '1 week'])
    axs[-1].xaxis.set_minor_formatter(NullFormatter())
    for ax in axs:
        ax.xaxis.set_minor_locator(NullLocator())
        ax.tick_params(axis='x', labelbottom=True)
    axs[-1].set_xlabel('Time of the consecutive releases')
    fig.legend(handles=[Line2D([], [], color='#666666', lw=1.6, marker='o', ms=6.5, mec='white',
                               label='releases of the size marked at the dot (the single release of that size)'),
                        Line2D([], [], color='#666666', lw=1.2, ls=(0, (4, 2)), label='releases → 0'),
                        Line2D([], [], color='#111111', lw=1.3, label='single release'),
                        Line2D([], [], color='#8a8882', lw=0.8, ls=(0, (1, 1.5)), label='payload')],
               loc='lower center', ncol=2, frameon=False)
    plot_style.normalizar_paper(fig)
    fig.tight_layout(rect=(0, 0.05, 1, 1), h_pad=0.8)
    plot_style.guardar_paper(fig, rutas.figura('exp22_carga_naves'))
    plt.close(fig)


def figura_carga_tiempo(d):
    """La carga por viaje contra el tiempo de las sueltas seguidas, LCA60T, una curva por tamaño de suelta (parte H)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.ticker import FixedLocator, NullFormatter
    from exp20_reparto_sueltas import _tono
    H = d['H']
    col = p.PLATAFORMAS[NAVE]['color']
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    tonos = np.linspace(0.65, -0.45, len(SUELTAS_H))                 # sueltas chicas claras, grandes oscuras
    fin = []                                                         # (N al final, rótulo): se juntan los que se tocan
    for (s, cv), w in zip(H['curvas'].items(), tonos):
        T, N = np.array(cv['T_h']), np.array(cv['N_t'])
        ok = T <= T_MAX
        ax.plot(T[ok], N[ok], color=_tono(col, w), lw=1.7)
        ax.plot(T[0], N[0], 'o', ms=4, color=_tono(col, w), mec='white', mew=0.6)        # la suelta única
        if fin and abs(fin[-1][0] - N[ok][-1]) < 1.5:
            fin[-1] = (fin[-1][0], fin[-1][1] + f', {s:g}')
        else:
            fin.append((N[ok][-1], f'{s:g}'))
    for y, lab in fin:
        ax.annotate(f'{lab} t', xy=(T_MAX, y), xytext=(3, 0), textcoords='offset points', fontsize=7.5, va='center',
                    color='#333333')
    # la suelta única de cada tamaño (donde nace cada curva) y su techo teórico (recuperación infinita)
    un = H['unica']
    ax.plot(un['T_h'], un['s_t'], color='#c0392b', lw=1.4, zorder=4)
    # la envolvente superior (s -> 0) y la región de operación entre las dos envolventes
    ev = H['envolvente']
    Te, Ne = np.array(ev['T_h']), np.array(ev['N_t'])
    ax.plot(Te, Ne, color='#111111', lw=1.2, ls=(0, (4, 2)), zorder=4)
    abajo = np.interp(Te, np.array(un['T_h']), np.array(un['s_t']), left=0.0, right=H['techo_t'])
    ax.fill_between(Te, abajo, Ne, color='#d9d6cf', alpha=0.35, lw=0, zorder=0)
    i = int(np.argmin(abs(Te - 2.5)))
    ax.annotate(f"upper envelope, releases s → 0:\nN = {H['C_menos_mc_t']:.1f} t · T / ((1 + w) T + a Q)\n"
                f"lower envelope: the single release (red)", xy=(Te[i], Ne[i]), xytext=(0.12, 45),
                textcoords='data', fontsize=7.5, ha='left', va='bottom', color='#111111',
                arrowprops=dict(arrowstyle='-', color='#888888', lw=0.7, shrinkA=0, shrinkB=2))
    ax.axhline(H['techo_t'], color='#c0392b', lw=0.9, ls=(0, (1, 1.5)))
    ax.annotate(f"largest single release: {H['techo_t']:.1f} t (recovery → ∞)", xy=(0.12, H['techo_t']), xytext=(0, 3),
                textcoords='offset points', fontsize=7.5, color='#c0392b', va='bottom')
    ax.axvline(H['aQ_h'] / (1.0 + H['w']), color='#8a8882', lw=0.9, ls=(0, (4, 3)))
    ax.annotate(f"a Q / (1 + w) = {H['aQ_h'] / (1.0 + H['w']):.1f} h:\nhalf of C′ / (1 + w)",
                xy=(H['aQ_h'] / (1.0 + H['w']), 2), xytext=(4, 0), textcoords='offset points',
                fontsize=7.5, color='#6b675f', va='bottom')
    ax.axhline(H['C_t'], color='#52514e', lw=0.9, ls=(0, (4, 3)))
    ax.annotate(f"payload of the LCA60T ({H['C_t']:.0f} t)", xy=(0.12, H['C_t']), xytext=(0, 3),
                textcoords='offset points', fontsize=7.5, color='#52514e', va='bottom')
    # los puntos marcados, en grupos de carga parecida (un solo rótulo por grupo, arriba a la izquierda)
    grupos = []
    for pt in sorted(H['puntos'], key=lambda x: x['N_t']):
        if grupos and abs(grupos[-1][-1]['N_t'] - pt['N_t']) < 1.5:
            grupos[-1].append(pt)
        else:
            grupos.append([pt])
    for g in grupos:
        lineas = []
        for pt in sorted(g, key=lambda x: x['T_h']):
            t = pt['t_s']
            lab = f"{pt['s_t']:g} t every {t / 60:.0f} min" if t < 3600 else f"{pt['s_t']:g} t every {t / 3600:.0f} h"
            lineas.append(f"{lab}: {pt['N_t']:.1f} t in {pt['T_h']:.1f} h")
            ax.plot(pt['T_h'], pt['N_t'], 'D', ms=5.5, color='#c0392b', mec='white', mew=0.7, zorder=5)
        x0 = min(pt['T_h'] for pt in g)
        y0 = max(pt['N_t'] for pt in g)
        ax.annotate('\n'.join(lineas), xy=(x0, y0), xytext=(-8, 14), textcoords='offset points', fontsize=7.5,
                    color='#222222', ha='right', va='bottom',
                    arrowprops=dict(arrowstyle='-', color='#888888', lw=0.7, shrinkA=0, shrinkB=3))
    ax.set_xscale('log')
    ax.set_xlim(0.1, T_MAX)
    ax.set_ylim(0, 66)
    ax.xaxis.set_major_locator(FixedLocator([0.1, 0.25, 0.5, 1, 2, 4, 8, 24, 72, 168]))
    ax.set_xticklabels(['6 min', '15 min', '30 min', '1 h', '2 h', '4 h', '8 h', '1 day', '3 days', '1 week'])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('Time of the consecutive releases')
    ax.set_ylabel('Delivered cargo per trip (t)')
    ax.grid(True, alpha=0.25)
    plot_style.normalizar_paper(fig)
    fig.tight_layout()
    fig.savefig(rutas.figura('exp22_carga_tiempo'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def figura_regen_signos(d):
    """La regeneración del H2 presurizado de una tonelada (una ventana: lo que cabe en los tanques), sola y solo en
    energía (autor, 2026-10-12): la electricidad acumulada dentro de la ventana al ritmo del modelo (t_regeneracion_h por
    tonelada), con la misma convención que la figura de la suelta (los gastos de cada equipo abajo de cero: la
    electrólisis del PEM y la compresión del H2 a los tanques; la energía química del H2 no entra como término, igual que
    en la suelta). Líneas verticales: dónde se habría acabado el neto de la suelta si se hubiera almacenado (con hélices a
    10 min y con tether; pesimista y optimista). Eje derecho (pedido del autor): la masa de H2 regenerado y guardado en
    los tanques. Franjas con la compresión pesimista. Es la misma con hélices o con tether."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.ticker import FixedLocator
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import fig_presupuesto_masa as fpm
    ep = {n.split(' ')[0]: v[0] for n, v in fpm.ESTILO_PAPER.items()}
    # el PEM y el compresor de H2 son cajas fijas (mostaza en la Fig. 5); el H2 de los tanques, el azul del H2 (Fig. 1)
    C = {'el': ep['Fixed'], 'comp': '#a8841c', 'neto': '#111111', 'h2': '#3a63a8', 'hel': '#CC79A7', 'amarre': '#555555'}
    k = d['d21']['filas'][1]['kWh']
    m_h2 = d['d21']['filas'][1]['masa_kg']['H2_regenerado']
    E_el = k['electrolisis']
    comp_pes, comp_opt = max(k['compresion_H2']), min(k['compresion_H2'])
    tr = p.t_regeneracion_h
    tau = np.linspace(0.0, tr, 200)                                 # [h] tiempo dentro de la ventana
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    n1 = -E_el * tau / tr
    n2 = n1 - comp_pes * tau / tr
    ax.fill_between(tau, 0, n1, color=C['el'], alpha=0.85, lw=0)
    ax.fill_between(tau, n1, n2, color=C['comp'], alpha=0.85, lw=0)
    ax.plot(tau, n1, color='white', lw=0.7)
    ax.plot(tau, n2, color=C['neto'], lw=1.8)
    ax.axhline(0, color='#444444', lw=0.8)
    # el neto de la suelta, si se hubiera almacenado, se acaba a las E_neto / E_regen x t_r
    netos = {'hel': d['E']['paso']['neto_kWh'], 'amarre': d['neto_kWh']}            # [pesimista, optimista] [kWh]
    E_reg = [E_el + comp_pes, E_el + comp_opt]
    acaba = {c: [netos[c][i] / E_reg[i] * tr for i in (0, 1)] for c in netos}       # [h]
    LS = {'hel': (0, (1, 1.2)), 'amarre': (0, (5, 1.5, 1, 1.5))}
    for c in netos:
        for x in acaba[c]:
            ax.axvline(x, color=C[c], lw=1.6, ls=LS[c], zorder=4)
    y_rot = -1350
    for i, caso in enumerate(('pessimistic', 'optimistic')):
        x = max(acaba['hel'][i], acaba['amarre'][i])
        ax.annotate(f"{caso}:\n{acaba['hel'][i]:.2f} h propellers\n{acaba['amarre'][i]:.2f} h tether", xy=(x, y_rot),
                    xytext=(5, 0), textcoords='offset points', fontsize=7.5, va='center',
                    bbox=dict(fc='white', ec='none', alpha=0.9, pad=1.5))
    resto = d['G']['resto_kWh']
    ax.annotate(f"−{E_reg[1]:.0f} to −{E_reg[0]:.0f} kWh in all;\nnot covered by the release: "
                f"{min(min(v) for v in resto.values()):.0f}–{max(max(v) for v in resto.values()):.0f} kWh",
                xy=(tr, -E_reg[0]), xytext=(-4, -6), textcoords='offset points', fontsize=7.5, ha='right', va='top')
    ax.set_xlim(0, tr)
    ax.set_ylim(-3000, 250)
    ax.xaxis.set_major_locator(FixedLocator([0, 3, 6, 9, 12, 15, 18, 21, 24]))
    ax.set_xlabel(f'Time into the regeneration, {tr:.0f} h per tonne (h)')
    ax.set_ylabel('Electricity used, cumulative (kWh)')
    ax.grid(True, alpha=0.25)
    # eje derecho: el H2 regenerado y guardado en los tanques
    bx = ax.twinx()
    bx.plot(tau, m_h2 * tau / tr, color=C['h2'], lw=2.0)
    bx.set_ylim(0, 1.15 * m_h2)
    bx.set_ylabel('Hydrogen stored in the tanks (kg)', color=C['h2'])
    bx.tick_params(axis='y', colors=C['h2'])
    bx.annotate(f'{m_h2:.1f} kg', xy=(tr, m_h2), xytext=(-4, 4), textcoords='offset points', fontsize=7.5,
                color=C['h2'], ha='right', va='bottom')
    fig.legend(handles=[Patch(color=C['el'], label='− electrolysis (regeneration unit, PEM)'),
                        Patch(color=C['comp'], label='− hydrogen compression to the tanks (350 bar)'),
                        Line2D([], [], color=C['neto'], lw=1.8, label='electricity in all'),
                        Line2D([], [], color=C['h2'], lw=2.0, label='hydrogen stored in the tanks (right axis)'),
                        Line2D([], [], color=C['hel'], lw=1.6, ls=LS['hel'],
                               label='net of the release with propellers (10 min), if stored: used up'),
                        Line2D([], [], color=C['amarre'], lw=1.6, ls=LS['amarre'],
                               label='net of the release with a tether, if stored: used up')],
               loc='lower center', ncol=2, frameon=False)
    plot_style.normalizar_paper(fig)
    for tx in (bx.yaxis.label,):
        tx.set_color(C['h2'])
    fig.tight_layout(rect=(0, 0.2, 1, 1))
    fig.savefig(rutas.figura('exp22_regen_signos'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def figura_potencia_paper(d, atm=1):
    """Figura para el paper (autor, 2026-10-13): las potencias de la suelta de 1 t contra el tiempo de recuperación, a
    1 atm, pesimista y optimista; las Figs. de la suelta con hélices y con tether juntas, solo en potencia. Arriba de
    cero el trabajo del motor; abajo los ventiladores y las hélices al soltar. Línea llena: el neto con hélices (punto:
    donde cruza cero); discontinua: el neto con tether, que no gasta energía. A 1 atm no hay compresor del escape; las
    hélices son las del barrido a 14 atm para el mismo tiempo (ver figura_suelta_signos)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.ticker import FixedLocator, NullFormatter
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import fig_presupuesto_masa as fpm
    i_arr = 1 if atm == 14 else 2
    k = d['d21']['filas'][i_arr]['kWh']
    etas, fs, A = d['d21']['etas'], d['d21']['fs'], p.freno_A_discos
    B = d['B']
    t = np.array([b['t_s'] for b in B])
    ep = {n.split(' ')[0]: v[0] for n, v in fpm.ESTILO_PAPER.items()}
    C = {'motor': ep['Conversion'], 'vent': ep['Heat'], 'comp': ep['Exhaust'], 'hel': '#CC79A7', 'neto': '#111111'}
    casos = ((f"Pessimistic\nengine efficiency {min(etas):.2f}, fans {100 * max(fs):.0f} % of the heat, "
              f"{min(A):.0f} m² disks", int(np.argmin(etas)), max(fs), np.array([max(b['P_hel_W']) for b in B]) / 1e6),
             (f"Optimistic\nengine efficiency {max(etas):.2f}, fans {100 * min(fs):.1f} % of the heat, "
              f"{max(A):.0f} m² disks", int(np.argmax(etas)), min(fs), np.array([min(b['P_hel_W']) for b in B]) / 1e6))
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.5), sharey=True)
    for ax, (tit, i, f, hel) in zip(axs, casos):
        W = k['trabajo_motor'][i] * KWH / t / 1e6                            # [MW] media en la recuperación
        V = f * k['calor_rechazado'][i] * KWH / t / 1e6
        Wc = k['compresor_escape'] * KWH / t / 1e6
        n1, n2 = -V, -V - Wc
        n3 = n2 - hel
        ax.fill_between(t, 0, W, color=C['motor'], alpha=0.85, lw=0)
        ax.fill_between(t, 0, n1, color=C['vent'], alpha=0.85, lw=0)
        if atm == 14:
            ax.fill_between(t, n1, n2, color=C['comp'], alpha=0.85, lw=0)
        ax.fill_between(t, n2, n3, color=C['hel'], alpha=0.85, lw=0)
        ax.plot(t, n2, color='white', lw=0.7)
        ax.plot(t, W + n3, color=C['neto'], lw=1.8, zorder=4)                # con hélices
        ax.plot(t, W + n2, color=C['neto'], lw=1.3, ls=(0, (4, 2)), zorder=4)  # con tether
        ax.axhline(0, color='#444444', lw=0.8)
        neto = W + n3
        s = np.where(np.diff(np.sign(neto)) != 0)[0]
        if len(s):
            i0 = s[0]
            tc = t[i0] * (t[i0 + 1] / t[i0]) ** (neto[i0] / (neto[i0] - neto[i0 + 1]))
            ax.plot(tc, 0, 'o', ms=4, color=C['neto'], zorder=5)
            lab = f'{tc / 60:.0f} min' if tc < 3600 else f'{tc / 3600:.1f} h'
            ax.annotate(f'net = 0 at {lab}', xy=(tc, 0), xytext=(4, 6), textcoords='offset points', fontsize=7.5)
        ax.axvline(S_REF * T_REF, color='#c0392b', lw=0.9, alpha=0.6)
        ax.set_xscale('log')
        ax.set_xlim(30, 2.0e5)
        ax.set_ylim(-1.5, 1.5)
        ax.xaxis.set_major_locator(FixedLocator([60, 600, 3600, 36000, 172800]))
        ax.set_xticklabels(['1 min', '10 min', '1 h', '10 h', '48 h'])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel('Recovery time of the tonne')
        ax.set_title(tit, fontsize=8, loc='left')
        ax.grid(True, alpha=0.25)
    axs[0].set_ylabel(f'Power for 1 t at {atm} atm (MW)')
    axs[0].annotate(f'{T_REF / 60:.0f} min', xy=(T_REF, 1.35), xytext=(4, 0), textcoords='offset points', fontsize=7.5,
                    color='#c0392b', va='center')
    manos = [Patch(color=C['motor'], label='+ engine work'), Patch(color=C['vent'], label='− fans')]
    if atm == 14:
        manos.append(Patch(color=C['comp'], label='− exhaust compressor'))
    manos += [Patch(color=C['hel'], label='− propellers at the release'),
              Line2D([], [], color=C['neto'], lw=1.8, label='net, propellers'),
              Line2D([], [], color=C['neto'], lw=1.3, ls=(0, (4, 2)), label='net, tether')]
    fig.legend(handles=manos, loc='lower center', ncol=3, frameon=False)
    plot_style.normalizar_paper(fig)
    for ax in axs:
        ax.title.set_fontsize(8)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    plot_style.guardar_paper(fig, rutas.figura(f'exp22_potencia_suelta_{atm}atm'))
    plt.close(fig)


def figura_suelta_signos(d, dispositivo='helices', atm=14):
    """La suelta de 1 t contra el tiempo de recuperación como contabilidad con signo (autor, 2026-10-12): arriba de cero
    lo que el motor da (su trabajo), abajo lo que se gasta (ventiladores, compresor del escape y, con hélices, las
    hélices); la curva negra es el neto a la planta y el punto, donde cruza cero (sin línea vertical). Filas: energía por
    tonelada y potencia (la media del motor en la recuperación contra el pico de las hélices al soltar: el criterio del
    límite de potencia), con la potencia en eje lineal recortado a +-1.5 MW (la symlog era casi lo mismo). Columnas: caso
    pesimista y optimista. dispositivo: 'helices' o 'amarre' (el tether sujeta sin gastar energía: sin esa franja).
    atm: 14 (caso base, HD 14 atm de EXP-21) o 1 (HD 1 atm: sin compresor del escape, phi_max; autor, 2026-10-13). A
    1 atm las hélices son las del barrido de 14 atm para el mismo tiempo de recuperación: su fuerza depende de lo que
    la quema absorbe en tiempo real, que cambia ~4 % entre las dos presiones (EXP-21 toma la misma aproximación)."""
    con_hel = dispositivo == 'helices'
    i_arr = 1 if atm == 14 else 2
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.ticker import FixedLocator, NullFormatter
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    k = d['d21']['filas'][i_arr]['kWh']
    etas, fs, A = d['d21']['etas'], d['d21']['fs'], p.freno_A_discos
    B = d['B']
    t = np.array([b['t_s'] for b in B])
    cero = np.zeros_like(t)
    discos = (f', {min(A):.0f} m² disks', f', {max(A):.0f} m² disks') if con_hel else ('', '')
    casos = (
        (f"Pessimistic, {atm} atm\nengine efficiency {min(etas):.2f}, fans {100 * max(fs):.0f} % of the heat"
         f"{discos[0]}",
         int(np.argmin(etas)), max(fs),
         np.array([max(b['E_hel_J']) for b in B]) / KWH if con_hel else cero,
         np.array([max(b['P_hel_W']) for b in B]) / 1e6 if con_hel else cero),
        (f"Optimistic, {atm} atm\nengine efficiency {max(etas):.2f}, fans {100 * min(fs):.1f} % of the heat"
         f"{discos[1]}",
         int(np.argmax(etas)), min(fs),
         np.array([min(b['E_hel_J']) for b in B]) / KWH if con_hel else cero,
         np.array([min(b['P_hel_W']) for b in B]) / 1e6 if con_hel else cero))
    # los colores de los equipos en la Fig. 5 del manuscrito (fig_presupuesto_masa.ESTILO_PAPER, variante A): motor
    # (conversión) verde, ventiladores (unidad de rechazo) azul, compresor del escape celeste; el freno no está en esa
    # figura: rosado
    import fig_presupuesto_masa as fpm
    ep = {n.split(' ')[0]: v[0] for n, v in fpm.ESTILO_PAPER.items()}
    C = {'motor': ep['Conversion'], 'vent': ep['Heat'], 'comp': ep['Exhaust'], 'hel': '#CC79A7', 'neto': '#111111'}
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.2), sharex=True)
    for j, (tit, i, f, hel_E, hel_P) in enumerate(casos):
        W, V, Wc = k['trabajo_motor'][i], f * k['calor_rechazado'][i], k['compresor_escape']
        for fila, ax in enumerate(axs[:, j]):
            esc = (lambda x: x * np.ones_like(t)) if fila == 0 else (lambda x: x * KWH / t / 1e6)   # kWh o MW
            hel = hel_E if fila == 0 else hel_P
            pos, n1 = esc(W), -esc(V)
            n2 = n1 - esc(Wc)
            n3 = n2 - hel
            ax.fill_between(t, 0, pos, color=C['motor'], alpha=0.85, lw=0)
            ax.fill_between(t, 0, n1, color=C['vent'], alpha=0.85, lw=0)
            ax.fill_between(t, n1, n2, color=C['comp'], alpha=0.85, lw=0)
            ax.fill_between(t, n2, n3, color=C['hel'], alpha=0.85, lw=0)
            for borde in (n1, n2):                                       # separadores blancos, como en la Fig. 5
                ax.plot(t, borde, color='white', lw=0.7)
            neto = pos + n3
            ax.plot(t, neto, color=C['neto'], lw=1.8)
            ax.axhline(0, color='#444444', lw=0.8)
            s = np.where(np.diff(np.sign(neto)) != 0)[0]
            if len(s):
                i0 = s[0]
                tc = t[i0] * (t[i0 + 1] / t[i0]) ** (neto[i0] / (neto[i0] - neto[i0 + 1]))
                ax.plot(tc, 0, 'o', ms=3.5, color=C['neto'], zorder=5)
                lab = f'{tc / 60:.0f} min' if tc < 3600 else f'{tc / 3600:.1f} h'
                ax.annotate(f'net = 0 at {lab}', xy=(tc, 0), xytext=(4, 6), textcoords='offset points', fontsize=7.5)
            ax.set_xscale('log')
            if fila == 1:
                ax.set_ylim(-1.5, 1.5)
            else:
                ax.set_ylim(-1300 if con_hel else -500, 800)
            ax.set_xlim(30, 2.0e5)
            ax.xaxis.set_major_locator(FixedLocator([60, 600, 3600, 36000, 172800]))
            ax.set_xticklabels(['1 min', '10 min', '1 h', '10 h', '48 h'])
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.axvline(S_REF * T_REF, color='#c0392b', lw=0.9, alpha=0.6)
            ax.grid(True, alpha=0.25)
        axs[0, j].set_title(tit, fontsize=8, loc='left')
    axs[0, 0].set_ylabel('Energy for 1 t (kWh)')
    axs[1, 0].set_ylabel('Power for 1 t (MW)')
    for ax in axs[1]:
        ax.set_xlabel('Recovery time of the tonne')
    axs[0, 0].annotate(f'reference delivery: {T_REF / 60:.0f} min', xy=(T_REF, 650), xytext=(4, 0),
                       textcoords='offset points', fontsize=7.5, color='#c0392b')
    manos = [Patch(color=C['motor'], label='+ engine work (conversion unit)'),
             Patch(color=C['vent'], label='− fans (heat-rejection unit)'),
             Patch(color=C['comp'], label='− exhaust compressor (14 atm)') if atm == 14 else
             Line2D([], [], ls='none', label='(no exhaust compressor at 1 atm)')]
    manos += ([Patch(color=C['hel'], label='− propellers (hold-down device)'),
               Line2D([], [], color=C['neto'], lw=1.8, label='net to the plant; dot: net = 0')] if con_hel else
              [Line2D([], [], color=C['neto'], lw=1.8, label='net to the plant (tether: holds without energy)')])
    fig.legend(handles=manos, loc='lower center', ncol=3 if con_hel else 2, frameon=False)
    plot_style.normalizar_paper(fig)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(rutas.figura(f'exp22_suelta_signos_{dispositivo}' + ('' if atm == 14 else '_1atm')), dpi=300,
                bbox_inches='tight')
    plt.close(fig)


def figura_potencia_tiempo(d):
    """Potencia de rechazo contra la duración de la quema, 1 t: toda forma de soltar la tonelada cae sobre P t = Q (la
    misma energía a cualquier ritmo); el freno la reparte en el tiempo y baja la potencia, y con ella la masa."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.ticker import FixedLocator, NullFormatter
    Q = d['Q_J']
    ap = {nv['corto']: nv for nv in rutas.leer_datos('fig_presupuesto_masa')['apilado']['naves']}
    col = {v['corto']: v['color'] for v in p.PLATAFORMAS.values()}
    b10 = rutas.leer_datos('exp21_balance_energia')['filas'][1]
    puntos = [('no hold-down device (ideal bound), LCA60T', ap['FW']['sin']['t_nat_s'], ap['FW']['sin']['P_MW'], '#c0392b'),
              ('hold-down, largest that fits (LCA60T)', ap['FW']['max']['t_rec_min'] * 60, ap['FW']['max']['P_MW'], col['FW']),
              ('hold-down, 10 MW', b10['t_s'], 10.0, col['FW']),
              ('hold-down, largest that fits (LZ 129)', ap['Hind']['max']['t_rec_min'] * 60, ap['Hind']['max']['P_MW'],
               col['Hind']),
              ('hold-down, largest that fits (ZR-3)', ap['ZR']['max']['t_rec_min'] * 60, ap['ZR']['max']['P_MW'], col['ZR'])]
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    ax.set_xscale('log')
    ax.set_yscale('log')
    t = np.geomspace(3, 2e4, 200)
    ax.plot(t, Q / t / 1e6, color='#555555', lw=1.4, ls=(0, (5, 2)), zorder=2)
    fig.canvas.draw()
    pp = ax.transData.transform([(5, Q / 5 / 1e6), (50, Q / 50 / 1e6)])
    ang = np.degrees(np.arctan2(pp[1, 1] - pp[0, 1], pp[1, 0] - pp[0, 0]))
    ax.annotate(f'P × t = {Q / 1e9:.2f} GJ per tonne: the same energy at any pace', xy=(6, Q / 6 / 1e6), xytext=(4, -8),
                textcoords='offset points', fontsize=7.5, color='#444444', rotation=ang, rotation_mode='anchor', va='top')
    OFF = [((10, 2), 'left'), ((10, 4), 'left'), ((-10, -2), 'right'), ((10, 6), 'left'), ((-10, -4), 'right')]
    for (nom, ts, P, c), ((dx, dy), ha) in zip(puntos, OFF):
        ax.plot(ts, P, 'o', ms=7, color=c, mec='white', mew=1.2, zorder=4)
        ax.annotate(f"{nom}\n{P:.3g} MW, {ts:.0f} s" if ts < 120 else f"{nom}\n{P:.3g} MW, {ts / 60:.1f} min",
                    xy=(ts, P), xytext=(dx, dy), textcoords='offset points', fontsize=7.5,
                    ha=ha, va='bottom' if dy > 0 else 'top', color='#222222')
    ax.set_xlim(3, 2e4)
    ax.set_ylim(0.3, 2000)
    ax.xaxis.set_major_locator(FixedLocator([10, 60, 600, 3600, 14400]))
    ax.set_xticklabels(['10 s', '1 min', '10 min', '1 h', '4 h'])
    ax.yaxis.set_major_locator(FixedLocator([1, 10, 100, 1000]))
    ax.set_yticklabels(['1', '10', '100', '1000'])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('Duration of the burn (recovery of the tonne)')
    ax.set_ylabel('Heat-rejection power (MW)')
    ax.grid(True, which='major', alpha=0.25)
    plot_style.normalizar_paper(fig)
    fig.tight_layout()
    fig.savefig(rutas.figura('exp22_potencia_tiempo'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def figura_por_tonelada(d):
    """Potencia de rechazo por tonelada liberada contra la suelta (log-log), las cuatro naves (parte F)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.ticker import FixedLocator, NullFormatter
    from matplotlib.lines import Line2D
    F = d['F']
    a, b, Q = F['a_t_por_MW'], F['b_t_por_t'], F['Q_J_por_t']
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.set_yscale('log')
    ax.set_xscale('log')
    for t_s, lab in ((60, '1 min'), (600, '10 min'), (3600, '1 h'), (36000, '10 h')):
        y = Q / t_s / 1e6
        ax.axhline(y, color='#b8b6b0', lw=0.9, ls=(0, (4, 3)), zorder=1)
        izq = t_s >= 3600                                   # los de abajo a la izquierda, donde no pasan curvas
        ax.annotate(f'hold-down device, {lab} per tonne: {y:.3g} MW/t', xy=(0.0022 if izq else 55, y), xytext=(0, 3),
                    textcoords='offset points', fontsize=7.5, color='#8a8882', ha='left' if izq else 'right',
                    va='bottom')
    orden = ('ZNT', 'ZR', 'FW', 'Hind')
    for c in orden:
        nv = F['naves'][c]
        C, ms = nv['carga_t'], np.array(nv['m_t'])
        ax.plot(ms, nv['P_sin_MW_t'], color=nv['color'], lw=1.8, zorder=3)
        ax.plot(min(1.0, C), nv['P_1t_MW_t'], 'o', ms=5, color=nv['color'], mec='white', mew=1.0, zorder=5)
        mf = np.geomspace(0.002, C / (1 + b) * 0.99999, 700)
        Pf = (C - (1 + b) * mf) / (a * mf)
        ok = Pf > 0.02
        ax.plot(mf[ok], Pf[ok], color=nv['color'], lw=1.4, ls=(0, (1, 1.5)), zorder=3)
        ax.plot(nv['cruce_t'], nv['cruce_MW_t'], 'D', ms=6, color=nv['color'], mec='#111111', mew=0.8, zorder=6)
        ax.annotate(f"{1000 * nv['cruce_t']:.0f} kg", xy=(nv['cruce_t'], nv['cruce_MW_t']), xytext=(6, 2),
                    textcoords='offset points', fontsize=7.5, color='#111111', fontweight='bold', va='bottom', zorder=7)
    ax.set_xlim(0.002, 60)
    ax.set_ylim(0.02, 2e4)
    ax.xaxis.set_major_locator(FixedLocator([0.002, 0.01, 0.1, 1, 10, 60]))
    ax.set_xticklabels(['0.002', '0.01', '0.1', '1', '10', '60'])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1000, 1e4]))
    ax.set_yticklabels(['0.1', '1', '10', '100', '1000', '10$^4$'])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('Release (t)')
    ax.set_ylabel('Heat-rejection power per tonne released (MW/t)')
    ax.grid(True, which='major', alpha=0.25)
    # leyendas abajo: las naves en una fila, los estilos debajo
    fig.legend(handles=[Line2D([], [], color=F['naves'][c]['color'], lw=1.8, label=F['naves'][c]['nombre'])
                        for c in orden], loc='lower center', bbox_to_anchor=(0.5, 0.1), ncol=4, frameon=False)
    fig.legend(handles=[Line2D([], [], color='#444444', lw=1.8, label='no hold-down device (ideal bound), ~ m$^{0.5}$; '
                                                                       'dots: 1 t'),
                        Line2D([], [], color='#444444', lw=1.4, ls=(0, (1, 1.5)), label='largest that fits the payload'),
                        Line2D([], [], ls='none', marker='D', ms=6, color='#999999', mec='#111111', mew=0.8,
                               label='largest release that fits without a hold-down device')],
               loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False)
    plot_style.normalizar_paper(fig)
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    fig.savefig(rutas.figura('exp22_potencia_por_tonelada'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def figura_entrega(d):
    """La entrega de referencia en el tiempo (parte E), en paneles apilados que comparten el tiempo (no hay doble eje):
    arriba, el H2 de los gases; abajo, la energía neta acumulada, lo que el motor entrega a la planta (ya descontados
    ventiladores, compresor del escape y hélices) menos la regeneración, entre el caso pesimista y el optimista. A la
    izquierda la entrega (horas, sueltas una tras otra); a la derecha la regeneración en la base con el PEM del equipo
    (días). El regreso a la base no se modela."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style
    from matplotlib.lines import Line2D
    E, tt = d['E'], d['E']['total']
    n, T = E['n'], E['t_paso_s']
    th = np.arange(n + 1) * T / 3600.0                                          # [h] fin de cada suelta
    gas = np.array(E['gas_kg']) / 1e3                                            # [t]
    dias = E['t_regen_pem_h'] / 24.0
    td = np.array([0.0, dias])
    neto, regen = [x / 1e3 for x in tt['neto_kWh']], [x / 1e3 for x in tt['regen_kWh']]   # [MWh]
    ent = {'pes': np.array([0.0, neto[0]]), 'opt': np.array([0.0, neto[1]])}
    reg = {'pes': np.array([neto[0], neto[0] - regen[1]]), 'opt': np.array([neto[1], neto[1] - regen[0]])}
    AZUL, NEGRO, GRIS = '#4169E1', '#2b2926', '#8a8882'
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.4), sharex='col', sharey='row',
                            gridspec_kw={'width_ratios': (1, 1), 'wspace': 0.05, 'hspace': 0.08})
    (ag, bg), (ae, be) = axs
    # ── arriba: el H2 de los gases ──
    ag.plot(th, gas, color=AZUL, lw=1.8, zorder=3)
    ag.plot(th[:-1], gas[:-1], 'o', ms=3.2, color=AZUL, mec='white', mew=0.5, zorder=4)
    bg.plot(td, [gas[-1], gas[-1] + tt['H2_regen_kg'] / 1e3], color=AZUL, lw=1.8, zorder=3)
    for ax in (ag, bg):
        ax.axhline(gas[0], color=GRIS, lw=0.9, ls=(0, (4, 3)), zorder=1)
    pct = 100.0 * (gas[0] - gas[-1]) / gas[0]
    ag.annotate(f"departure: {gas[0]:.2f} t, with {n * E['s_t']:.0f} t of cargo\nand {E['equipo_t']:.1f} t of equipment",
                xy=(th[0], gas[0]), xytext=(6, 6), textcoords='offset points', fontsize=7.5, va='bottom')
    ag.text(0.03, 0.05, f"at the end: {gas[-1]:.2f} t; {tt['H2_kg'] / 1e3:.2f} t burned ({pct:.0f} %),\n"
            f"{tt['agua_kg'] / 1e3:.1f} t of water kept as ballast", transform=ag.transAxes, fontsize=7.5, va='bottom')
    bg.text(0.03, 0.72, f"regenerated: {tt['H2_regen_kg'] / 1e3:.2f} t", transform=bg.transAxes, fontsize=7.5,
            va='top')
    bg.text(0.97, 0.05, 'after the return to base (not modeled)', transform=bg.transAxes, fontsize=7.5, color=GRIS,
            ha='right', va='bottom')
    ag.set_ylabel('Hydrogen in the gas cells (t)')
    ag.set_ylim(gas[-1] - 0.3, gas[0] + 0.45)
    # ── abajo: la energía neta acumulada ──
    for ax, x, yy in ((ae, th[[0, -1]], ent), (be, td, reg)):
        ax.fill_between(x, yy['pes'], yy['opt'], color=NEGRO, alpha=0.18, lw=0, zorder=2)
        ax.plot(x, yy['opt'], color=NEGRO, lw=1.3, zorder=3)
        ax.plot(x, yy['pes'], color=NEGRO, lw=1.3, ls=(0, (4, 2)), zorder=3)
        ax.axhline(0.0, color=GRIS, lw=0.9, zorder=1)
    ae.annotate(f"engine, net to the plant:\n+{neto[0]:.1f} to +{neto[1]:.1f} MWh", xy=(th[-1], neto[1]), xytext=(-4, 6),
                textcoords='offset points', fontsize=7.5, ha='right', va='bottom')
    be.text(0.97, 0.95, f"regeneration (electrolysis\nand compression): {regen[0]:.1f}–{regen[1]:.1f} MWh",
            transform=be.transAxes, fontsize=7.5, ha='right', va='top')
    be.annotate(f"net from the ground:\n{-reg['opt'][1]:.0f}–{-reg['pes'][1]:.0f} MWh", xy=(dias, reg['pes'][1]),
                xytext=(-4, -4), textcoords='offset points', fontsize=7.5, ha='right', va='top')
    ae.set_ylabel('Net energy, cumulative (MWh)')
    ae.set_ylim(1.25 * reg['pes'][1], 0.5 * regen[1])
    ae.set_xlim(0, th[-1])
    be.set_xlim(0, dias)
    ae.set_xlabel(f"Delivery: one {E['s_t']:g}-t release every {T / 60:.0f} min (h)")
    be.set_xlabel(f"At base: regeneration, onboard stack of {E['P_pem_kW']:.0f} kW (days)")
    # eje cortado: sin los bordes interiores, con marcas de corte
    for ax in (ag, ae):
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='y', right=False)
    for ax in (bg, be):
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', left=False)
    k = dict(color='#444444', lw=0.9, clip_on=False)
    for ax, x0 in ((ag, 1.0), (ae, 1.0), (bg, 0.0), (be, 0.0)):
        for y0 in (0.0, 1.0):
            ax.plot([x0 - 0.012, x0 + 0.012], [y0 - 0.025, y0 + 0.025], transform=ax.transAxes, **k)
    for ax in axs.flat:
        ax.grid(True, alpha=0.25)
    etas, fs, A = d['d21']['etas'], d['d21']['fs'], p.freno_A_discos
    fig.legend(handles=[Line2D([], [], color=AZUL, lw=1.8, marker='o', ms=3.2, mec='white', mew=0.5,
                               label=f"hydrogen in the gas cells; dots: each {E['s_t']:g}-t release"),
                        Line2D([], [], color=GRIS, lw=0.9, ls=(0, (4, 3)), label='hydrogen at departure'),
                        Line2D([], [], color=NEGRO, lw=1.3, label=f"optimistic: engine efficiency {max(etas):.2f}, fans "
                               f"{100 * min(fs):.1f} % of the heat, {max(A):.0f} m² disks"),
                        Line2D([], [], color=NEGRO, lw=1.3, ls=(0, (4, 2)), label=f"pessimistic: {min(etas):.2f}, "
                               f"{100 * max(fs):.0f} %, {min(A):.0f} m²")],
               loc='lower center', ncol=2, frameon=False)
    plot_style.normalizar_paper(fig)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.97, bottom=0.2)
    fig.savefig(rutas.figura('exp22_entrega_tiempo'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def informe(d, ch):
    h = lambda s: (f'{s / 60:.1f} min' if s < 3600 else f'{s / 3600:.1f} h') if s else '—'
    n = d['neto_kWh']
    L = ['EXP-22 — ENERGÍA Y POTENCIA DE LAS DOS ETAPAS, POR TONELADA (modelo/exp22_energia_potencia.py)', '',
         'A. CICLO COMPLETO POR TONELADA (kWh)  [suelta: EXP-21; regeneración: electrólisis + compresión]']
    for fl in d['d21']['filas']:
        k = fl['kWh']
        ext = [k['electrolisis'] + k['compresion_H2'][0] - k['neto_planta'][1],
               k['electrolisis'] + k['compresion_H2'][1] - k['neto_planta'][0]]
        quim = fl['masa_kg']['H2_regenerado'] * p.Q_comb_H2 / KWH
        L.append(f"  {fl['nombre']:20s} quema {k['quimica_PCI']:.0f}, neto {k['neto_planta'][0]:.0f}-{k['neto_planta'][1]:.0f}"
                 f" | regenera: electrólisis {k['electrolisis']:.0f} + compresión {k['compresion_H2'][0]:.0f}-"
                 f"{k['compresion_H2'][1]:.0f}, devuelve {quim:.0f} | externa neta {ext[0]:.0f}-{ext[1]:.0f}")
    L += ['', f"B. SUELTA CON FRENO, LCA60T, 1 t, 14 atm  (neto del motor sin hélices: {n[0]:.0f}-{n[1]:.0f} kWh por tonelada)",
          f"  límite de POTENCIA (hélices al soltar = potencia neta del motor): hasta {h(d['limites_s']['potencia_pes'])}"
          f" (pesimista) / {h(d['limites_s']['potencia_opt'])} (optimista)",
          f"  límite de ENERGÍA (hélices en toda la recuperación = neto del motor): hasta {h(d['limites_s']['energia_pes'])}"
          f" (pesimista) / {h(d['limites_s']['energia_opt'])} (optimista)"]
    for P in (2e6, 10e6, 39e6):
        b = min(d['B'], key=lambda x: abs(x['P_W'] - P))
        L.append(f"  {b['P_W'] / 1e6:5.1f} MW: {h(b['t_s'])}, la quema absorbe {b['m_e_kg']:.0f} kg; hélices "
                 f"{min(b['P_hel_W']) / 1e6:.2f}-{max(b['P_hel_W']) / 1e6:.2f} MW al soltar, "
                 f"{min(b['E_hel_J']) / KWH:.0f}-{max(b['E_hel_J']) / KWH:.0f} kWh en la recuperación")
    L += ['', 'C. TAMAÑO DE LA SUELTA, LCA60T: energía de las hélices por tonelada (kWh/t)']
    for c in d['C']:
        L.append(f"  {c['m_t']:5g} t | con 10 MW: {h(c['t10_s'])}, {c['E10_kWh_t'][0]:.0f}-{c['E10_kWh_t'][1]:.0f} kWh/t"
                 f" | con 10 min por tonelada: {c['ET_kWh_t'][0]:.0f}-{c['ET_kWh_t'][1]:.0f} kWh/t, hélices "
                 f"{min(c['PT_hel_W']) / 1e6:.2f}-{max(c['PT_hel_W']) / 1e6:.2f} MW")
    L += ['', f"D. REGENERACIÓN: {d['E_el_kWh']:.0f} kWh de electrólisis por tonelada"]
    for x in d['D']:
        L.append(f"  en {x['t_h']:5g} h: electrolizador {x['P_el_kW']:.0f} kW por tonelada, PEM {x['masa_PEM_t'][0]:.2f}-"
                 f"{x['masa_PEM_t'][1]:.2f} t")
    E, ps, tt = d['E'], d['E']['paso'], d['E']['total']
    r2 = lambda x, dec=0, f=1.0: f"{x[0] / f:.{dec}f}-{x[1] / f:.{dec}f}"
    L += ['', f"E. ENTREGA DE REFERENCIA: {E['nave']}, {E['n']} sueltas de {E['s_t']:g} t, {E['t_paso_s'] / 60:.0f} min "
              f"cada una (carga de negocio {E['N_t']:.2f} t, equipo {E['equipo_t']:.2f} t; [pesimista, optimista])",
          f"  {'':44s} {'por suelta':>18s} {'entrega':>18s}",
          f"  {'duración':44s} {E['t_paso_s'] / 60:15.0f} min {E['n'] * E['t_paso_s'] / 3600:16.1f} h",
          f"  {'rechazo instalado [MW]':44s} {ps['P_rech_W'] / 1e6:18.2f}",
          f"  {'calor rechazado con el motor [kWh | MWh]':44s} {r2(ps['calor_kWh']):>18s} {r2(tt['calor_kWh'], 1, 1e3):>18s}",
          f"  {'H2 quemado [kg | t]':44s} {ps['H2_kg']:18.1f} {tt['H2_kg'] / 1e3:18.2f}",
          f"  {'agua retenida como lastre [kg | t]':44s} {ps['agua_kg']:18.0f} {tt['agua_kg'] / 1e3:18.1f}",
          f"  {'hélices, potencia al soltar [MW]':44s} {r2(ps['P_hel_W'], 2, 1e6):>18s}",
          f"  {'hélices, energía [kWh | MWh]':44s} {r2(ps['hel_kWh'], 1):>18s} {r2(tt['hel_kWh'], 2, 1e3):>18s}",
          f"  {'neto del motor sin hélices [kWh | MWh]':44s} {r2(ps['motor_kWh']):>18s} {r2(tt['motor_kWh'], 1, 1e3):>18s}",
          f"  {'neto a la planta [kWh | MWh]':44s} {r2(ps['neto_kWh']):>18s} {r2(tt['neto_kWh'], 1, 1e3):>18s}",
          f"  {'regeneración: electrólisis + compresión':44s} {r2(ps['regen_kWh']):>18s} {r2(tt['regen_kWh'], 1, 1e3):>18s}",
          f"  {'electricidad de tierra neta':44s} {r2(ps['externa_kWh']):>18s} {r2(tt['externa_kWh'], 1, 1e3):>18s}",
          f"  regeneración en la base con el PEM del equipo ({E['P_pem_kW']:.0f} kW, una suelta en "
          f"{p.t_regeneracion_h:g} h): {E['t_regen_pem_h'] / 24:.0f} días; toda la entrega en {p.t_regeneracion_h:g} h"
          f" pide {r2(E['P_base_MW'], 2)} MW de electrólisis en tierra",
          f"  H2 de los gases: {E['gas_kg'][0]:.0f} kg al salir ({E['carga0_t']:.2f} t a bordo), {E['gas_kg'][-1]:.0f} kg"
          f" al terminar: se quema el {100 * (E['gas_kg'][0] - E['gas_kg'][-1]) / E['gas_kg'][0]:.1f} %",
          '', 'F. POTENCIA DE RECHAZO POR TONELADA (MW/t); la suelta más grande que cabe sin freno (cruce)']
    for c, nv in d['F']['naves'].items():
        L.append(f"  {nv['nombre']:16s} sin freno, 1 t (o la carga útil): {nv['P_1t_MW_t']:6.0f} MW/t | cruce: "
                 f"{1000 * nv['cruce_t']:5.0f} kg a {nv['cruce_MW_t']:4.0f} MW/t ({nv['cruce_t'] * nv['cruce_MW_t']:.1f} MW)")
    L += ['', 'CHEQUEOS'] + [f"  {'ok  ' if ok else 'FALLA'} {n_}: {a:.6g} vs {b:.6g}" for n_, a, b, ok in ch] + \
         [f"{sum(ok for *_, ok in ch)}/{len(ch)} chequeos ok"]
    return '\n'.join(L)


# ── Salidas LaTeX para la bitácora (bitacora/EXP22_energia_potencia.tex): macros \valVeintidos... y tablas en
# resultados/exploracion/, para que la bitácora no lleve cifras a mano. Fuera de numeros.tex y del sello.
SALIDA = os.path.join(rutas.RESULTADOS, 'exploracion')


def _f(x, dec):
    s = f'{x:,.{dec}f}' if abs(x) >= 10000 else f'{x:.{dec}f}'
    return s.replace(',', '{,}').replace('-', '$-$') if x < 0 else s.replace(',', '{,}')


def _r(a, b, dec):
    ta, tb = _f(min(a, b), dec), _f(max(a, b), dec)
    return ta if ta == tb else f'{ta}--{tb}'


def _dur(s):
    return f'{s:.0f}~s' if s < 120 else (f'{s / 60:.1f}~min' if s < 5400 else f'{s / 3600:.1f}~h')


def _tabla(nombre, columnas, cabecera, filas):
    texto = ('\\begin{tabular}{' + columnas + '}\n\\toprule\n' + cabecera + ' \\\\\n\\midrule\n'
             + ''.join((f + '\n') if f.startswith(('\\midrule', '\\addlinespace')) else (f + ' \\\\\n') for f in filas)
             + '\\bottomrule\n\\end{tabular}\n')
    with open(os.path.join(SALIDA, f'exp22_tabla_{nombre}.tex'), 'w', encoding='utf-8') as fh:
        fh.write('% Generado por exp22_energia_potencia.py (EXPLORATORIO, fuera del sello). No editar a mano.\n' + texto)


def exportar_tex(d):
    os.makedirs(SALIDA, exist_ok=True)
    M = []

    def val(nombre, texto, comentario):
        M.append(f'\\newcommand{{\\valVeintidos{nombre}}}{{{texto}}}  % {comentario}')

    d21, E, ps, tt, F = d['d21'], d['E'], d['E']['paso'], d['E']['total'], d['F']
    filas21 = d21['filas']
    k14, m14 = filas21[1]['kWh'], filas21[1]['masa_kg']
    r20 = rutas.leer_datos('exp20_reparto_sueltas')
    # ── parámetros y por tonelada ──
    val('Eta', _r(*d21['etas'], 2), 'eficiencia del motor, pesimista--optimista')
    val('F', _r(*d21['fs'], 3), 'fracción del calor que piden los ventiladores, optimista--pesimista')
    val('FM', _f(p.freno_FM, 1), 'figura de mérito de las hélices')
    val('A', _r(*p.freno_A_discos, 0), 'área total de discos de las hélices [m2]')
    val('Rho', _f(p.rho_aire_std, 3), 'densidad del aire [kg/m3]')
    val('PCI', _f(p.Q_comb_H2 / 1e6, 1), 'PCI del H2 [MJ/kg]')
    val('Phi', _f(p.phi_objetivo, 2), 'fracción del agua retenida a 14 atm')
    val('EEl', _f(p.E_H2_electrolisis, 1), 'electricidad de electrólisis [kWh/kg]')
    val('EComp', _r(*d21['e_comp_kWhkg'], 2), 'electricidad de compresión del H2 regenerado [kWh/kg]')
    val('Wc', _f(d21['w_comp_Jkg'] / 1e6, 2), 'trabajo isotérmico mínimo del compresor del escape [MJ/kg H2]')
    val('Qcero', _f(filas21[1]['q0_Jkg'] / 1e6, 1), 'calor a rechazar por kg de H2 sin motor, q0 [MJ/kg]')
    val('QGJ', _f(d['Q_J'] / 1e9, 2), 'calor por tonelada, Q [GJ]')
    val('QkWh', _f(d['Q_J'] / KWH, 0), 'calor por tonelada, Q [kWh]')
    val('Hdos', _f(m14['H2_quemado'], 1), 'H2 quemado por tonelada [kg]')
    val('Agua', _f(m14['agua_retenida'], 0), 'agua retenida por tonelada [kg]')
    val('Quim', _f(k14['quimica_PCI'], 0), 'energía química quemada por tonelada [kWh]')
    val('NetoMotor', _r(*d['neto_kWh'], 0), 'neto del motor por tonelada, sin hélices [kWh]')
    val('NetoHD', _r(*k14['neto_planta'], 0), 'neto a la planta por tonelada, HD 14 atm 10 MW [kWh]')
    val('ElT', _f(k14['electrolisis'], 0), 'electrólisis por tonelada, 14 atm [kWh]')
    val('CompT', _r(*k14['compresion_H2'], 0), 'compresión del H2 regenerado por tonelada, 14 atm [kWh]')
    val('RegenT', _r(*[k14['electrolisis'] + c for c in k14['compresion_H2']], 0), 'regeneración por tonelada [kWh]')
    ext = [k14['electrolisis'] + k14['compresion_H2'][0] - k14['neto_planta'][1],
           k14['electrolisis'] + k14['compresion_H2'][1] - k14['neto_planta'][0]]
    val('ExtT', _r(*ext, 0), 'electricidad de tierra neta por tonelada, HD 14 atm [kWh]')
    val('PagaPct', _r(100 * k14['neto_planta'][0] / (k14['electrolisis'] + k14['compresion_H2'][1]),
                      100 * k14['neto_planta'][1] / (k14['electrolisis'] + k14['compresion_H2'][0]), 0),
        'parte de la regeneración que paga el neto del motor, HD 14 atm [%]')
    val('TRegen', _f(p.t_regeneracion_h, 0), 'tiempo de regeneración de una tonelada con el PEM del equipo [h]')
    quim14 = m14['H2_regenerado'] * p.Q_comb_H2 / KWH
    val('PerdElPct', _f(100 * (k14['electrolisis'] - quim14) / k14['electrolisis'], 0),
        'parte de la electricidad de la electrólisis que sale como calor [%]')
    k1, m1 = filas21[2]['kWh'], filas21[2]['masa_kg']
    ext1 = [k1['electrolisis'] + k1['compresion_H2'][0] - k1['neto_planta'][1],
            k1['electrolisis'] + k1['compresion_H2'][1] - k1['neto_planta'][0]]
    val('ExtUna', _r(*ext1, 0), 'electricidad de tierra neta por tonelada, HD 1 atm [kWh]')
    # (2026-10-12, autor) por ahora el análisis no considera la reposición del H2 del vapor no retenido
    val('PhiUna', _f(filas21[2]['phi'], 3), 'fracción del agua retenida a 1 atm (phi_max)')
    val('RegenHdos', _f(m14['H2_regenerado'], 1), 'H2 regenerado por tonelada, 14 atm [kg]')
    val('RegenHdosUna', _f(m1['H2_regenerado'], 1), 'H2 regenerado por tonelada, 1 atm [kg]')
    val('HelPct', _r(100 * min(d['E']['paso']['hel_kWh']) / max(d['neto_kWh']),
                     100 * max(d['E']['paso']['hel_kWh']) / min(d['neto_kWh']), 0),
        'energía de las hélices / neto del motor, sueltas de 1 t de 10 min [%]')
    pem24 = next(x for x in d['D'] if abs(x['t_h'] - p.t_regeneracion_h) < 1e-9)['masa_PEM_t']
    val('PemDiaT', _r(*[d['E']['n'] * x for x in pem24], 1), 'PEM para regenerar toda la entrega en un día [t]')
    # ── B ──
    lim = d['limites_s']
    val('LimPotPes', _f(lim['potencia_pes'] / 60, 0), 'hélices al soltar <= potencia neta del motor hasta, pesimista [min]')
    val('LimPotOpt', _f(lim['potencia_opt'] / 3600, 1), 'lo mismo, optimista [h]')
    val('LimEnePes', _f(lim['energia_pes'] / 3600, 1), 'energía de las hélices <= neto del motor hasta, pesimista [h]')
    val('LimEneOpt', _f(lim['energia_opt'] / 3600, 1), 'lo mismo, optimista [h]')
    val('PBarrido', f"{d['B'][0]['P_W'] / 1e6:g}--{d['B'][-1]['P_W'] / 1e6:g}", 'barrido del rechazo instalado [MW]')
    val('TBarrido', f"{_dur(d['B'][-1]['t_s'])}--{_dur(d['B'][0]['t_s'])}", 'barrido del tiempo de recuperación')
    val('Hdiez', _r(*d['B_tabla'][1]['E_hel_kWh'], 0), 'energía de las hélices, 1 t en 10 min [kWh]')
    val('CompEsc', _f(k14['compresor_escape'], 0), 'compresor del escape por tonelada, 14 atm [kWh]')
    val('EtaEl', _f(p.eta_electrolisis, 2), 'eficiencia de la electrólisis sobre el PCI')
    val('NetoPaso', _r(*d['E']['paso']['neto_kWh'], 0), 'neto a la planta, 1 t en 10 min (crédito de la suelta) [kWh]')
    val('QuimRegen', _f(quim14, 0), 'electrólisis: vuelve como H2 a los tanques, por tonelada [kWh]')
    val('PerdEl', _f(k14['electrolisis'] - quim14, 0), 'electrólisis: sale como calor, por tonelada [kWh]')
    pem24t = next(x for x in d['D'] if abs(x['t_h'] - p.t_regeneracion_h) < 1e-9)['masa_PEM_t']
    val('PemMasa', _r(*pem24t, 2), 'PEM para una tonelada en t_regeneracion_h [t]')
    val('PRegen', _r(*[(k14['electrolisis'] + c) / p.t_regeneracion_h for c in k14['compresion_H2']], 0),
        'electrolizador + compresor en t_regeneracion_h, por tonelada [kW]')
    perdR = [k14['electrolisis'] + c - quim14 for c in k14['compresion_H2']]
    val('PerdRegen', _r(*perdR, 0), 'regeneración: lo que no queda como H2 (calor), por tonelada [kWh]')
    val('GuardaPct', _r(*[100 * quim14 / (k14['electrolisis'] + c) for c in k14['compresion_H2']], 0),
        'regeneración: parte de la electricidad que queda como H2 (PCI) [%]')
    val('CalorRegenKW', _r(*[x / p.t_regeneracion_h for x in perdR], 0),
        'calor a sacar del PEM y del compresor con t_regeneracion_h por tonelada [kW]')
    # G. el resto de la regeneración con sol (paneles del casco del LCA60T, Lima)
    Gs = d['G']
    todos = lambda f: [v for c in Gs['resto_kWh'] for v in f(c)]
    val('SolA', _f(Gs['A_m2'], 0), 'área de paneles del casco, LCA60T [m2]')
    val('SolGInv', _f(Gs['G_kWh_m2_dia']['invierno'], 1), 'irradiancia efectiva sobre el casco, invierno [kWh/m2/día]')
    val('SolGVer', _f(Gs['G_kWh_m2_dia']['verano'], 1), 'lo mismo, verano [kWh/m2/día]')
    val('Resto', _r(min(todos(lambda c: Gs['resto_kWh'][c])), max(todos(lambda c: Gs['resto_kWh'][c])), 0),
        'regeneración no cubierta por el neto de la suelta almacenado, hélices o tether [kWh/t]')
    for t_, suf in (('a-Si', 'ASi'), ('CIGS', 'Cigs')):
        x = Gs['tec'][t_]
        val('SolP' + suf, _f(x['P_pleno_kW'], 0), f'paneles {t_} del casco a sol pleno [kW]')
        val('SolEf' + suf, _f(100 * x['ef'], 1), f'eficiencia de los paneles {t_} [%]')
        val('SolHoras' + suf, _r(min(todos(lambda c: x['horas'][c])), max(todos(lambda c: x['horas'][c])), 1),
            f'horas de sol pleno para el resto, paneles {t_} [h]')
        for est, se in (('verano', 'Ver'), ('invierno', 'Inv')):
            val('SolDias' + suf + se, _r(min(todos(lambda c: x['dias'][c][est])), max(todos(lambda c: x['dias'][c][est])),
                                         2), f'días de Lima para el resto, paneles {t_}, {est}')
    filas_s = []
    for c, nomc in (('hel', 'Hélices, 10~min'), ('amarre', 'Tether')):
        for i, caso in enumerate(('pesimista', 'optimista')):
            celdas = [f"{nomc}, {caso}", _f(Gs['resto_kWh'][c][i], 0)]
            for t_ in ('a-Si', 'CIGS'):
                x = Gs['tec'][t_]
                celdas += [_f(x['horas'][c][i], 1), _f(x['dias'][c]['verano'][i], 2), _f(x['dias'][c]['invierno'][i], 2)]
            filas_s.append(' & '.join(celdas))
    _tabla('solar', 'l r r r r r r r',
           ' & Resto (kWh) & \\multicolumn{3}{c}{a-Si} & \\multicolumn{3}{c}{CIGS} \\\\\n'
           ' & & sol pleno (h) & verano (d) & invierno (d) & sol pleno (h) & verano (d) & invierno (d)', filas_s)
    # H. carga por viaje contra tiempo de entrega
    H = d['H']
    val('TechoUnica', _f(H['techo_t'], 1), 'suelta única más grande, compresor en unidades, recuperación infinita [t]')
    val('TechoPaper', _f(H['techo_paper_t'], 1), 'lo mismo en la carta de diseño (EXP-18; igual desde 2026-10-13) [t]')
    val('Rmax', _f(H['C_t'] / H['aQ_h'], 1), 'ritmo de entrega al que el rechazo se come la carga útil, C / aQ [t/h]')
    # I. el módulo de regeneración, independiente de la nave
    I = d['I']
    val('ModKW', _r(*I['kW_por_kgph'], 1), 'módulo a tope: potencia por kg de H2 por hora (PEM + compresor) [kW]')
    val('ModPemKW', _f(I['e_el_kWh_kg'], 1), 'PEM: potencia por kg de H2 por hora [kW]')
    val('ModCompKW', _r(*I['e_comp_kWh_kg'], 2), 'compresor de H2: potencia por kg de H2 por hora [kW]')
    val('ModPemKg', _r(*I['pem_kg_por_kgph'], 0), 'PEM: masa por kg de H2 por hora [kg]')
    val('ModCompCap', _f(I['cap_comp_kgph'], 2), 'una unidad de compresor de catálogo [kg de H2 por hora]')
    val('ModCompKg', _f(I['comp_kg_unidad'], 0), 'masa de esa unidad [kg]')
    for P in ('350', '700'):
        x = I['tanques'][P]
        val('Tanque' + {'350': 'Tres', '700': 'Siete'}[P] + 'Kg', _f(x['kg_por_kg'], 1),
            f'tanques a {P} bar: masa por kg de H2 guardado [kg/kg]')
        val('Tanque' + {'350': 'Tres', '700': 'Siete'}[P] + 'L', _f(x['L_por_kg'], 0),
            f'tanques a {P} bar: volumen por kg de H2 guardado [L/kg]')
    # I2. el módulo alimentado solo con paneles, en % de la carga útil
    J = d['I2']
    val('GAnual', _f(J['G_h_dia'], 1), 'horas de sol pleno por día, Lima, promedio anual [h]')
    val('TanqueTonPorT', _f(J['tanque_t_por_t'], 2), 'tanque vacío para el H2 de una tonelada liberada [t por t]')
    for tec, suf in (('a-Si', 'ASi'), ('CIGS', 'Cigs')):
        val('CapCasco' + suf, _f(J['cap_casco_kgph'][tec], 0), f'capacidad que dan los paneles {tec} del casco del LCA60T [kg/h]')
        m = J['tabla'][tec]['12|7']
        val('Ejemplo' + suf, _r(*m['total_t'], 1), f'12 kg/h por una semana con {tec}: masa total [t]')
        val('EjemploPct' + suf, _r(*[100 * x / J['C_t'] for x in m['total_t']], 0), f'lo mismo en % de 60 t')
    # I4. costo por tonelada recuperada a bordo
    C4 = d['I4']
    val('CostoPiso', _f(C4['piso_t_por_t'], 2), 'piso del costo por tonelada recuperada: los tanques [t por t]')
    nomf = {'a-Si': 'ASi', 'CIGS': 'Cigs', 'red': 'Red'}
    for f, x in C4['fuentes'].items():
        val('CostoA' + nomf[f], _r(*x['a'], 1 if x['a'][1] > 1 else 2), f'a del costo por tonelada, {f} [t día por t]')
        val('CostoEquil' + nomf[f], _r(*x['T_equilibrio_d'], 0 if x['T_equilibrio_d'][0] > 5 else 1),
            f'desde cuántos días cuesta menos de 1 t por t, {f}')
        val('CostoSemana' + nomf[f], _r(*x['costo']['7'], 2), f'costo por tonelada en una semana, {f} [t por t]')
    filas4 = []
    for f, nom in (('red', 'Red de tierra (24~h al día)'), ('a-Si', 'Paneles a-Si (Lima, anual)'),
                   ('CIGS', 'Paneles CIGS (Lima, anual)')):
        x = C4['fuentes'][f]
        dec = 1 if x['a'][1] > 1 else 2
        filas4.append(f"{nom} & {_r(*x['a'], dec)} & " + ' & '.join(_r(*x['costo'][f'{T:g}'], 2) for T in (1, 2, 7, 30))
                      + f" & {_r(*x['T_equilibrio_d'], 0 if x['T_equilibrio_d'][0] > 5 else 1)}")
    _tabla('costo', 'l r r r r r r',
           'Fuente & $a$ (t·día/t) & 1~día & 2~días & 1~semana & 30~días & $<$1~t/t desde (días)', filas4)
    # I2 universal
    U = d['I2U']
    for t_, suf in (('a-Si', 'ASi'), ('CIGS', 'Cigs')):
        x = U['tec'][t_]
        val('UnivK' + suf, _f(x['K_kg_por_kgph'], 0), f'parte de potencia por kg/h, {t_}, PEM pesado [kg]')
        val('UnivX' + suf, _f(x['x_star_h'], 0), f'x* = K / tanque por kg, {t_} [h de sol pleno]')
        val('UnivXDias' + suf, _f(x['x_star_dias_anual'], 1 if x['x_star_dias_anual'] < 10 else 0),
            f'lo mismo en días de Lima (promedio anual), {t_}')
    val('UnivKComp', _f(U['k_comp'], 0), 'compresor continuo: masa de la unidad / su caudal [kg por kg/h]')
    # I2, lectura inversa: presupuesto y días de viaje -> capacidad que entra
    K = d['I3']
    val('AguaPorKg', _f(K['agua_por_kg'], 2), 'agua electrolizada por kg de H2 regenerado [kg]')
    val('OdosPorKg', _f(K['O2_por_kg'], 2), 'oxígeno venteado por kg de H2 regenerado [kg]')
    m_h2_t = d['d21']['filas'][1]['masa_kg']['H2_regenerado']
    for f in K['filas']:
        if f['presupuesto_t'] == 6.0 and f['dias'] == 2.0:
            suf = 'ASi' if f['tec'] == 'a-Si' else 'Cigs'
            val('InvCap' + suf, _f(f['c_kgph'], 1), f"6 t y 2 días de viaje, {f['tec']}: capacidad que entra [kg/h]")
            val('InvHdos' + suf, _f(f['H2_kg'], 0), f"lo mismo: H2 regenerado [kg]")
            val('InvSuelta' + suf, _f(f['H2_kg'] / m_h2_t, 2), f"lo mismo: toneladas de suelta recuperadas [t]")
    filas_k = []
    for tec in ('a-Si', 'CIGS'):
        for f in [x for x in K['filas'] if x['tec'] == tec]:
            filas_k.append(f"{tec} & {f['presupuesto_t']:g} & {f['dias']:g} & {_f(f['c_kgph'], 1)} & "
                           f"{_f(f['paneles_t'], 2)} & {_f(f['pem_t'], 2)} & {_f(f['comp_t'], 2)} & "
                           f"{_f(f['tanques_t'], 2)} & {_f(f['total_t'], 2)} & {_f(f['H2_kg'], 0)} & "
                           f"{_f(f['H2_kg'] / m_h2_t, 2)}")
            if f['dias'] == 7.0 and not (tec == 'CIGS' and f['presupuesto_t'] == 20.0):
                filas_k.append('\\midrule' if f['presupuesto_t'] == 20.0 else '\\addlinespace')
    _tabla('inversa', 'l r r r r r r r r r r',
           'Paneles & Presupuesto (t) & Días & Capacidad (kg/h) & Paneles (t) & PEM (t) & Compresor (t) & Tanques (t) & '
           'Total (t) & H$_2$ (kg) & Suelta (t)', filas_k)
    _tabla('modulo', 'l r r',
           ' & Por kg de H$_2$ por hora & Fuente', [
               f"Potencia del PEM & {_f(I['e_el_kWh_kg'], 1)}~kW & electrólisis, {_f(p.E_H2_electrolisis, 1)}~kWh/kg",
               f"Potencia del compresor de H$_2$ & {_r(*I['e_comp_kWh_kg'], 2)}~kW & 20 a 350~bar, DOE típica--medida",
               f"\\textbf{{Potencia del módulo}} & \\textbf{{{_r(*I['kW_por_kgph'], 1)}~kW}} & ",
               f"Masa del PEM & {_r(*I['pem_kg_por_kgph'], 0)}~kg & {_r(*p.kWkg_PEM, 1)}~kW/kg",
               f"Masa del compresor & {_f(I['comp_kg_unidad'], 0)}~kg por cada {_f(I['cap_comp_kgph'], 2)}~kg/h & "
               f"unidad de catálogo",
               '\\midrule',
               f"Tanques a 350~bar & {_f(I['tanques']['350']['kg_por_kg'], 1)}~kg y {_f(I['tanques']['350']['L_por_kg'], 0)}~L "
               f"por kg guardado & Tipo IV del DOE, escalado",
               f"Tanques a 700~bar & {_f(I['tanques']['700']['kg_por_kg'], 1)}~kg y {_f(I['tanques']['700']['L_por_kg'], 0)}~L "
               f"por kg guardado & Tipo IV del DOE (2023)"])
    # el neto de la suelta, si se hubiera almacenado, se acaba a las E_neto / E_regen x t_regeneracion_h de la ventana
    E_reg = [k14['electrolisis'] + max(k14['compresion_H2']), k14['electrolisis'] + min(k14['compresion_H2'])]
    for nom, net in (('AcabaHel', d['E']['paso']['neto_kWh']), ('AcabaAmarre', d['neto_kWh'])):
        val(nom, _r(*[net[i] / E_reg[i] * p.t_regeneracion_h for i in (0, 1)], 2),
            f"neto de la suelta ({'hélices, 10 min' if nom == 'AcabaHel' else 'tether'}) almacenado: se acaba a las "
            '... h de la ventana de regeneración, pesimista--optimista [h]')
    val('Ventila', _r(*k14['ventiladores'], 0), 'ventiladores por tonelada, optimista--pesimista [kWh]')
    val('TrabMotor', _r(*k14['trabajo_motor'], 0), 'trabajo del motor por tonelada [kWh]')
    # ── EXP-20 (reparto) ──
    val('a', _f(r20['a_t_por_MW'], 2), 'masa de la unidad de rechazo por MW [t/MW]')
    val('b', _f(r20['b_t_por_t'], 2), 'cajas fijas por tonelada [t/t]')
    val('bp', _f(r20['b_sin_comp_t_por_t'], 2), "b' sin el compresor de H2 [t/t]")
    val('mc', _f(r20['comp_t'], 3), 'compresor de H2, una unidad de catálogo [t]')
    val('C', _f(E['carga_util_t'], 0), 'carga útil del LCA60T [t]')
    sin = filas21[0]                                                  # NHD, 14 atm: el ritmo sin freno, LCA60T 1 t
    val('PsinMW', _f(sin['P_W'] / 1e6, 0), 'rechazo sin freno, LCA60T 1 t [MW]')
    val('TsinS', _f(sin['t_s'], 0), 'duración de la quema sin freno, LCA60T 1 t [s]')
    val('aQ', _f(r20['aQ_s'] / 3600, 1), 'a Q [h]')
    # ── E ──
    val('Pasos', f"{E['n']}", 'sueltas de 1 t de 10 min en el LCA60T')
    val('N', _f(E['N_t'], 2), 'carga de negocio de la ecuación del reparto [t]')
    val('Equipo', _f(E['equipo_t'], 1), 'equipo de la entrega [t]')
    val('CargaCero', _f(E['carga0_t'], 1), 'a bordo al salir [t]')
    val('Horas', _f(E['n'] * E['t_paso_s'] / 3600, 0), 'duración de la entrega [h]')
    val('PrechMW', _f(ps['P_rech_W'] / 1e6, 1), 'rechazo instalado [MW]')
    val('HdosTotal', _f(tt['H2_kg'] / 1e3, 2), 'H2 quemado en la entrega [t]')
    val('AguaTotal', _f(tt['agua_kg'] / 1e3, 1), 'agua retenida en la entrega [t]')
    val('NetoTotal', _r(*[x / 1e3 for x in tt['neto_kWh']], 1), 'neto a la planta en la entrega [MWh]')
    val('RegenTotal', _r(*[x / 1e3 for x in tt['regen_kWh']], 1), 'regeneración de la entrega [MWh]')
    val('ExtTotal', _r(*[x / 1e3 for x in tt['externa_kWh']], 0), 'electricidad de tierra neta de la entrega [MWh]')
    val('PemKW', _f(E['P_pem_kW'], 0), 'PEM del equipo [kW]')
    val('Dias', _f(E['t_regen_pem_h'] / 24, 0), 'días de regeneración con el PEM del equipo')
    val('BaseMW', _r(*E['P_base_MW'], 1), 'electrólisis en tierra para regenerar la entrega en un día [MW]')
    val('GasIni', _f(E['gas_kg'][0] / 1e3, 2), 'H2 de los gases al salir [t]')
    val('GasFin', _f(E['gas_kg'][-1] / 1e3, 2), 'H2 de los gases al terminar [t]')
    val('GasPct', _f(100 * (E['gas_kg'][0] - E['gas_kg'][-1]) / E['gas_kg'][0], 1), 'H2 de los gases quemado [%]')
    val('HelPasoMW', _r(*[x / 1e6 for x in ps['P_hel_W']], 2), 'hélices al soltar, por suelta [MW]')
    val('HelPaso', _r(*ps['hel_kWh'], 1), 'energía de las hélices por suelta [kWh]')
    # ── F ──
    for c, nv in F['naves'].items():
        val('Cabe' + {'Hind': 'Hind', 'FW': 'FW', 'ZR': 'ZR', 'ZNT': 'NT'}[c], _f(1000 * nv['cruce_t'], 0),
            f"{nv['nombre']}: suelta más grande sin freno [kg]")
    with open(os.path.join(SALIDA, 'exp22_numeros.tex'), 'w', encoding='utf-8') as fh:
        fh.write('% Generado por exp22_energia_potencia.py (EXPLORATORIO, fuera del sello). No editar a mano.\n'
                 + '\n'.join(M) + '\n')

    # ── tabla A: balance por tonelada, tres arreglos ──
    def fila(nombre, fn):
        return nombre + ' & ' + ' & '.join(fn(f) for f in filas21)

    rr = lambda v, dec=0: _r(v[0], v[1], dec) if isinstance(v, list) else _f(v, dec)
    quim = lambda f: f['masa_kg']['H2_regenerado'] * p.Q_comb_H2 / KWH
    extf = lambda f: [f['kWh']['electrolisis'] + f['kWh']['compresion_H2'][0] - f['kWh']['neto_planta'][1],
                      f['kWh']['electrolisis'] + f['kWh']['compresion_H2'][1] - f['kWh']['neto_planta'][0]]
    _tabla('balance', 'l r r r', ' & ' + ' & '.join(f['nombre'].replace(', 10 MW', '') for f in filas21), [
        fila('$\\varphi$', lambda f: _f(f['phi'], 3)),
        fila('H$_2$ quemado (kg)', lambda f: _f(f['masa_kg']['H2_quemado'], 1)),
        fila('Agua retenida (kg)', lambda f: _f(f['masa_kg']['agua_retenida'], 0)),
        fila('H$_2$ regenerado (kg)', lambda f: _f(f['masa_kg']['H2_regenerado'], 1)),
        fila('Duración', lambda f: _dur(f['t_s'])),
        fila('Potencia de rechazo (MW)', lambda f: _f(f['P_W'] / 1e6, 0)),
        '\\midrule',
        fila('\\textit{Entra}: química, PCI (kWh)', lambda f: _f(f['kWh']['quimica_PCI'], 0)),
        fila('\\quad latente del agua condensada', lambda f: _f(f['kWh']['latente'], 0)),
        fila('\\quad sensible neto', lambda f: _f(f['kWh']['sensible_neto'], 0)),
        fila('\\quad compresor del escape', lambda f: _f(f['kWh']['compresor_escape'], 0)),
        fila('\\textit{Sale}: trabajo del motor', lambda f: rr(f['kWh']['trabajo_motor'])),
        fila('\\quad calor rechazado', lambda f: rr(f['kWh']['calor_rechazado'])),
        fila('\\textit{Usa}: ventiladores', lambda f: rr(f['kWh']['ventiladores'])),
        fila('\\quad hélices', lambda f: rr(f['kWh']['helices'])),
        fila('\\textbf{Neto a la planta}', lambda f: '\\textbf{' + rr(f['kWh']['neto_planta']) + '}'),
        '\\midrule',
        fila('\\textit{Regenera}: electrólisis', lambda f: _f(f['kWh']['electrolisis'], 0)),
        fila('\\quad compresión del H$_2$', lambda f: rr(f['kWh']['compresion_H2'])),
        fila('\\quad vuelve a los tanques (PCI)', lambda f: _f(quim(f), 0)),
        fila('\\quad pérdidas de la electrólisis', lambda f: _f(f['kWh']['electrolisis'] - quim(f), 0)),
        fila('\\textbf{Electricidad de tierra neta}', lambda f: '\\textbf{' + rr(extf(f)) + '}'),
    ])
    # ── tabla B: la suelta de 1 t contra el tiempo de recuperación ──
    _tabla('suelta', 'r r r r r r r',
           'Recuperación & Rechazo (MW) & Absorbe la quema (kg) & Hélices (MW) & Neto motor (MW) & Hélices (kWh) & '
           'Neto motor (kWh)',
           [f"{_dur(b['t_s'])} & {_f(b['P_W'] / 1e6, 2)} & {_f(b['m_e_kg'], 0)} & {_r(*[x / 1e6 for x in b['P_hel_W']], 2)} & "
            f"{_r(*[x / 1e6 for x in b['P_neto_W']], 2)} & {_r(*b['E_hel_kWh'], 0)} & {_r(*d['neto_kWh'], 0)}"
            for b in d['B_tabla']])
    # ── tabla C: tamaño de la suelta ──
    _tabla('tamano', 'r r r r r',
           'Suelta (t) & Con 10~MW: recuperación & hélices (kWh/t) & Con 10~min/t: hélices (MW) & hélices (kWh/t)',
           [f"{c['m_t']:g} & {_dur(c['t10_s'])} & {_r(*c['E10_kWh_t'], 0)} & {_r(*[x / 1e6 for x in c['PT_hel_W']], 2)} & "
            f"{_r(*c['ET_kWh_t'], 0)}" for c in d['C']])
    # ── tabla D: regeneración contra su duración ──
    _tabla('regeneracion', 'r r r', 'Duración (h) & Electrolizador (kW por t) & PEM (t por t)',
           [f"{x['t_h']:g} & {_f(x['P_el_kW'], 0)} & {_r(*x['masa_PEM_t'], 2)}" for x in d['D']])
    # ── tabla E: la entrega de referencia ──
    M_ = lambda v, dec=1: _r(v[0] / 1e3, v[1] / 1e3, dec)
    _tabla('entrega', 'l r r', ' & Por suelta & Entrega completa', [
        f"Duración & {E['t_paso_s'] / 60:.0f}~min & {E['n'] * E['t_paso_s'] / 3600:.0f}~h",
        f"Rechazo instalado & {_f(ps['P_rech_W'] / 1e6, 1)}~MW & --",
        f"Calor rechazado con el motor & {_r(*ps['calor_kWh'], 0)}~kWh & {M_(tt['calor_kWh'])}~MWh",
        f"H$_2$ quemado & {_f(ps['H2_kg'], 1)}~kg & {_f(tt['H2_kg'] / 1e3, 2)}~t",
        f"Agua retenida como lastre & {_f(ps['agua_kg'], 0)}~kg & {_f(tt['agua_kg'] / 1e3, 1)}~t",
        f"Hélices al soltar & {_r(*[x / 1e6 for x in ps['P_hel_W']], 2)}~MW & --",
        f"Hélices, energía & {_r(*ps['hel_kWh'], 1)}~kWh & {M_(tt['hel_kWh'], 2)}~MWh",
        f"Neto del motor, sin hélices & {_r(*ps['motor_kWh'], 0)}~kWh & {M_(tt['motor_kWh'])}~MWh",
        f"\\textbf{{Neto a la planta}} & \\textbf{{{_r(*ps['neto_kWh'], 0)}~kWh}} & \\textbf{{{M_(tt['neto_kWh'])}~MWh}}",
        '\\midrule',
        f"Electrólisis & {_f(ps['electrolisis_kWh'], 0)}~kWh & {_f(tt['electrolisis_kWh'] / 1e3, 1)}~MWh",
        f"Compresión del H$_2$ & {_r(*ps['compresion_kWh'], 0)}~kWh & {M_(tt['compresion_kWh'])}~MWh",
        f"\\textbf{{Electricidad de tierra neta}} & \\textbf{{{_r(*ps['externa_kWh'], 0)}~kWh}} & "
        f"\\textbf{{{M_(tt['externa_kWh'])}~MWh}}",
        f"Regeneración con el PEM del equipo ({_f(E['P_pem_kW'], 0)}~kW) & {p.t_regeneracion_h:.0f}~h & "
        f"{E['t_regen_pem_h'] / 24:.0f}~días",
    ])
    # ── tabla F: potencia por tonelada y la suelta más grande sin freno ──
    _tabla('por_tonelada', 'l r r r r r',
           'Nave & Carga útil (t) & Sin freno, 1~t (MW/t) & Cabe sin freno (kg) & En el cruce (MW/t) & (MW)',
           [f"{nv['nombre']} & {_f(nv['carga_t'], 1)} & {_f(nv['P_1t_MW_t'], 0)} & {_f(1000 * nv['cruce_t'], 0)} & "
            f"{_f(nv['cruce_MW_t'], 0)} & {_f(nv['cruce_t'] * nv['cruce_MW_t'], 1)}" for nv in F['naves'].values()])


def main():
    d = datos()
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    ch = chequeos(d)
    figuras(d)
    exportar_tex(d)
    txt = informe(d, ch)
    rutas.guardar_datos('exp22_energia_potencia', {k: v for k, v in d.items() if k != 'd21'})
    with open(rutas.tabla('exp22_energia_potencia').replace('.tex', '.txt'), 'w', encoding='utf-8') as fh:
        fh.write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
