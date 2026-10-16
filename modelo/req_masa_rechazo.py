# -*- coding: utf-8 -*-
# [EN] Mass requirement of the heat-rejection unit and gap with existing classes (Section 5.1, Tables S7-S8).
# =============================================================================
# req_masa_rechazo.py — requisito de masa de la unidad de rechazo y brecha con lo que existe (paper 1).
# Ruta nueva del autor (2026-09-28): el paper no elige equipo para la caja de rechazo; deriva el requisito y lo
# compara con las clases que existen (sección 5.1 del paper; tablas S7 y S8).
#
# Presupuesto: la carga útil del caso base, menos las cajas fijas (regeneración, tanques, compresor), repartida
# entre las cajas que escalan con el flujo, POR kg/s de H2 quemado (el conjunto, no una caja sola):
#     B = (carga - fijas) / m*                                   [kg por kg/s]
# Arreglo con quemador (eta = 0): la batería (buffer) mueve los ventiladores, que gastan la fracción f del calor;
# con la clase de batería más liviana, lo que queda es para la unidad de rechazo:
#     kW/kg pedido (f) = (q m* / m*) / (B - m_buffer(f))          (curva; no hay reparto si m_buffer(f) >= B)
# Con motor no hay requisito: la clase de motor más liviana ya pesa más que B por kg/s.
# Otra fuente para los ventiladores (I-2, auditoría 3): una turbina de la clase CT7 que SOLO los mueve, con su escape
# fuera del condensador; pesa mucho menos que la batería y da el extremo bajo del requisito (quema H2 sin retener su
# agua). La tabla 7 da el requisito de la turbina a la batería y agrupa las clases por su temperatura de rechazo: sobre
# la banda de condensación del ciclo (las de avión, 80-160 C) o dentro de ella (I-3).
# Candidatos (cada uno con fuente): FlyZero (sistema y radiador, 2026 y 2035), el condensador enfriado por aire de
# una central (Zhang et al. 2020) y la recuperación de agua de los rígidos (Goodyear 1975, convertida con el calor
# por kg de agua de este ciclo).
# Salidas: tabla 6 del paper (art_masas: todas las cajas por kg/s de H2, el rechazo sin elegir; antes la generaba
# balance_masa.py en toneladas), tabla de candidatos (art_brecha_rechazo) y datos. La figura (req_masa_rechazo) no va
# al paper: el autor eligió la tabla (2026-09-29).
# =============================================================================

import sys
import numpy as np
import parametros as p
import rutas
import fuentes
from balance_masa import masa_buffer, f_caida_presion, FIJAS
from exp05_condensacion_tiempo_real import T_rocio, phi_max, AIRE_ESTEQ, q_rech
from exp02_quema_centrada import L_ef

V_CONSUMO_KN = 50.0                       # nota *** de Goodyear 1975, p. 19: consumo a 50 nudos
AGUA_POR_KG_METANO = 2.0 * 18.015 / 16.043   # cota: ningún hidrocarburo da más agua por kg que el metano (CH4 -> 2 H2O)


def main():
    bm = rutas.leer_datos('balance_masa')['ideal']
    e6 = rutas.leer_datos('exp06_requisitos_minimos')
    nave = p.NAVE_BASE
    mdot = bm['mdot_kgps']
    carga = nave['carga_max']
    fijo = bm['cabe_en_carga']['quemador_liviano']['fijo_kg']
    # (2026-10-14, objeción 1.8) la tonelada que se suelta va a bordo hasta la suelta: sale del presupuesto, como en el
    # desglose apilado (fig_presupuesto_masa) y en el flujo que cabe (balance_masa)
    m_suelta = bm['m_kg']
    B = (carga - fijo - m_suelta) / mdot                   # [kg por kg/s]
    q_q = bm['P_rej_eta0_W'] / mdot                        # calor por kg/s con quemador [W por kg/s]
    kw_bat = max(p.kWkg_bateria)                           # clase de batería más liviana (paquete)
    t_buf = bm['t_buf_s']

    def m_buf(f):
        return masa_buffer(f * q_q, t_buf, kw_bat)         # [kg por kg/s]

    def k_req(f):
        resto = B - m_buf(f)
        return q_q / 1e3 / resto if resto > 0 else float('inf')

    # sensibilidad (ítem 8, autor 2026-09-29): con la mejor celda encontrada para 30 s (Ionblox), llevada a paquete
    kw_mejor = p.kWkg_celda_30s_mejor * p.f_celda_paquete
    wh_mejor = p.Whkg_celda_30s_mejor * p.f_celda_paquete

    def k_req_mejor(f):
        m_b = max(f * q_q * t_buf / 3600.0 / wh_mejor, f * q_q / 1e3 / kw_mejor)
        return q_q / 1e3 / (B - m_b) if B > m_b else float('inf')

    # ventiladores movidos por una turbina chica de la clase CT7 que SOLO los mueve (I-2, auditoría 3): su escape no pasa
    # por el condensador, así que no diluye el vapor, pero quema H2 que no deja agua. Con la de mayor kW/kg de eje
    # (la más liviana) el requisito del rechazo es el más bajo de todos los arreglos
    kw_turb = sorted(w / 1e3 / m for w, m in p.TURBOEJE_TCDS)     # [kW/kg] de eje, de menor a mayor

    def m_turb(f, kw):
        return f * q_q / 1e3 / kw                                  # [kg por kg/s]

    def k_req_turb(f, kw):
        resto = B - m_turb(f, kw)
        return q_q / 1e3 / resto if resto > 0 else float('inf')

    # f físico (R-03, revisión 2026-10-01): una unidad en estacionario gasta solo la caída de presión de su
    # intercambiador, f = dp / (eta rho cp dT), con el aire de enfriamiento a la temperatura exterior
    rho_ext = p.P_atmos / (p.R_aire * (p.T_aire_exterior + 273.15))
    f_fis = f_caida_presion()
    aire_kgps = q_q / (p.cp_aire * p.dT_aire_enfr) * mdot           # aire de enfriamiento del caso base [kg/s]

    pendiente = max(q_q * t_buf / 3600.0 / p.Whkg_bateria, q_q / 1e3 / kw_bat)   # kg por kg/s por unidad de f
    f_max = B / pendiente
    conv_min = min(c['min_kg'] for c in bm['cajas'] if c['caja'] == 'conversion') / mdot

    # ── candidatos ──
    cand = []
    for i, anio in enumerate(('2026', '2035')):
        cand.append({'clave': f'FZ_sistema_{anio}', 'familia': 'FlyZero, system',
                     'nombre': f'Aircraft TMS, {anio}', 'T': '100', 'temp_ciclo': False,
                     'aire': 'ram', 'f': p.perdida_rad[i], 'kWkg': (p.kWkg_TMS_sistema[i],) * 2, 'cita': 'FlyZero2022'})
    for i, (anio, T) in enumerate((('2026', '80'), ('2035', '160'))):
        cand.append({'clave': f'FZ_radiador_{anio}', 'familia': 'FlyZero, radiator',
                     'nombre': f'Aircraft air radiator only, {anio}', 'T': T, 'temp_ciclo': False,
                     'aire': 'ram', 'f': p.perdida_rad[i], 'kWkg': (p.kWkg_radiador[i],) * 2, 'cita': 'FlyZero2022'})
    ms = p.ACC_masas_t
    M_acc = sum(ms.values()) * 1e3
    M_haces = (ms['haces_acero'] + ms['haces_aluminio']) * 1e3
    f_acc = p.ACC_vent_MWe / p.ACC_Q_MWth
    T_acc = f'{p.ACC_T_aire_C + p.ACC_ITD_C:.0f}'
    for clave, nombre, M in (('ACC_sistema', 'Power-plant ACC, whole system', M_acc),
                             ('ACC_haces', 'Power-plant ACC, finned-tube bundles', M_haces)):
        cand.append({'clave': clave, 'familia': 'Power-plant ACC', 'nombre': nombre, 'T': T_acc, 'aire': 'fans',
                     'temp_ciclo': p.ACC_T_aire_C + p.ACC_ITD_C <= T_rocio(p.P_atmos)[0],
                     'f': f_acc, 'kWkg': (p.ACC_Q_MWth * 1e3 / M,) * 2, 'cita': 'ZhangSpatari2020', 'masa_t': M / 1e3})
    # rígidos: masa del equipo por kg/s de agua (agua ~ combustible a 50 nudos; cota favorable: potencia máxima y el
    # agua del metano), convertida a kW/kg con el calor por kg de agua de ESTE ciclo
    gy = fuentes.FUENTES['GY75']['datos']
    agua, cons, vmax = gy['agua_lb'][0], gy['consumo_lbh'][0], gy['vmax_kn'][0]
    agua_kgps = mdot * (1.0 + p.ratio_O2_H2) * p.phi_objetivo
    q_agua = bm['P_rej_eta0_W'] / agua_kgps                  # [J por kg de agua]
    t_por_kgps, t_fav = {}, {}
    for n in agua:
        base = agua[n] * p.LB / (cons[n] * p.LB / 3600.0)                  # kg de equipo por kg/s de agua
        t_por_kgps[n] = base
        t_fav[n] = base / ((vmax[n] / V_CONSUMO_KN) ** 3 * AGUA_POR_KG_METANO)
    kw_hist = (q_agua / 1e3 / max(t_por_kgps.values()), q_agua / 1e3 / min(t_por_kgps.values()))
    kw_hist_fav = q_agua / 1e3 / min(t_fav.values())
    cand.append({'clave': 'rigidos', 'familia': 'Historical rigids', 'nombre': 'Rigid airships, water recovery',
                 'T': 'exhaust', 'aire': 'ram', 'f': None, 'kWkg': kw_hist, 'kWkg_favorable': kw_hist_fav,
                 'temp_ciclo': True,                  # condensaban el agua del escape: bajo su punto de rocío
                 'cita': 'Goodyear1975', 't_por_kgps': t_por_kgps, 't_por_kgps_favorable': t_fav})
    for c in cand:
        fpide = 0.0 if c['f'] is None else c['f']            # rígidos: f no reportado -> el caso más favorable
        c['pide'] = k_req(fpide)
        c['cabe'] = c['pide'] != float('inf')
        c['brecha'] = (c['pide'] / max(c.get('kWkg_favorable', 0.0), c['kWkg'][1]),
                       c['pide'] / c['kWkg'][0]) if c['cabe'] else None
        c['cumple'] = c['cabe'] and max(c.get('kWkg_favorable', 0.0), c['kWkg'][1]) >= c['pide']
        c['pide_bateria_mejor'] = k_req_mejor(fpide)
        c['cumple_bateria_mejor'] = max(c.get('kWkg_favorable', 0.0), c['kWkg'][1]) >= c['pide_bateria_mejor']
        # con la turbina más liviana para los ventiladores: el requisito más bajo (I-2)
        c['pide_turbina'] = k_req_turb(fpide, kw_turb[-1])
        c['cumple_turbina'] = max(c.get('kWkg_favorable', 0.0), c['kWkg'][1]) >= c['pide_turbina']
        tiene = (c['kWkg'][0], max(c.get('kWkg_favorable', 0.0), c['kWkg'][1]))
        c['brecha_rango'] = (c['pide_turbina'] / tiene[1], c['pide'] / tiene[0] if c['cabe'] else float('inf'))

    f35, f26 = p.perdida_rad[1], p.perdida_rad[0]
    turb = {'kWkg': (kw_turb[0], kw_turb[-1]),
            't_por_kgps_f2035': (m_turb(f35, kw_turb[-1]) / 1e3, m_turb(f35, kw_turb[0]) / 1e3),
            'kWkg_pide_f2035': (k_req_turb(f35, kw_turb[-1]), k_req_turb(f35, kw_turb[0])),
            'kWkg_pide_f2026': (k_req_turb(f26, kw_turb[-1]), k_req_turb(f26, kw_turb[0])),
            # H2 que quema la turbina por kg/s del quemador principal, con la eficiencia de conversión del paper
            'H2_extra_kgps_por_kgps': f35 * q_q / (p.rend_gen * p.Q_comb_H2),
            'cumplen': [c['clave'] for c in cand if c['cumple_turbina']]}
    ciclo = [c for c in cand if c['temp_ciclo']]
    brecha_ciclo = (min(c['brecha_rango'][0] for c in ciclo), max(c['brecha_rango'][1] for c in ciclo))

    # (2026-10-13) el requisito a cada presión (auditoría, objeción 5.2). Arriba va con el calor de phi objetivo (0.99),
    # que pide el escape a ~14 atm, pero sin su compresor (que solo ya no cabe: PilaCompEscSinFW). Sin compresor, a
    # 1 atm, se retiene phi_max: el calor por kg es q(0, phi_max) y la tasa sube como L(1) / L(phi_max) (ec. 7; s y
    # el tiempo del buffer no cambian, porque L m* tampoco). evaluar() rehace el requisito y la brecha de las clases a
    # la temperatura del ciclo con (tasa, calor, phi): con los de arriba reproduce lo de siempre (se comprueba aquí)
    def evaluar(mdot_x, q_x, phi_x):
        B_x = (carga - fijo - m_suelta) / mdot_x

        def kreq(f):
            resto = B_x - masa_buffer(f * q_x, t_buf, kw_bat)
            return q_x / 1e3 / resto if resto > 0 else float('inf')

        def kturb(f, kw):
            resto = B_x - f * q_x / 1e3 / kw
            return q_x / 1e3 / resto if resto > 0 else float('inf')

        def kmejor(f):
            m_b = max(f * q_x * t_buf / 3600.0 / wh_mejor, f * q_x / 1e3 / kw_mejor)
            return q_x / 1e3 / (B_x - m_b) if B_x > m_b else float('inf')
        q_agua_x = q_x / ((1.0 + p.ratio_O2_H2) * phi_x)          # los rígidos se convierten con el calor por kg de agua
        # (2026-10-14, suplementario A2) cada clase, no solo la brecha de las que condensan: las tablas S7 y S8 van a
        # 1 atm, el caso sin dispositivo del manuscrito (objeción 5.2)
        filas = []
        for c in cand:
            f = 0.0 if c['f'] is None else c['f']
            if c['clave'] == 'rigidos':
                kw_c = (q_agua_x / 1e3 / max(t_por_kgps.values()), q_agua_x / 1e3 / min(t_por_kgps.values()))
                fav_c = q_agua_x / 1e3 / min(t_fav.values())
            else:
                kw_c, fav_c = tuple(c['kWkg']), c.get('kWkg_favorable')
            tiene = (kw_c[0], max(fav_c or 0.0, kw_c[1]))
            pide, pide_t = kreq(f), kturb(f, kw_turb[-1])
            cabe = pide != float('inf')
            filas.append({'clave': c['clave'], 'kWkg': kw_c, 'kWkg_favorable': fav_c, 'pide': pide,
                          'pide_turbina': pide_t, 'cabe': cabe,
                          'brecha_rango': (pide_t / tiene[1], pide / tiene[0] if cabe else float('inf'))})
        br = [fl['brecha_rango'] for fl, c in zip(filas, cand) if c['temp_ciclo']]
        pend = max(q_x * t_buf / 3600.0 / p.Whkg_bateria, q_x / 1e3 / kw_bat)   # kg por kg/s por unidad de f
        return {'mdot_kgps': mdot_x, 'phi': phi_x, 'q_MW_por_kgps': q_x / 1e6, 'P_MW': mdot_x * q_x / 1e6,
                'presupuesto_t_por_kgps': B_x / 1e3, 'kWkg_pide_ffis': tuple(kreq(f) for f in f_fis),
                'brecha_temp_ciclo': (min(b[0] for b in br), max(b[1] for b in br)),
                'kWkg_pide_f0': kreq(0.0), 'kWkg_pide_f2035': kreq(p.perdida_rad[1]), 'f_max': B_x / pend,
                'kWkg_pide_bat_mejor_f2035': kmejor(p.perdida_rad[1]), 'q_agua_MJkg': q_agua_x / 1e6, 'filas': filas}
    ev = evaluar(mdot, q_q, p.phi_objetivo)
    _dif = [abs(a - b) for a, b in zip(ev['kWkg_pide_ffis'], [k_req(f) for f in f_fis])]
    _dif += [abs(a - b) for a, b in zip(ev['brecha_temp_ciclo'], brecha_ciclo)]
    _dif += [abs(ev['kWkg_pide_f0'] - k_req(0.0)), abs(ev['kWkg_pide_f2035'] - k_req(p.perdida_rad[1])),
             abs(ev['f_max'] - f_max), abs(ev['kWkg_pide_bat_mejor_f2035'] - k_req_mejor(p.perdida_rad[1])),
             abs(ev['q_agua_MJkg'] - q_agua / 1e6)]
    for fl, c in zip(ev['filas'], cand):
        _dif += [abs(fl['pide'] - c['pide']) if c['cabe'] else 0.0, abs(fl['pide_turbina'] - c['pide_turbina']),
                 abs(fl['brecha_rango'][0] - c['brecha_rango'][0]),
                 abs(fl['brecha_rango'][1] - c['brecha_rango'][1]) if c['cabe'] else 0.0]
    if max(_dif) > 1e-9:
        sys.exit('req_masa_rechazo: evaluar() no reproduce el requisito con el calor de phi objetivo')
    phi1 = phi_max(p.T_salida_cond)
    una_atm = evaluar(mdot * L_ef(1.0) / L_ef(phi1), q_rech(0.0, phi1), phi1)

    # ── requisitos y cajas por kg/s (tabla 6 del paper: todas las cajas; el rechazo, sin elegir) ──
    por_kgps = lambda cajas, nom: tuple(x / mdot for x in next((c['min_kg'], c['max_kg']) for c in cajas
                                                                if c['caja'] == nom))
    buf_q = por_kgps(bm['quemador']['cajas'], 'buffer')            # batería que mueve los ventiladores del quemador
    req = {'mdot_kgps': mdot, 'escape_kgps': mdot * (1.0 + AIRE_ESTEQ), 'agua_kgps': agua_kgps,
           'phi': p.phi_objetivo, 'P_quemador_MW': bm['P_rej_eta0_W'] / 1e6, 'P_motor_MW': bm['P_rej_W'] / 1e6,
           'eta_motor': p.rend_gen, 'T_aire_C': p.T_aire_exterior, 'T_rocio_1atm_C': T_rocio(p.P_atmos)[0],
           'T_salida_C': p.T_salida_cond, 'phi_1atm': phi_max(p.T_salida_cond), 'P_99_atm': e6['P_req_atm'],
           'T_rocio_99_C': T_rocio(e6['P_req_atm'] * p.P_atmos)[0], 'tau_s': p.tau_condensador,
           't_arranque_s': p.n_tau_arranque * p.tau_condensador, 'carga_t': carga / 1e3, 'fijas_t': fijo / 1e3, 'suelta_t': m_suelta / 1e3,
           'presupuesto_t_por_kgps': B / 1e3, 'q_MW_por_kgps': q_q / 1e6, 'kWkg_bateria': kw_bat,
           'kWkg_pide_f0': k_req(0.0), 'kWkg_pide_fACC': k_req(f_acc), 'kWkg_pide_f2035': k_req(p.perdida_rad[1]),
           'f_2035': p.perdida_rad[1], 'f_max': f_max, 'conversion_min_t_por_kgps': conv_min / 1e3,
           'q_agua_MJkg': q_agua / 1e6, 'v_consumo_kn': V_CONSUMO_KN, 'turbina_ventiladores': turb,
           'brecha_temp_ciclo': brecha_ciclo, 'una_atm': una_atm,
           'f_fisico': f_fis, 'kWkg_pide_ffis': tuple(k_req(f) for f in f_fis),
           'buffer_ffis_t_por_kgps': tuple(m_buf(f) / 1e3 for f in f_fis),
           'rho_aire_ext': rho_ext, 'aire_enfr_kgps': aire_kgps, 'aire_enfr_m3s': aire_kgps / rho_ext,
           'area_entrada_m2': aire_kgps / rho_ext / p.v_aire_entrada,
           'seccion_casco_m2': 3.141592653589793 * nave['D'] ** 2 / 4.0,
           'sens_bateria': {'kWkg_celda': p.kWkg_celda_30s_mejor, 'kWkg_paquete': kw_mejor,
                            'kWkg_pide_f2035': k_req_mejor(p.perdida_rad[1]),
                            'kWkg_pide_f2026': k_req_mejor(p.perdida_rad[0]),
                            'cumplen': [c['clave'] for c in cand if c['cumple_bateria_mejor']]},
           'conversion_t_por_kgps': tuple(x / 1e3 for x in por_kgps(bm['cajas'], 'conversion')),
           'buffer_quemador_t_por_kgps': tuple(x / 1e3 for x in buf_q),
           'queda_rechazo_t_por_kgps': (B - buf_q[0]) / 1e3,
           'turbina_t_por_kgps': tuple(x / mdot / 1e3 for x in bm['turbina_kg']),
           'fijas_kg': {c['caja']: (c['min_kg'], c['max_kg']) for c in bm['cajas']
                        if c['caja'] in FIJAS}}
    f_malla = np.linspace(0.0, 0.999 * f_max, 400)
    curva = {'f': list(f_malla), 'kWkg': [k_req(x) for x in f_malla], 'buffer_t_por_kgps': [m_buf(x) / 1e3 for x in f_malla]}
    # (2026-10-14, suplementario A2) las tablas S7 y S8 a 1 atm, el caso sin dispositivo del manuscrito: tasa, calor y
    # phi de una_atm; el buffer por kg/s escala con el calor (dimensionado por potencia o por energía, lineal en P)
    u_tab = {**req, 'mdot_kgps': una_atm['mdot_kgps'], 'presupuesto_t_por_kgps': una_atm['presupuesto_t_por_kgps'],
             'buffer_quemador_t_por_kgps': tuple(x * una_atm['q_MW_por_kgps'] / (q_q / 1e6)
                                                 for x in req['buffer_quemador_t_por_kgps'])}
    u_tab['queda_rechazo_t_por_kgps'] = u_tab['presupuesto_t_por_kgps'] - u_tab['buffer_quemador_t_por_kgps'][0]
    cand_1 = [{**c, **{k: v for k, v in fl.items() if not (k == 'kWkg_favorable' and v is None)}}
              for c, fl in zip(cand, una_atm['filas'])]
    una_atm['tabla_s7'] = {k: u_tab[k] for k in ('mdot_kgps', 'presupuesto_t_por_kgps', 'buffer_quemador_t_por_kgps',
                                                 'queda_rechazo_t_por_kgps')}
    rutas.guardar_datos('req_masa_rechazo', {'requisitos': req, 'candidatos': cand, 'curva': curva,
                                             'acc_masa_t': M_acc / 1e3, 'acc_haces_t': M_haces / 1e3, 'f_acc': f_acc})
    tabla_masas(u_tab)
    tabla_brecha(cand_1, u_tab)
    graficar(curva, cand, req)

    print(f"presupuesto {B/1e3:.2f} t por kg/s (carga {carga/1e3:.0f} t - fijas {fijo/1e3:.2f} t, m* {mdot:.3f} kg/s)")
    print(f"calor {q_q/1e6:.1f} MW por kg/s; batería {kw_bat:.3f} kW/kg; f_max {100*f_max:.1f} %; "
          f"motor más liviano {conv_min/1e3:.1f} t por kg/s")
    for f in (0.0, f_acc, 0.05, 0.10, p.perdida_rad[1], p.perdida_rad[0]):
        print(f"  f {100*f:5.1f} %: buffer {m_buf(f)/1e3:6.2f} t por kg/s -> pide {k_req(f):8.2f} kW/kg")
    for c in cand:
        kw = f"{c['kWkg'][0]:.2f}" + (f"-{c['kWkg'][1]:.2f}" if c['kWkg'][1] != c['kWkg'][0] else '')
        fc = '--' if c['f'] is None else f"{100 * c['f']:.1f}"
        print(f"  {c['nombre']:<52} f {fc:>5} % tiene {kw:>10} pide {c['pide']:7.2f} cumple {c['cumple']}")


def _r(a, b, dec):
    """Rango 'a--b' con el redondeo de la tabla; un solo número si coinciden."""
    ta, tb = f'{a:.{dec}f}', f'{b:.{dec}f}'
    return ta if ta == tb else f'{ta}--{tb}'


def tabla_masas(r):
    """Tabla 6 del paper: todas las cajas, las que escalan con el flujo en t por kg/s de H2 y las fijas en t; la unidad
    de rechazo no se elige (recibe lo que deja el presupuesto; candidatos en la tabla de brecha)."""
    cv, bf, fj = r['conversion_t_por_kgps'], r['buffer_quemador_t_por_kgps'], r['fijas_kg']
    # (autor, 2026-10-08) nombres cortos en "Representative class"; antes: "Li-ion pack, sized by power, driving the fans
    # of a burner", "Not chosen: what the budget leaves (Table 7)", "PEM, one tonne recovered per day", "Type IV, all the
    # hydrogen of a tonne" y "Electric, regenerated hydrogen to the tanks" (lo que decían va en el pie de la tabla)
    L = [r'\begin{tabular}{@{} l l r l @{}}', r'\toprule',
         r'\textbf{Black box} & \textbf{Representative class} & \textbf{Mass} & \textbf{Source} \\', r'\midrule',
         r'\multicolumn{4}{@{}l}{\textit{Scale with the burn rate (t per kg/s of hydrogen)}} \\',
         r'Energy conversion unit & H\textsubscript{2} reciprocating engine, $\lambda = 1$ & '
         f"{_r(cv[0], cv[1], 0)} & " + r'\cite{Rotax2023, DEUTZ2024} \\',
         r'Energy buffer & Li-ion pack (power-sized) & '
         f"{_r(bf[0], bf[1], 0)} & " + r'\cite{Pattanayak2025} \\',
         r'Heat-rejection unit & Not chosen (budget remainder) & '
         f"$\\le$ {r['queda_rechazo_t_por_kgps']:.1f} & -- " + r'\\',
         r'\addlinespace',
         r'\multicolumn{4}{@{}l}{\textit{Fixed, per tonne released (t)}} \\',
         r'Regeneration unit & Aviation PEM stack & '
         f"{_r(fj['regeneracion'][0] / 1e3, fj['regeneracion'][1] / 1e3, 2)} & " + r'\cite{Zhang2025PEMFC} \\',
         r'H\textsubscript{2} tanks & Type IV vessels & '
         f"{_r(fj['tanques_por_t'][0] / 1e3, fj['tanques_por_t'][1] / 1e3, 2)} & " + r'\cite{DOE_H2storage} \\',
         r'H\textsubscript{2} compressor & Electric industrial unit & '
         f"{_r(fj['compresor'][0] / 1e3, fj['compresor'][1] / 1e3, 2)} & " + r'\cite{HydroPac2011} \\',
         # D (2026-10-13): el agua retenida de la tonelada, en un tanque de tela colapsable
         r'Water tank & Collapsible fabric tank & '
         f"{_r(fj['tanque_agua'][0] / 1e3, fj['tanque_agua'][1] / 1e3, 3)} & " + r'\cite{MILT53029C} \\',
         r'\midrule',
         r'\multicolumn{2}{@{}l}{Budget for the boxes that scale: '
         # ANTES (revisión 2026-09-30, P-04: las cajas fijas más livianas suman 1.45 t; con 1.5 el lector no lo reconoce):
         # f"({r['carga_t']:.0f}~t $-$ {r['fijas_t']:.1f}~t) / {r['mdot_kgps']:.2f}~kg/s" + r'} & '
         f"({r['carga_t']:.0f}~t $-$ {r['suelta_t']:.0f}~t $-$ {r['fijas_t']:.2f}~t) / {r['mdot_kgps']:.2f}~kg/s"
         + r'} & '
         f"{r['presupuesto_t_por_kgps']:.1f} & " + r'\\',
         r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('art_masas'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/req_masa_rechazo.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


def _n(x):
    # ANTES (revisión 2026-09-30, P-03: el ACC de 0.135 kW/kg salía 0.14 y la brecha 57-60 no se podía rehacer):
    # return f'{x:.2f}' if x < 1 else (f'{x:.1f}' if x < 10 else f'{x:.0f}')
    return f'{x:.3g}' if x < 1 else (f'{x:.1f}' if x < 10 else f'{x:.0f}')     # bajo 1: tres cifras significativas


def _r2(a, b):
    """Rango con el redondeo de _n; un solo número si coinciden (rígidos, f = 0: sin ventiladores que mover)."""
    return _n(a) if _n(a) == _n(b) else f'{_n(a)}--{_n(b)}'


def _g(x):
    """Brecha (necesita / tiene): una cifra significativa bajo 1 (entonces cabe), entera arriba. (2026-10-14) Cerca de
    1 (0.95-1.05), dos decimales, para que la tabla distinga la clase que alcanza (0.96) de la que no (1.02)."""
    if 0.95 <= x < 1.05:
        return f'{x:.2f}'
    return f'{x:.1f}' if x < 10 else f'{x:.0f}'


DOS_LINEAS = {'Power-plant ACC, finned-tube bundles': r'\makecell[l]{Power-plant ACC,\\finned-tube bundles}',
              'Rigid airships, water recovery': r'\makecell[l]{Rigid airships,\\water recovery}'}


def tabla_brecha(cand, r):
    """Tabla 7: cada clase con lo que rinde y lo que el ciclo le pide a su f, de la turbina más liviana para los
    ventiladores a la batería (I-2, auditoría 3); agrupadas por dónde rechazan el calor: sobre la banda de
    condensación del ciclo o a su temperatura (I-3: el radiador de 2035 rechaza a 160 C)."""
    L = [r'\begin{tabular}{@{} l r l r r r r l @{}}', r'\toprule',
         r'\textbf{Candidate} & \textbf{$T$ (\textdegree C)} & \textbf{Air} & \textbf{$f$ (\%)} & '
         r'\textbf{Has} & \textbf{Needs} & \textbf{Gap} & \textbf{Source} \\',
         r' & & & & \multicolumn{2}{c}{(kW/kg)} & ($\times$) & \\', r'\midrule']
    grupos = ((False, r'\textit{Rejecting above the condensation band}'),
              (True, r'\textit{Condensing at the temperature of the cycle}'))
    for i, (en_ciclo, titulo) in enumerate(grupos):
        if i:
            L.append(r'\addlinespace')
        L.append(r'\multicolumn{8}{@{}l}{' + titulo + r'} \\')
        for c in (x for x in cand if x['temp_ciclo'] == en_ciclo):
            f = (r'--\textsuperscript{a}' if c['f'] is None else f"{100*c['f']:.0f}" if c['f'] >= 0.1
                 else f"{100*c['f']:.1f}")
            has = _n(c['kWkg'][0]) + ('' if c['kWkg'][1] == c['kWkg'][0] else f"--{_n(c['kWkg'][1])}")
            if 'kWkg_favorable' in c:
                # (autor, 2026-10-08) en dos líneas: el rango y debajo la cota favorable
                has = r'\makecell[r]{' + has + r'\\' + f"($\\le$ {_n(c['kWkg_favorable'])})" + '}'
            lo, hi = c['brecha_rango']
            if c['cabe']:
                needs = _r2(c['pide_turbina'], c['pide'])
                gap = _g(lo) if _g(lo) == _g(hi) else f'{_g(lo)}--{_g(hi)}'
            else:
                needs = f"{_n(c['pide_turbina'])}--none\\textsuperscript{{b}}"
                gap = f'{_g(lo)}--' + r'$\infty$'
            # (autor, 2026-10-08) los dos nombres largos de la última parte, en dos líneas
            nombre = DOS_LINEAS.get(c['nombre'], c['nombre'])
            L.append(f"{nombre} & {c['T']} & {c['aire']} & {f} & {has} & {needs} & {gap} & "
                     f"\\cite{{{c['cita']}}}" + r' \\')
    # las notas (a: f de los rígidos; b: sin reparto con batería) van en el pie de la tabla del manuscrito, con macros
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('art_brecha_rechazo'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/req_masa_rechazo.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


def graficar(curva, cand, r):
    import matplotlib.pyplot as plt
    import plot_style  # noqa
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    f = np.array(curva['f']) * 100
    k = np.array(curva['kWkg'])
    y_max = 300.0
    ax.fill_between(f, k, y_max, color='#009E73', alpha=0.12, lw=0, label='Fits in the payload')
    ax.plot(f, k, color='#009E73', lw=2, label='Required (burner, lightest battery)')
    ax.axvline(100 * r['f_max'], color='gray', ls='--', lw=1)
    ax.axvspan(100 * r['f_max'], 24, color='gray', alpha=0.08, lw=0)
    ax.text(100 * r['f_max'] + 0.3, 0.13, 'the battery alone\nfills the payload', fontsize=8, color='gray', va='bottom')
    estilo = {'FlyZero, system': ('s', '#0072B2'), 'FlyZero, radiator': ('D', '#56B4E9'),
              'Power-plant ACC': ('o', '#D55E00'), 'Historical rigids': ('^', '#E69F00')}
    vistos = set()
    for c in cand:
        m, col = estilo[c['familia']]
        x = 0.0 if c['f'] is None else 100 * c['f']
        lab = c['familia'] if c['familia'] not in vistos else None
        vistos.add(c['familia'])
        if c['familia'] == 'Historical rigids':
            ax.plot([x, x], [c['kWkg'][0], c['kWkg_favorable']], color=col, lw=3, alpha=0.6)
            ax.plot(x, c['kWkg'][0], m, color=col, ms=8, label=lab, clip_on=False)
            ax.plot(x, c['kWkg_favorable'], m, color=col, ms=8, mfc='white', clip_on=False)
            ax.annotate('rigid airships\n(f not reported)', (x, c['kWkg_favorable']), xytext=(8, 4),
                        textcoords='offset points', fontsize=8, color=col)
            continue
        ax.plot(x, c['kWkg'][0], m, color=col, ms=8, label=lab, clip_on=False)
        texto = {'FZ_sistema_2026': 'system 2026', 'FZ_sistema_2035': 'system 2035',
                 'FZ_radiador_2026': 'radiator 2026', 'FZ_radiador_2035': 'radiator 2035',
                 'ACC_sistema': 'whole system', 'ACC_haces': 'bundles only'}[c['clave']]
        izq = c['clave'].endswith('2035')                   # a la izquierda: no cruzar la línea de f_max
        ax.annotate(texto, (x, c['kWkg'][0]), xytext=(-8 if izq else 8, -3), textcoords='offset points', fontsize=8,
                    color=col, ha='right' if izq else 'left')
    ax.set_yscale('log')
    ax.set_ylim(0.1, y_max)
    ax.set_xlim(-1.0, 24)                                  # margen a la izquierda: los rígidos (f = 0) no quedan sobre el eje
    ax.set_xlabel('Fan power as a fraction of the heat rejected, $f$ (%)')
    ax.set_ylabel('Heat rejected per unit mass (kW/kg)')
    ax.legend(fontsize=8, loc='lower center', bbox_to_anchor=(0.45, 0.02))
    ax.grid(True, which='both', alpha=0.2)
    fig.tight_layout()
    fig.savefig(rutas.figura('req_masa_rechazo'), dpi=200, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
