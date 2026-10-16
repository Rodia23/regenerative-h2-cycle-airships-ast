# [EN] Mass and energy balance per tonne, by stage (Section 6).
# =============================================================================
# exp21_balance_energia.py — EXP-21: balance de masa y energía POR TONELADA liberada en los tres arreglos de la
# sección 5 (autor, 2026-10-12): NHD a 14 atm (sin freno, a la tasa de la cota ideal), HD a 14 atm y HD a 1 atm (con
# freno y la potencia de referencia de 10 MW). APARTE: no está en el paper ni en correr_todo.py.
#
# Por tonelada la energía no depende del ritmo: la misma masa de H2 se quema y el mismo calor se rechaza, en segundos
# (NHD) o en minutos (HD); el ritmo solo cambia la potencia. Cambian la retención (phi 0.99 a 14 atm, phi_max a 1 atm),
# el trabajo del compresor del escape (solo a 14 atm) y la energía de las hélices (solo con freno).
# Fuentes, las mismas del modelo: q_rech (EXP-05, el único balance de q), exp12_helio (rutas por tonelada), exp06 (trabajo
# isotérmico mínimo del compresor del escape), freno (hélices: energía con la fuerza en rampa a cero), exp16 (energía de
# compresión del H2 regenerado), parametros (PCI, L_v, eta del motor, electrólisis), balance_masa (f por caída de
# presión), fig_presupuesto_masa (potencia y duración sin freno).
# Chequeos: cierre del balance de energía, q por sus términos = q_rech, calor por tonelada = exp12, P x t = Q (con y
# sin freno), electrólisis a 1 atm = balance_masa, hélices = freno, masa de agua = r m_H2, neto = DoceNetaCiclo x t.
# Salidas: resultados/datos/exp21_balance_energia.json, resultados/tablas/exp21_balance_energia.txt
# =============================================================================
import sys

import parametros as p
import rutas

KWH = 3.6e6                                   # [J/kWh]


def terminos_q(phi):
    """Los términos de q_rech(0, phi) por kg de H2 [J/kg]: PCI, latente, sensible que entra y sensible que sale."""
    from exp05_condensacion_tiempo_real import AIRE_ESTEQ
    r_agua = 1.0 + p.ratio_O2_H2
    entra = (p.cp_H2 + AIRE_ESTEQ * p.cp_aire) * (p.T_aire_exterior - p.T_ref_comb)
    sale = ((AIRE_ESTEQ - p.ratio_O2_H2) * p.cp_N2
            + r_agua * (phi * p.cp_agua_liq + (1.0 - phi) * p.cp_vapor)) * (p.T_salida_cond - p.T_ref_comb)
    return {'PCI': p.Q_comb_H2, 'latente': phi * r_agua * p.L_v_agua, 'entra': entra, 'sale': sale}


def datos():
    from exp05_condensacion_tiempo_real import q_rech
    from balance_masa import f_caida_presion
    rt = rutas.leer_datos('exp12_helio')['rutas']
    r99 = next(r for r in rt if 'phi = 0.99' in r['ruta'])
    r1 = next(r for r in rt if '1 atm' in r['ruta'])
    w_comp = rutas.leer_datos('exp06_requisitos_minimos')['W_comp_MJkg'] * 1e6          # [J/kg H2], mínimo isotérmico
    e_comp = rutas.leer_datos('exp16_balance_neto')['e_comp_kWhkg']                      # [kWh/kg] regenerado a 350 bar
    etas = list(p.eta_motor_H2)
    fs = [min(f_caida_presion()), max(p.perdida_rad)]
    fr = rutas.leer_datos('freno')
    hel10 = next(x for x in fr['filas'] if x['corto'] == 'FW' and x['m_kg'] == 1000.0 and abs(x['P_W'] - 1e7) < 1.0)
    ap = {nv['corto']: nv for nv in rutas.leer_datos('fig_presupuesto_masa')['apilado']['naves']}
    sinFW = ap['FW']['sin']
    r_agua = 1.0 + p.ratio_O2_H2

    def arreglo(nombre, ruta, phi, comp, P_W, t_s, E_hel_J):
        m = ruta['kg_por_t']                                    # H2 quemado por tonelada [kg]
        tq = terminos_q(phi)
        q0 = q_rech(0.0, phi)
        W = w_comp if comp else 0.0
        trab = [e * p.Q_comb_H2 * m for e in etas]             # trabajo del motor [J]
        calor = [(q0 - e * p.Q_comb_H2 + W) * m for e in etas]  # calor rechazado con trabajo exportado [J]
        vent = [f * c for f, c in zip(fs, calor)]               # ventiladores (f mínima con eta mínima, como DoceNeta)
        vent_todo = sorted(f * c for f in fs for c in calor)
        neto = sorted(e * p.Q_comb_H2 * m - f * (q0 - e * p.Q_comb_H2 + W) * m - W * m - E for e in etas for f in fs
                      for E in E_hel_J)
        regen = phi * m                                         # H2 que la electrólisis devuelve [kg]
        repos = (1.0 - phi) * m                                 # H2 que se pierde y se repone desde los tanques [kg]
        return {
            'nombre': nombre, 'phi': phi, 'P_W': P_W, 't_s': t_s,
            'masa_kg': {'H2_quemado': m, 'agua_retenida': phi * r_agua * m, 'vapor_venteado': (1 - phi) * r_agua * m,
                        'H2_regenerado': regen, 'H2_reposicion': repos},
            'q_terminos_Jkg': tq, 'q0_Jkg': q0,
            'kWh': {'quimica_PCI': tq['PCI'] * m / KWH, 'latente': tq['latente'] * m / KWH,
                    'sensible_neto': (tq['entra'] - tq['sale']) * m / KWH, 'calor_eta0': q0 * m / KWH,
                    'compresor_escape': W * m / KWH, 'trabajo_motor': [x / KWH for x in trab],
                    'calor_rechazado': [x / KWH for x in calor], 'ventiladores': [vent_todo[0] / KWH, vent_todo[-1] / KWH],
                    'helices': [min(E_hel_J) / KWH, max(E_hel_J) / KWH], 'neto_planta': [neto[0] / KWH, neto[-1] / KWH],
                    'electrolisis': regen * p.E_H2_electrolisis, 'compresion_H2': [regen * e for e in e_comp],
                    'reposicion_tierra': repos * p.E_H2_electrolisis},
        }

    t14, t1 = hel10['t_rec_s'], r1['calor_GJ_por_t'] * 1e9 / 1e7
    E14 = hel10['E_hel_J']
    filas = [arreglo('NHD, 14 atm', r99, p.phi_objetivo, True, sinFW['P_MW'] * 1e6, sinFW['t_nat_s'], [0.0, 0.0]),
             arreglo('HD, 14 atm, 10 MW', r99, p.phi_objetivo, True, 1e7, t14, E14),
             # a 1 atm la misma fuerza durante más tiempo: la energía de las hélices escala con la recuperación
             arreglo('HD, 1 atm, 10 MW', r1, (r1['L'] - r99['L'] + (1 + p.ratio_O2_H2) * p.phi_objetivo) /
                     (1 + p.ratio_O2_H2), False, 1e7, t1, [e * t1 / t14 for e in E14])]
    return {'etas': etas, 'fs': fs, 'w_comp_Jkg': w_comp, 'e_comp_kWhkg': e_comp, 'filas': filas, 'r99': r99, 'r1': r1,
            'hel10': hel10}


def chequeos(d):
    """Cierres e identidades; devuelve [(nombre, obtenido, esperado, ok)]."""
    from exp05_condensacion_tiempo_real import q_rech
    out = []

    def ch(nombre, a, b, tol):
        out.append((nombre, a, b, abs(a - b) <= tol * max(abs(b), 1e-12)))

    for f in d['filas']:
        n, k, m = f['nombre'], f['kWh'], f['masa_kg']
        tq = f['q_terminos_Jkg']
        ch(f'{n}: q por sus términos = q_rech(0, phi)', tq['PCI'] + tq['latente'] + tq['entra'] - tq['sale'],
           q_rech(0.0, f['phi']), 1e-12)
        for i, e in enumerate(d['etas']):
            entra = k['quimica_PCI'] + k['latente'] + k['sensible_neto'] + k['compresor_escape']
            sale = k['calor_rechazado'][i] + k['trabajo_motor'][i]
            ch(f'{n}, eta {e}: entra (PCI + latente + sensible + compresor) = sale (calor + trabajo)', entra, sale, 1e-12)
        ch(f'{n}: agua retenida + venteada = r m_H2', m['agua_retenida'] + m['vapor_venteado'],
           (1 + p.ratio_O2_H2) * m['H2_quemado'], 1e-12)
        ch(f'{n}: H2 regenerado = agua retenida / r', m['H2_regenerado'], m['agua_retenida'] / (1 + p.ratio_O2_H2), 1e-12)
        ch(f'{n}: potencia x duración = calor por tonelada (eta 0)', f['P_W'] * f['t_s'] / KWH, k['calor_eta0'],
           2e-2 if n.startswith('NHD') else 1e-9)
    f14, f1 = d['filas'][1], d['filas'][2]
    ch('HD 14 atm: calor por tonelada = exp12 (ruta phi 0.99)', f14['kWh']['calor_eta0'] * KWH / 1e9,
       d['r99']['calor_GJ_por_t'], 1e-12)
    ch('HD 1 atm: calor por tonelada = exp12 (ruta 1 atm)', f1['kWh']['calor_eta0'] * KWH / 1e9, d['r1']['calor_GJ_por_t'],
       1e-9)
    ch('HD 1 atm: electrólisis por tonelada = balance_masa (1 atm)', f1['kWh']['electrolisis'],
       rutas.leer_datos('balance_masa')['E_el_por_t_kWh'], 1e-3)
    # el neto de la tabla de rutas (DoceNetaCiclo, en MW a 10 MW) no resta el compresor del escape ni las hélices
    import math
    from balance_masa import f_caida_presion
    from exp05_condensacion_tiempo_real import q_rech as _q
    _n = sorted(e * p.Q_comb_H2 - f * (_q(0.0, p.phi_objetivo) - e * p.Q_comb_H2) for e in d['etas'] for f in d['fs'])
    m14 = f14['masa_kg']['H2_quemado']
    ch('HD 14 atm: neto sin compresor ni hélices, mínimo = DoceNetaCiclo x t (kWh)', _n[0] * m14 / KWH,
       _n[0] * m14 / f14['t_s'] * f14['t_s'] / KWH, 1e-12)
    ch('HD 14 atm: neto de la tabla menos compresor y hélices = neto de esta tabla (mínimo)',
       _n[0] * m14 / KWH - f14['kWh']['compresor_escape'] - d['fs'][1] * f14['kWh']['compresor_escape']
       - f14['kWh']['helices'][1], f14['kWh']['neto_planta'][0], 1e-9)
    ch('HD 14 atm: hélices = freno (1 t, 10 MW, LCA60T)', f14['kWh']['helices'][1] * KWH, max(d['hel10']['E_hel_J']), 1e-12)
    return out


def tabla_txt(d):
    f = d['filas']
    rng = lambda v, fmt='{:.0f}': (fmt.format(v[0]) + '–' + fmt.format(v[1])) if isinstance(v, list) else fmt.format(v)
    filas = [('H2 quemado (kg)', lambda x: rng(x['masa_kg']['H2_quemado'], '{:.1f}')),
             ('agua retenida (kg)', lambda x: rng(x['masa_kg']['agua_retenida'])),
             ('vapor venteado (kg)', lambda x: rng(x['masa_kg']['vapor_venteado'])),
             ('H2 de reposición, de los tanques (kg)', lambda x: rng(x['masa_kg']['H2_reposicion'], '{:.2f}')),
             ('duración', lambda x: f"{x['t_s']:.0f} s" if x['t_s'] < 120 else f"{x['t_s'] / 60:.1f} min"),
             ('potencia de rechazo (MW)', lambda x: rng(x['P_W'] / 1e6)),
             ('ENTRA: energía química, PCI (kWh)', lambda x: rng(x['kWh']['quimica_PCI'])),
             ('ENTRA: latente del agua condensada (kWh)', lambda x: rng(x['kWh']['latente'])),
             ('ENTRA: sensible neto (kWh)', lambda x: rng(x['kWh']['sensible_neto'])),
             ('ENTRA: compresor del escape, mín. (kWh)', lambda x: rng(x['kWh']['compresor_escape'])),
             ('= TOTAL que entra (kWh)', lambda x: rng(x['kWh']['quimica_PCI'] + x['kWh']['latente']
                                                      + x['kWh']['sensible_neto'] + x['kWh']['compresor_escape'])),
             ('SALE: trabajo del motor (kWh)', lambda x: rng(x['kWh']['trabajo_motor'])),
             ('SALE: calor rechazado (kWh)', lambda x: rng(sorted(x['kWh']['calor_rechazado']))),
             ('USO: ventiladores (kWh)', lambda x: rng(x['kWh']['ventiladores'])),
             ('USO: hélices (kWh)', lambda x: rng(x['kWh']['helices'])),
             ('NETO a la planta (kWh)', lambda x: rng(x['kWh']['neto_planta'])),
             ('DESPUÉS: electrólisis (kWh)', lambda x: rng(x['kWh']['electrolisis'])),
             ('DESPUÉS: compresión del H2 a 350 bar (kWh)', lambda x: rng(x['kWh']['compresion_H2'])),
             ('DESPUÉS: reposición desde tierra (kWh)', lambda x: rng(x['kWh']['reposicion_tierra']))]
    L = ['EXP-21 — BALANCE POR TONELADA LIBERADA (modelo/exp21_balance_energia.py; aparte, no está en el paper)',
         f"eta del motor {d['etas']}, f de los ventiladores {[round(x, 3) for x in d['fs']]}", '',
         f"{'':44s}" + ''.join(f"{x['nombre']:>22s}" for x in f)]
    for nom, fn in filas:
        L.append(f'{nom:44s}' + ''.join(f'{fn(x):>22s}' for x in f))
    return '\n'.join(L)


def main():
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    d = datos()
    ch = chequeos(d)
    txt = tabla_txt(d) + '\n\nCHEQUEOS\n' + '\n'.join(
        f"  {'ok  ' if ok else 'FALLA'} {n}: {a:.6g} vs {b:.6g}" for n, a, b, ok in ch) + \
        f"\n{sum(ok for *_, ok in ch)}/{len(ch)} chequeos ok"
    rutas.guardar_datos('exp21_balance_energia', {k: v for k, v in d.items() if k not in ('r99', 'r1', 'hel10')})
    with open(rutas.tabla('exp21_balance_energia').replace('.tex', '.txt'), 'w', encoding='utf-8') as fh:
        fh.write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
