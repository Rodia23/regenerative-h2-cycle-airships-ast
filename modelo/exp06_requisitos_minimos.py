# -*- coding: utf-8 -*-
# [EN] Real-time requirements when the water is condensed (Section 5.1).
# =============================================================================
# exp06_requisitos_minimos.py — EXP-06 (análisis, SIN simular): requisitos
# mínimos de la maniobra base con la condición de CONDENSAR TODO.
#
# Base (acordada con el autor): el minimax del artículo, que en su forma cerrada
# solo quema (sin reinyección). Restricción: el agua se condensa al ritmo en que
# se produce (φ = phi_objetivo, en tiempo real). Se reportan, por plataforma,
# carga y umbral h:
#   flujos     ṁ* (minimax EXACTO, con la corrección ε del Ap. E), H2 por entrega B = m/L, potencia térmica ṁ*·PCI,
#              aire estequiométrico, agua producida
#   condensar  P_cond = ṁ*·q (q = PCI + 9 φ L_v, sin exportar trabajo),
#              aire de enfriamiento P/(cp ΔT_aire), UA = P/ΔT_m,
#              potencia específica P/m (kW por kg de carga)
#   manera     presión del escape necesaria para φ a T_salida (con aire exterior
#              no alcanza a 1 atm) y su trabajo mínimo de compresión (interno:
#              no cambia q si no se exporta trabajo)
# D_lim = D_lim_oper (artículo); columna aparte con D_lim = h/3 (mejor verificado en
# EXP-01): todo lo térmico escala con ṁ*.
# Se evalúan fórmulas (instantáneo): no es una simulación.
# =============================================================================

import io
import sys
import csv
import math
import parametros as p
import rutas
import minimax as mm
from exp02_quema_centrada import L_ef, cargas_eval
from exp05_condensacion_tiempo_real import p_sat, phi_max, N_INERTE, AIRE_ESTEQ, q_rech, w_comp_escape

H_UMBRALES = (p.h_max_precision, p.h_max_seguridad)
P_INSTALADA_MW = (1.0, 10.0, 100.0)   # potencias de condensación instaladas de referencia


def paso_maximo(nave, h, P_cap, q):
    """Mayor masa por paso [kg] cuyo minimax (±h, D_lim oper.) se condensa en tiempo real con P_cap [W].
    nave: nombre o dict de parametros.PLATAFORMAS (el paso no puede pasar de lo que la nave levanta)."""
    d = mm.nave(nave)
    lo, hi = 1.0, min(60000.0, d['sustentacion'] - d['m_vacio'])
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if mm.resolver(d, mid, h, p.D_lim_oper)['dm_star'] * q <= P_cap else (lo, mid)
    return lo


def presion_para_phi(phi, T_C):
    """Presión total [Pa] a la que, enfriando a T_C, condensa la fracción phi."""
    return p_sat(T_C) * (1.0 + N_INERTE / (1.0 - phi))


def p_sat_buck(T_C):
    """Buck (1996) sobre agua líquida; solo para ubicar temperaturas bajo 0 °C."""
    return 611.21 * math.exp((18.678 - T_C / 234.5) * (T_C / (257.14 + T_C)))


def T_para_phi(phi, P=p.P_atmos):
    """Temperatura [°C] a la que hay que enfriar a presión P para condensar phi."""
    ps = (1.0 - phi) * P / (N_INERTE + 1.0 - phi)
    lo, hi = -40.0, 100.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if p_sat_buck(mid) < ps else (lo, mid)
    return 0.5 * (lo + hi)


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    phi = p.phi_objetivo
    L = L_ef(phi)
    q = q_rech(0.0, phi)
    P_req = presion_para_phi(phi, p.T_salida_cond)
    T_1atm = T_para_phi(phi)
    W_comp = w_comp_escape(P_req, p.T_salida_cond)     # mínimo: solo la fase gaseosa (M-1, auditoría 3)
    f_ma = 1.0 / math.sqrt(1.0 + p.NAVE_BASE['k'])      # factor de la masa añadida en el caso base

    filas = []
    for n, d in p.PLATAFORMAS.items():
        for m_t in cargas_eval(n):
            m = m_t * 1000.0
            for h in H_UMBRALES:
                r = mm.resolver(n, m, h, p.D_lim_oper)
                dm = r['dm_star']                                   # con masa añadida (inercia.py)
                dm_sin = mm.resolver(n, m, h, p.D_lim_oper, k_ma=0.0)['dm_star']
                P = dm * q
                dm3 = mm.resolver(n, m, h, h / 3.0)['dm_star']     # mejor D_lim verificado (EXP-01)
                filas.append({
                    'plataforma': n, 'corto': d['corto'], 'carga_t': m_t, 'h_m': h, 'M_total_kg': r['M_total'],
                    'dm_star_kgps': dm, 'dm_star_aprox_kgps': r['dm_star_aprox'], 'dm_star_sin_ma_kgps': dm_sin,
                    'H2_kg': m / L, 'P_termica_MW': dm * p.Q_comb_H2 / 1e6,
                    'aire_comb_kgps': dm * AIRE_ESTEQ, 'agua_kgps': dm * (1.0 + p.ratio_O2_H2) * phi,
                    'P_cond_MW': P / 1e6, 'P_cond_sin_ma_MW': dm_sin * q / 1e6,
                    'aire_enfr_tps': P / (p.cp_aire * p.dT_aire_enfr) / 1e3,
                    'ratio_aire_enfr_comb': (P / (p.cp_aire * p.dT_aire_enfr)) / (dm * AIRE_ESTEQ),
                    'UA_MWK': P / p.dT_media_cond / 1e6,
                    'P_especifica_kWkg': P / m / 1e3,
                    'Q_GJ': m / L * q / 1e9,
                    'dm_star_Dh3_kgps': dm3, 'P_cond_Dh3_MW': dm3 * q / 1e6,
                })

    # Liberación fraccionada con el minimax por paso: cada paso vuelve a y0 (±h por paso), así que el
    # flujo por paso escala ~ (m/n)^(3/2). Con una potencia instalada P_cap: paso máximo, número de
    # pasos y duración = max(n · intervalo_multidrop, Q/P_cap) (el calor total no cambia).
    pasos = []
    for n, d in p.PLATAFORMAS.items():
        m_tot = d['carga_max']       # la carga útil COMPLETA (I-1, auditoría 3: antes cargas_eval, con tope de 30 t)
        for h in H_UMBRALES:
            for P_mw in P_INSTALADA_MW:
                ms = min(paso_maximo(n, h, P_mw * 1e6, q), m_tot)
                n_p = math.ceil(m_tot / ms - 1e-9)
                Q = m_tot / L * q
                pasos.append({'plataforma': n, 'corto': d['corto'], 'carga_t': m_tot / 1e3, 'h_m': h,
                              'P_instalada_MW': P_mw, 'paso_max_kg': ms, 'n_pasos': n_p,
                              't_pasos_min': n_p * p.intervalo_multidrop / 60.0, 't_calor_min': Q / (P_mw * 1e6) / 60.0,
                              't_total_min': max(n_p * p.intervalo_multidrop, Q / (P_mw * 1e6)) / 60.0})
    with open(rutas.dato('exp06_pasos', 'csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(pasos[0].keys()))
        w.writeheader()
        for r in pasos:
            w.writerow({kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in r.items()})

    with open(rutas.dato('exp06_requisitos_minimos', 'csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        for r in filas:
            w.writerow({kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in r.items()})

    def buscar(nombre, m_t, h):
        return next(f for f in filas if f['plataforma'] == nombre and f['carga_t'] == m_t and f['h_m'] == h)
    pf1, pf5 = buscar(p.nombre_dirigible, 1.0, 1.0), buscar(p.nombre_dirigible, 1.0, 5.0)
    # R-04 (revisión 2026-10-01): la cota ideal junta phi = 1 en la tasa (L_eff) con el calor de condensar phi
    # objetivo, que a 1 atm pide 14 atm o -2 °C. Caso coherente a 1 atm: el agua que condensa ahí (phi_max a
    # T_salida) en la tasa y en el calor.
    phi_1 = phi_max(p.T_salida_cond)
    r_1 = mm.resolver(p.nombre_dirigible, 1000.0, p.h_max_precision, p.D_lim_oper, rho_gas=mm.rho_H2,
                      L_gas=L_ef(phi_1))
    base_1atm = {'phi': phi_1, 'L': L_ef(phi_1), 'dm_star_kgps': r_1['dm_star'], 'q_MJkg': q_rech(0.0, phi_1) / 1e6,
                 'P_cond_MW': r_1['dm_star'] * q_rech(0.0, phi_1) / 1e6}
    fw1, fw5 = buscar('Flying Whales LCA60T', 20.0, 1.0), buscar('Flying Whales LCA60T', 20.0, 5.0)
    rutas.guardar_datos('exp06_requisitos_minimos', {
        'phi_objetivo': phi, 'L_phi': L, 'q_MJkg': q / 1e6, 'P_req_atm': P_req / p.P_atmos,
        'T_1atm_C': T_1atm, 'W_comp_MJkg': W_comp / 1e6, 'W_comp_frac_PCI': W_comp / p.Q_comb_H2,
        'factor_masa_anadida': f_ma, 'base_1t_h1': pf1, 'base_1t_h5': pf5, 'FW_20t_h1': fw1, 'FW_20t_h5': fw5,
        'base_1t_1atm': base_1atm,
        'pasos': pasos,
        # amortiguador térmico: agua a bordo que absorba el calor de 1 t calentándose 50 K
        'buffer_agua_t_por_t': 1000.0 / L * q / (p.cp_agua_liq * p.dT_buffer_agua) / 1e3,
    })
    tabla(filas)
    tabla_pasos(pasos)

    print(f'Condensar φ = {phi}: a {p.T_salida_cond:g} °C hace falta {P_req/p.P_atmos:.1f} atm '
          f'(a 1 atm habría que enfriar a {T_1atm:.1f} °C); trabajo mínimo de compresión {W_comp/1e6:.1f} MJ/kg '
          f'({100*W_comp/p.Q_comb_H2:.0f} % del PCI). q = {q/1e6:.1f} MJ/kg, L = {L:.2f}')
    print(f'Caso coherente a 1 atm (φ = φ_max = {base_1atm["phi"]:.3f}): ṁ* {base_1atm["dm_star_kgps"]:.2f} kg/s, '
          f'P_cond {base_1atm["P_cond_MW"]:.0f} MW (cota ideal: {pf1["dm_star_kgps"]:.2f} kg/s, {pf1["P_cond_MW"]:.0f} MW)\n')
    print(f'{"plataforma":<22}{"carga":>6}{"h":>4}{"ṁ*":>7}{"B kg":>7}{"P_th MW":>9}{"aire":>7}{"agua":>7}'
          f'{"P_cond":>8}{"enfr t/s":>9}{"x comb":>7}{"UA MW/K":>8}{"kW/kg":>7}')
    for f in filas:
        print(f'{f["plataforma"]:<22}{f["carga_t"]:>6g}{f["h_m"]:>4g}{f["dm_star_kgps"]:>7.2f}{f["H2_kg"]:>7.0f}'
              f'{f["P_termica_MW"]:>9.0f}{f["aire_comb_kgps"]:>7.0f}{f["agua_kgps"]:>7.1f}{f["P_cond_MW"]:>8.0f}'
              f'{f["aire_enfr_tps"]:>9.1f}{f["ratio_aire_enfr_comb"]:>7.0f}{f["UA_MWK"]:>8.1f}{f["P_especifica_kWkg"]:>7.0f}')


def tabla(filas):
    """Tabla de requisitos mínimos (resultados/tablas/)."""
    def n(x, dec):
        s = f'{x:,.{dec}f}' if abs(x) >= 10000 else f'{x:.{dec}f}'
        return s.replace(',', '{,}')
    L = [r'\begin{tabular}{llrrrrrrrrr}', r'\toprule',
         r'Plataforma & Carga & $h$ & $\dot m^*$ & $B$ & $P_\text{th}$ & Aire & $P_\text{cond}$ '
         r'& Aire enfr. & $UA$ & $P_\text{cond}/m$ \\',
         r' & [t] & [m] & [kg/s] & [kg] & [MW] & [kg/s] & [MW] & [t/s] & [MW/K] & [kW/kg] \\', r'\midrule']
    previo = None
    for f in filas:
        if previo is not None and f['plataforma'] != previo:
            L.append(r'\addlinespace')
        nombre = f['plataforma'] if f['plataforma'] != previo else ''
        previo = f['plataforma']
        L.append(f"{nombre} & {f['carga_t']:g} & {f['h_m']:g} & {n(f['dm_star_kgps'], 2)} & {n(f['H2_kg'], 0)} & "
                 f"{n(f['P_termica_MW'], 0)} & {n(f['aire_comb_kgps'], 0)} & {n(f['P_cond_MW'], 0)} & "
                 f"{n(f['aire_enfr_tps'], 1)} & {n(f['UA_MWK'], 1)} & {n(f['P_especifica_kWkg'], 0)} \\\\")
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('exp06_requisitos_minimos'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/exp06_requisitos_minimos.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')



def tabla_pasos(pasos):
    """Paso máximo y duración de la entrega según la potencia de condensación instalada."""
    enc = ' & '.join(rf'\multicolumn{{3}}{{c}}{{{P:g}~MW}}' for P in P_INSTALADA_MW)
    sub = ' & '.join(['paso [kg] & pasos & [min]'] * len(P_INSTALADA_MW))
    L = [r'\begin{tabular}{llr' + 'rrr' * len(P_INSTALADA_MW) + '}', r'\toprule',
         r' &  &  & ' + enc + r' \\',
         r'Plataforma & Carga [t] & $h$ [m] & ' + sub + r' \\', r'\midrule']
    claves = sorted({(r['plataforma'], r['h_m']) for r in pasos},
                    key=lambda k: (list(p.PLATAFORMAS).index(k[0]), k[1]))
    previo = None
    for nombre, h in claves:
        fs = [r for r in pasos if r['plataforma'] == nombre and r['h_m'] == h]
        celdas = ' & '.join(f"{r['paso_max_kg']:.0f} & {r['n_pasos']} & {r['t_total_min']:.0f}" for r in fs)
        L.append(f"{nombre if nombre != previo else ''} & {fs[0]['carga_t']:g} & {h:g} & {celdas}" + r' \\')
        previo = nombre
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('exp06_pasos'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/exp06_requisitos_minimos.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
