# -*- coding: utf-8 -*-
# [EN] Tables of the article and of the supplementary material, in LaTeX.
# =============================================================================
# tablas_articulo.py — tablas del ARTÍCULO (paper 1), en inglés, generadas desde
# los resultados del pipeline. Ninguna cifra se escribe a mano: cada celda sale de
# parametros.py o de resultados/datos/. El .tex del artículo hace \input de
# resultados/tablas/art_*.tex (aplanar.py las pega en el envío).
#
# Corre después de los experimentos (lee sus datos) y antes de verificacion.py.
# =============================================================================

import io
import sys
import csv
import parametros as p
import rutas
import minimax as mm
import inercia

NOMBRE = {'Hind': 'LZ~129 Hindenburg', 'FW': 'Flying Whales LCA60T', 'ZR': 'ZR-3 Los Angeles', 'ZNT': 'Zeppelin~NT'}
ORDEN = ('Hind', 'FW', 'ZR', 'ZNT')


def n(x, dec):
    """Mismo redondeo que exportar_numeros.py: separador de miles {,} desde 10 000."""
    if x != x:
        return '--'
    s = f'{x:,.{dec}f}' if abs(x) >= 10000 else f'{x:.{dec}f}'
    return s.replace(',', '{,}')


def guardar(nombre, lineas):
    with open(rutas.tabla(nombre), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/tablas_articulo.py. NO EDITAR A MANO.\n' + '\n'.join(lineas) + '\n')


def plataformas():
    """Tabla 2: cada nave por su casco (volumen total, esbeltez, k) y su peso vacío publicado; la inercia es
    (1 + k)·rho_aire·V_casco (inercia.py). Peso vacío no publicado (NT): derivado, marcado con a; carga no
    publicada (NT, "reserved" en el TCDS; M-3 de la auditoría 3): la liberación de referencia, marcada con b."""
    L = [r'\begin{tabular}{@{} l r r r r r r r r @{}}', r'\toprule',
         r'\textbf{Airship} & \makecell{\boldmath$V_\text{hull}$\\\textbf{(m\textsuperscript{3})}} & '
         # ANTES (revisión 2026-09-30, R-02: L es también la pesadez neta de la Ec. 3; la longitud del casco pasa a \ell):
         # r'\boldmath$L/D$ & \boldmath$k$ & '
         r'\boldmath$\ell/D$ & \boldmath$k$ & '
         r'\makecell{\textbf{Empty mass}\\\textbf{(t)}} & \makecell{\textbf{Max payload}\\\textbf{(t)}} & '
         r'\makecell{\boldmath$A$\\\textbf{(m\textsuperscript{2})}} & '
         r'\makecell{\boldmath$M_I$\\\textbf{(t)}} & \boldmath$\varepsilon$ \\', r'\midrule']
    for c in ORDEN:
        nom, d = next((k, v) for k, v in p.PLATAFORMAS.items() if v['corto'] == c)
        r = mm.resolver(d, p.m_carga, p.h_max_precision, p.D_lim_oper)
        vacio = n(d['m_vacio'] / 1e3, 1) + ('' if d['m_vacio_publicado'] else r'\textsuperscript{a}')
        L.append(f"{NOMBRE[c]} & {n(d['V_casco'], 0)} & {n(d['esbeltez'], 1)} & {n(d['k'], 2)} & {vacio} & "
                 f"{n(d['carga_max'] / 1e3, 1)}{'' if d['carga_publicada'] else chr(92) + 'textsuperscript{b}'} & "
                 f"{n(d['Area'], 0)} & {n(r['M_inercial'] / 1e3, 1)} & "
                 f"{n(r['epsilon'], 3)}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    guardar('art_plataformas', L)
    # revisión 2026-10-08 (punto B): la copia del artículo va sin epsilon, que no entra en ningún cálculo;
    # art_plataformas queda igual para el del 2026-10-01
    # (autor, 2026-10-08) encabezados con símbolo y la unidad debajo: masa vacía M_str (como en la Ec. 1) y carga
    # máxima m_load,max (m_load solo es la carga que se suelta, 1 t en el caso base)
    simbolos = ((r'\makecell{\textbf{Empty mass}\\\textbf{(t)}}', r'\makecell{\boldmath$M_\text{str}$\\\textbf{(t)}}'),
                (r'\makecell{\textbf{Max payload}\\\textbf{(t)}}',
                 r'\makecell{\boldmath$m_\text{load,max}$\\\textbf{(t)}}'),
                (r' & \boldmath$\varepsilon$ \\', r' \\'))
    cab = []
    for x in L[1:4]:
        for viejo, nuevo in simbolos:
            x = x.replace(viejo, nuevo)
        cab.append(x)
    R = [r'\begin{tabular}{@{} l r r r r r r r @{}}'] + cab
    for x in L[4:-2]:
        celdas = x[:-len(r' \\')].split(' & ')
        # (autor, 2026-10-08) masa vacía no publicada (NT): sin la cifra derivada, solo la marca; el supuesto (casco lleno
        # de H2) va en la nota a del pie
        if celdas[4].endswith(r'\textsuperscript{a}'):
            celdas[4] = r'--\textsuperscript{a}'
        R.append(' & '.join(celdas[:-1]) + r' \\')
    R += [r'\bottomrule', r'\end{tabular}']
    assert all('varepsilon' not in x for x in R) and len(R) == len(L)
    guardar('art_plataformas_rev', R)


def cota_ideal():
    m = rutas.leer_datos('minimax_plataformas')['single_1t']
    # revisión 2026-10-01: dos columnas con el control reactivo puro (reactivo_vs_limpia.py): la tasa instalada que
    # necesita un controlador que responde solo al soltar para no subir más de 1 m, o más que la quema limpia
    rv = rutas.leer_datos('reactivo_vs_limpia')['naves']
    L = [r'\begin{tabular}{@{} l r r r r r @{}}', r'\toprule',
         r' & \multicolumn{3}{c}{\textbf{Preconditioned burn}} & \multicolumn{2}{c}{\textbf{Reactive only}} \\',
         r'\cmidrule(lr){2-4}\cmidrule(l){5-6}',
         r'\textbf{Airship} & \makecell{\boldmath$\dot m^*_{H_2}$\\\textbf{(kg/s)}} & '
         r'\makecell{\boldmath$t_\text{anticip}$\\\textbf{(s)}} & \makecell{\boldmath$|\Delta h|_\text{peak}$\\\textbf{(m)}} & '
         r'\makecell{\textbf{1-m rise}\\\textbf{(kg/s)}} & \makecell{\textbf{Same }\boldmath$|\Delta h|$\\\textbf{(kg/s)}} \\',
         r'\midrule']
    for c in ORDEN:
        nom = next(k for k, v in p.PLATAFORMAS.items() if v['corto'] == c)
        x = m[nom]
        r = rv[nom]['reactivo_sim']
        L.append(f"{NOMBRE[c]} & {n(x['dm_opt'], 2)} & {n(x['t_anticipo'], 2)} & "
                 f"{n(x['dh_peak_m'], 2)} & {n(r['precision'], 1)} & {n(r['igual'], 0)}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    guardar('art_cota_ideal', L)
    # revisión 2026-10-06 (recorte): la copia del artículo lleva solo la quema con precondicionamiento; las columnas
    # del control reactivo pasan al suplementario (art_reactivo). art_cota_ideal queda igual para el del 2026-10-01.
    L2 = [r'\begin{tabular}{@{} l r r r @{}}', r'\toprule',
          r'\textbf{Airship} & \makecell{\boldmath$\dot m^*_{H_2}$\\\textbf{(kg/s)}} & '
          r'\makecell{\boldmath$t_\text{anticip}$\\\textbf{(s)}} & \makecell{\boldmath$|\Delta h|_\text{peak}$\\\textbf{(m)}} \\',
          r'\midrule']
    L3 = [r'\begin{tabular}{@{} l r r r r @{}}', r'\toprule',
          r' & \multicolumn{2}{c}{\textbf{Installed rate (kg/s)}} & \multicolumn{2}{c}{\textbf{Times }\boldmath$\dot m^*_{H_2}$} \\',
          r'\cmidrule(lr){2-3}\cmidrule(l){4-5}',
          r'\textbf{Airship} & \textbf{1-m rise} & \textbf{Same }\boldmath$|\Delta h|$ & \textbf{1-m rise} & '
          r'\textbf{Same }\boldmath$|\Delta h|$ \\', r'\midrule']
    for c in ORDEN:
        nom = next(k for k, v in p.PLATAFORMAS.items() if v['corto'] == c)
        x, r = m[nom], rv[nom]['reactivo_sim']
        L2.append(f"{NOMBRE[c]} & {n(x['dm_opt'], 2)} & {n(x['t_anticipo'], 2)} & {n(x['dh_peak_m'], 2)}" + r' \\')
        L3.append(f"{NOMBRE[c]} & {n(r['precision'], 1)} & {n(r['igual'], 0)} & {n(r['precision'] / x['dm_opt'], 1)} & "
                  f"{n(r['igual'] / x['dm_opt'], 0)}" + r' \\')
    L2 += [r'\bottomrule', r'\end{tabular}']
    L3 += [r'\bottomrule', r'\end{tabular}']
    guardar('art_cota_ideal_rev', L2)
    guardar('art_reactivo', L3)
    # PROTOTIPO 2026-10-07 (pedido del autor, pendiente H): tabla 3 enriquecida para la copia, todavía sin usar en el
    # .tex. Sale la columna |dh|_peak (0.20 m en todas: (1+s) D_lim, va al pie) y entran, por nave: el tiempo en que la
    # quema limpia procesa el escalón, m/(L dm*); la velocidad de bajada al soltar, 3 D_lim/t_anticip; la tasa que
    # necesita el controlador sin quema previa, en kg/s (el autor prefiere el valor bruto; la razón con dm* va al pie),
    # para 1 m y para 5 m; y el tiempo en que los tanques
    # devuelven el H2 quemado antes de soltar si la liberación falla: L dm* t_anticip / (rho_a/rho_H2) al ritmo dm*.
    import minimax as _mm
    # (autor, 2026-10-08) sin la columna "Release fails / Tanks return (s)": solo vale para la quema limpia (el control
    # reactivo no quema antes de soltar) y confundía junto al control reactivo; el dato va al pie (\valAbortoS)
    # (autor, 2026-10-15) entra P_rej = ṁ* q, el rechazo de calor en tiempo real al retener phi_objetivo del agua (el
    # mismo q que la sección 5.1); el del LCA60T debe ser el de \valReqPcondBaseUno (exp06, base_1t_h1)
    r6t = rutas.leer_datos('exp06_requisitos_minimos')
    q6 = r6t['q_MJkg']
    L4 = [r'\begin{tabular}{@{} l r r r r r r r @{}}', r'\toprule',
          # (autor, 2026-10-08) antes "Controller only, no burn before release (kg/s)" y "Sink at release (m/s)"
          r' & \multicolumn{5}{c}{\textbf{Clean burn}} & \multicolumn{2}{c}{\makecell{\textbf{Reactive control}\\\textbf{(kg/s)}}} \\',
          r'\cmidrule(lr){2-6}\cmidrule(l){7-8}',
          r'\textbf{Airship} & \makecell{\boldmath$\dot m^*_{H_2}$\\\textbf{(kg/s)}} & '
          r'\makecell{\boldmath$P_\text{rej}$\\\textbf{(MW)}} & '
          r'\makecell{\boldmath$t_\text{anticip}$\\\textbf{(s)}} & \makecell{\boldmath$t_\text{burn}$\\\textbf{(s)}} & '
          r'\makecell{\textbf{Sink rate}\\\textbf{(m/s)}} & \makecell{\textbf{Rise}\\\boldmath$\leq$\textbf{ 1 m}} & '
          r'\makecell{\textbf{Rise}\\\boldmath$\leq$\textbf{ 5 m}} \\', r'\midrule']
    for c in ORDEN:
        nom = next(k for k, v in p.PLATAFORMAS.items() if v['corto'] == c)
        x, r = m[nom], rv[nom]
        dm, ta = x['dm_opt'], x['t_anticipo']
        t_paso = 1000.0 / (_mm.L_eff * dm)
        v0 = 3.0 * p.D_lim_oper / ta
        p_rej = dm * q6
        if c == 'FW' and abs(p_rej - r6t['base_1t_h1']['P_cond_MW']) > 0.01 * r6t['base_1t_h1']['P_cond_MW']:
            sys.exit(f"tabla 3: P_rej del LCA60T ({p_rej:.1f} MW) distinto del de la sección 5.1 "
                     f"({r6t['base_1t_h1']['P_cond_MW']:.1f} MW)")
        L4.append(f"{NOMBRE[c]} & {n(dm, 2)} & {n(p_rej, 0)} & {n(ta, 1)} & {n(t_paso, 1)} & {n(v0, 3)} & "
                  f"{n(r['reactivo_sim']['precision'], 1)} & {n(r['reactivo_sim']['seguridad'], 1)}" + r' \\')
    L4 += [r'\bottomrule', r'\end{tabular}']
    guardar('art_cota_ideal_rev2', L4)


def requisitos():
    with open(rutas.dato('exp06_requisitos_minimos', 'csv'), encoding='utf-8') as f:
        filas = list(csv.DictReader(f))
    sel = []
    for c in ORDEN:
        nom, d = next((k, v) for k, v in p.PLATAFORMAS.items() if v['corto'] == c)
        cap = min(d['carga_max'] / 1000.0, 30.0)
        for carga, h in ((1.0, 1.0), (1.0, 5.0)) + (((cap, 5.0),) if cap > 1.0 else ()):
            r = next(x for x in filas if x['corto'] == c and abs(float(x['carga_t']) - carga) < 1e-3
                     and float(x['h_m']) == h)          # el CSV va redondeado (capacidad del LZ 129: 19.046 t)
            sel.append((c, r))
    L = [r'\begin{tabular}{@{} l r r r r r r r @{}}', r'\toprule',
         r'\textbf{Airship} & \makecell{\boldmath$m_\text{load}$\\\textbf{(t)}} & '
         r'\makecell{\boldmath$h_\text{max}$\\\textbf{(m)}} & \makecell{\boldmath$\dot m^*_{H_2}$\\\textbf{(kg/s)}} & '
         r'\makecell{\textbf{Thermal}\\\textbf{(MW)}} & \makecell{\boldmath$P_\text{rej}$\\\textbf{(MW)}} & '
         r'\makecell{\textbf{Cooling air}\\\textbf{(t/s)}} & \makecell{\boldmath$Q$\\\textbf{(GJ)}} \\', r'\midrule']
    previo = None
    for c, r in sel:
        if previo is not None and c != previo:
            L.append(r'\addlinespace')
        nom = NOMBRE[c] if c != previo else ''
        previo = c
        L.append(f"{nom} & {n(float(r['carga_t']), 1 if float(r['carga_t']) % 1 else 0)} & {n(float(r['h_m']), 0)} & "
                 f"{n(float(r['dm_star_kgps']), 2)} & {n(float(r['P_termica_MW']), 0)} & {n(float(r['P_cond_MW']), 0)} & "
                 f"{n(float(r['aire_enfr_tps']), 1)} & {n(float(r['Q_GJ']), 1)}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    guardar('art_requisitos', L)


def escalones():
    e = rutas.leer_datos('exp10_escalones_condensador')
    casos = ('ideal', 'ventear_tau0s', 'ventear_tau5s', 'ventear_tau30s')
    L = [r'\begin{tabular}{@{} l r r r r r r r @{}}', r'\toprule',
         r'& & & & \multicolumn{4}{c}{\textbf{H\textsubscript{2} burned per tonne released (kg/t)}} \\',
         r'\cmidrule(l){5-8}',
         r'\textbf{Airship} & \makecell{\boldmath$P_\text{inst}$\\\textbf{(MW)}} & \makecell{\textbf{Step}\\\textbf{(kg)}} & '
         r'\makecell{\textbf{Burn}\\\textbf{(s)}} & \textbf{Ideal} & \boldmath$\tau = 0$ & \boldmath$\tau = 5$\textbf{~s} & '
         r'\boldmath$\tau = 30$\textbf{~s} \\', r'\midrule']
    previo = None
    for c in ORDEN:
        for P in (10.0, 100.0):
            xs = {x['caso']: x for x in e['casos'] if x['corto'] == c and x['P_instalada_MW'] == P}
            if previo is not None and c != previo:
                L.append(r'\addlinespace')
            nom = NOMBRE[c] if c != previo else ''
            previo = c
            celdas = []
            for k in casos:
                x = xs[k]
                celdas.append(n(x['H2_por_t_kg'], 0) + ('' if x['cumple_h_max'] else r'$^\dagger$'))
            ide = xs['ideal']
            L.append(f"{nom} & {n(P, 0)} & {n(ide['escalon_kg'], 0)} & {n(ide['dur_quema_s'], 0)} & " + ' & '.join(celdas)
                     + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    guardar('art_escalones', L)


PASOS_T = (1.0, 2.0, 5.0, 10.0, 20.0)


def potencia_escalon():
    """P_rej = ṁ* q para liberaciones de 1-20 t (cota ideal, umbral de seguridad), hasta la carga de cada nave."""
    from exp05_condensacion_tiempo_real import q_rech
    q = q_rech(0.0, p.phi_objetivo)
    L = [r'\begin{tabular}{@{} l ' + 'r ' * len(PASOS_T) + '@{}}', r'\toprule',
         r'& \multicolumn{' + str(len(PASOS_T)) + r'}{c}{\textbf{Release (t)}} \\',
         r'\cmidrule(l){2-' + str(len(PASOS_T) + 1) + '}',
         r'\textbf{Airship} & ' + ' & '.join(f'\\textbf{{{m:g}}}' for m in PASOS_T) + r' \\', r'\midrule']
    datos = {}
    for c in ORDEN:
        nom, d = next((k, v) for k, v in p.PLATAFORMAS.items() if v['corto'] == c)
        celdas = []
        for m in PASOS_T:
            if m * 1000.0 > d['carga_max']:
                celdas.append('--')
                continue
            r = mm.resolver(d, m * 1000.0, p.h_max_seguridad, p.D_lim_oper)   # h_max no dimensiona
            P = r['dm_star'] * q / 1e6
            datos[(c, m)] = P
            celdas.append(n(P, 0))
        L.append(f'{NOMBRE[c]} & ' + ' & '.join(celdas) + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    guardar('art_potencia_escalon', L)
    return datos


AMBIENTES_C = (15.0, 20.0, 25.0, 30.0, 35.0, 40.0)


def ambiente():
    """Fracción que condensa y H2 por tonelada según la temperatura exterior (1 atm, salida = exterior + 10 °C)."""
    from exp05_condensacion_tiempo_real import phi_max
    k_rho = p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1)) - 1.0
    acerc = p.T_salida_cond - p.T_aire_exterior
    fil = []
    for Ta in AMBIENTES_C:
        phi = phi_max(Ta + acerc)
        fil.append((Ta, Ta + acerc, phi, 1000.0 / (k_rho + (1.0 + p.ratio_O2_H2) * phi)))
    L = [r'\begin{tabular}{@{} l ' + 'r ' * len(fil) + '@{}}', r'\toprule',
         r'Outside air (\textdegree C) & ' + ' & '.join(n(x[0], 0) for x in fil) + r' \\',
         r'Exit of the heat-rejection unit (\textdegree C) & ' + ' & '.join(n(x[1], 0) for x in fil) + r' \\', r'\midrule',
         r'Water retained, $\phi_\text{max}$ & ' + ' & '.join(n(x[2], 2) for x in fil) + r' \\',
         r'H\textsubscript{2} per tonne (kg) & ' + ' & '.join(n(x[3], 1) for x in fil) + r' \\',
         r'\bottomrule', r'\end{tabular}']
    guardar('art_ambiente', L)
    return fil


def helio():
    e = rutas.leer_datos('exp12_helio')
    nombres = ['H\\textsubscript{2} burned, water retained ($\\phi = PHI$)',
               'H\\textsubscript{2} burned, water retained at 1~atm ($\\phi = PHI$)',
               'H\\textsubscript{2} vented', 'He compressed on board (COSH, 20~MPa)', 'He vented']
    L = [r'\begin{tabular}{@{} l r r r r r @{}}', r'\toprule',
         r'\textbf{Route} & \makecell{\boldmath$L$\\\textbf{(kg/kg)}} & \makecell{\textbf{Gas}\\\textbf{(kg/t)}} & '
         r'\makecell{\textbf{Volume}\\\textbf{(m\textsuperscript{3}/t)}} & \makecell{\textbf{Heat}\\\textbf{(GJ/t)}} & '
         r'\makecell{\textbf{Lost}\\\textbf{(kg/t)}} \\', r'\midrule']
    phis = [p.phi_objetivo, 1.0 - e['rutas'][1]['perdida_frac']]
    for i, r in enumerate(e['rutas']):
        nom = nombres[i].replace("PHI", f"{phis[i]:.2f}") if i < 2 else nombres[i]
        L.append(f"{nom} & {n(r['L'], 2)} & {n(r['kg_por_t'], 1)} & {n(r['m3_por_t'], 0)} & "
                 f"{n(r['calor_GJ_por_t'], 2)} & {n(r['perdido_kg_por_t'], 1)}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    guardar('art_helio', L)

    # ── Revisión 2026-10-06: la comparación en DOS ETAPAS para la copia (art_helio_rev; art_helio queda igual) ──
    # Rutas en columnas, magnitudes en filas. Liberación: L, calor, potencia de golpe (misma forma cerrada: el tiempo no
    # depende del gas) y en el tiempo de recuperación del ciclo con el rechazo de referencia (freno), potencia neta
    # hacia la planta en ese tiempo (+ genera, - consume) y cómo se guarda lo que se quita. Devolución de la
    # sustentación (corrección del autor, paso [5]): gas sacado de los tanques a bordo (el ciclo, los B kg que quemó;
    # COSH, todo su He), la parte que devuelve al instante, el compuesto de esos tanques (mismo escalado P·V para los
    # dos, cota inferior), con qué vuelve el resto, la energía de esta etapa y el gas perdido.
    import minimax as mm
    from exp05_condensacion_tiempo_real import q_rech
    from exp02_quema_centrada import L_ef
    from exp16_balance_neto import e_comp
    from balance_masa import f_caida_presion
    t_rec = (1000.0 / L_ef(p.phi_objetivo)) * q_rech(0.0, p.phi_objetivo) / p.P_ref_masa
    R = e['rutas']
    dm_h2 = e['base']['H2']['dm_kgps']
    ref = f"{e['P_ref_MPa']:.1f}"
    alm = next(a for a in e['almacen'] if abs(a['P_MPa'] - e['P_ref_MPa']) < 1e-9)
    golpe = [e['base']['H2']['potencia_MW'], dm_h2 * mm.L_eff / R[1]['L'] * R[1]['q_Jkg'] / 1e6, None,
             e['base']['He']['potencia_MW'][ref], None]
    f_min, f_max = min(f_caida_presion()), max(p.perdida_rad)
    ec = e_comp()
    col = []
    for i, r in enumerate(R):
        quema = r['gas'] == 'H2' and r['q_Jkg'] > 0.0
        c = {'L': n(r['L'], 2), 'calor': n(r['calor_GJ_por_t'], 2),
             'golpe': '--' if golpe[i] is None else n(golpe[i], 0),
             'rec': '--' if r['calor_GJ_por_t'] == 0.0 else n(r['calor_GJ_por_t'] * 1e9 / t_rec / 1e6, 1)}
        if quema:
            dm = r['kg_por_t'] / t_rec
            # (2026-10-12) a 14 atm (phi 0.99) el motor también paga el compresor del escape, que sale como calor
            Wc = rutas.leer_datos('exp06_requisitos_minimos')['W_comp_MJkg'] * 1e6 if '1 atm' not in r['ruta'] else 0.0
            neta = [eta * p.Q_comb_H2 * dm - f * (r['q_Jkg'] - eta * p.Q_comb_H2 + Wc) * dm - Wc * dm
                    for eta in p.eta_motor_H2 for f in (f_min, f_max)]
            c['neta'] = '+' + n(min(neta) / 1e6, 1) + '--' + n(max(neta) / 1e6, 1)
            c['guarda'] = r'water, ' + n(r['agua_m3_por_t'], 2) + r'~m\textsuperscript{3}'
            Rg = r['kg_por_t'] * (1.0 - r['perdida_frac'])
            c['energia'] = n(Rg * (p.E_H2_electrolisis + min(ec)) / 1e3, 1) + '--' + n(Rg * (p.E_H2_electrolisis + max(ec)) / 1e3, 1)
            c['resto'] = r'\makecell{electrolysis,\\O\textsubscript{2} to the air}'
            c['t_resto'] = n(r['t_resto_h'], 0)
        elif r['gas'] == 'He' and r['q_Jkg'] > 0.0:
            c['neta'] = r'$-$' + n(r['calor_GJ_por_t'] * 1e9 / t_rec / 1e6, 1)
            c['guarda'] = r'gas, ' + n(alm['m3_almacen_por_t'], 1) + r'~m\textsuperscript{3}'
            c['energia'] = '0'
            c['resto'] = '--'
            c['t_resto'] = '--'
        else:
            c.update({'neta': '0', 'guarda': 'none', 'energia': '--', 'resto': 'ground resupply', 't_resto': '--'})
        if r['reserva_kg'] is None:
            c.update({'reserva': '--', 'instante': '--', 'tanque': '--'})
        else:
            c['reserva'] = n(r['reserva_kg'], 1)
            c['instante'] = n(r['al_instante_kg'] / 10.0, 0)
            c['tanque'] = n(r['reserva_compuesto_t'], 2)
        c['perdido'] = n(r['perdido_kg_por_t'], 1)
        col.append(c)
    # Auditoría 2 (M1): H2 comprimido a bordo, a la presión de los tanques del ciclo, como tercera columna
    hc = e['h2_comprimido']
    col.insert(2, {'L': n(hc['L'], 2), 'calor': n(hc['calor_GJ_por_t'], 2), 'golpe': n(hc['potencia_MW'], 0),
                   'rec': n(hc['calor_GJ_por_t'] * 1e9 / t_rec / 1e6, 1),
                   'neta': r'$-$' + n(hc['calor_GJ_por_t'] * 1e9 / t_rec / 1e6, 1),
                   'guarda': r'gas, ' + n(hc['m3_almacen_por_t'], 1) + r'~m\textsuperscript{3}',
                   'reserva': n(hc['kg_por_t'], 1), 'instante': n(hc['al_instante_kg'] / 10.0, 0),
                   'tanque': n(hc['reserva_compuesto_t'], 2), 'resto': '--', 't_resto': '--', 'energia': '0',
                   'perdido': n(hc['perdido_kg_por_t'], 1)})
    cab = [r'\makecell{H\textsubscript{2} burned\\($\phi = ' + f"{p.phi_objetivo:.2f}" + r'$)}',
           r'\makecell{H\textsubscript{2} burned\\at 1~atm ($\phi = ' + f"{1.0 - R[1]['perdida_frac']:.2f}" + r'$)}',
           r'\makecell{H\textsubscript{2} compressed\\(' + f"{hc['P_MPa']:.0f}" + r'~MPa)}',
           r'\makecell{H\textsubscript{2}\\vented}',
           r'\makecell{He compressed\\(COSH, ' + f"{e['P_ref_MPa']:.0f}" + r'~MPa)}', r'\makecell{He\\vented}']
    fila = lambda nom, k: nom + ' & ' + ' & '.join(c[k] for c in col) + r' \\'
    L2 = [r'\begin{tabular}{@{} l c c c c c c @{}}', r'\toprule', r'Per tonne released & ' + ' & '.join(cab) + r' \\',
          r'\midrule', r'\multicolumn{7}{@{}l}{\textit{Release}} \\',
          fila(r'Heaviness per kg of gas, $L$ (kg/kg)', 'L'),
          fila(r'Heat to reject (GJ)', 'calor'),
          fila(r'Power, at once (MW)', 'golpe'),
          fila(r'Power, in ' + f"{t_rec / 60:.1f}" + r'~min with a hold-down device (MW)', 'rec'),
          fila(r'Net power to the plant in that time (MW)', 'neta'),
          fila(r'Lift removed, stored as', 'guarda'),
          r'\addlinespace', r'\multicolumn{7}{@{}l}{\textit{Restoring the lift}} \\',
          fila(r'Gas returned from onboard tanks (kg)', 'reserva'),
          fila(r'Lift restored at once (\%)', 'instante'),
          fila(r'Composite of those tanks, at least (t)', 'tanque'),
          fila(r'Rest of the lift restored by', 'resto'),
          fila(r'Time to restore the rest (h)', 't_resto'),
          fila(r'Energy in this stage (MWh)', 'energia'),
          fila(r'Gas lost (kg)', 'perdido'),
          r'\bottomrule', r'\end{tabular}']
    guardar('art_helio_rev', L2)

    # PROTOTIPO 2026-10-07 (pedido del autor, pendiente I): la misma tabla, más legible, todavía sin usar en el .tex.
    # Las unidades van en las celdas, no en las etiquetas; columnas agrupadas por gas; sale la fila "Power in 10.6 min"
    # (es el calor entre 634 s, y para el ciclo 10 MW por definición: la lleva la fila de potencia neta); entran el gas
    # procesado por tonelada y, para el H2 venteado, la electrólisis del H2 que se repone en tierra.
    t_nat = 1000.0 / (mm.L_eff * dm_h2)                      # el mismo para todas las rutas (Sección 8)

    def u(x, unidad):
        return x if x in ('--', '0') else x + '~' + unidad

    def cero(x):
        return '0' if x in ('0.00', '0.0', '0') else x

    def pila(a, b):
        return r'\makecell{' + a + r'\\' + b + '}'

    orden = [R[0], R[1], 'hc', R[2], R[3], R[4]]
    c2 = []
    for i, r in enumerate(orden):
        c = col[i]
        hc_ = r == 'hc'
        kg = hc['kg_por_t'] if hc_ else r['kg_por_t']
        x = {'L': c['L'], 'gas': u(n(kg, 1), 'kg'), 'calor': u(cero(c['calor']), 'GJ'),
             'golpe': u('0' if c['golpe'] == '--' else c['golpe'], 'MW'), 'neta': u(c['neta'], 'MW'),
             'reserva': u(c['reserva'], 'kg'), 'instante': u(c['instante'], r'\%'), 'tanque': u(c['tanque'], 't'),
             't_resto': u(c['t_resto'], 'h'), 'energia': u(cero(c['energia']), 'MWh'),
             'perdido': u(cero(c['perdido']), 'kg')}
        if hc_:
            x['guarda'] = pila('gas', n(hc['m3_almacen_por_t'], 1) + r'~m\textsuperscript{3}')
            x['resto'] = '--'
        elif r['gas'] == 'H2' and r['q_Jkg'] > 0.0:
            x['guarda'] = pila('water', n(r['agua_m3_por_t'], 2) + r'~m\textsuperscript{3}')
            x['resto'] = pila('electrolysis,', r'O\textsubscript{2} to the air')
        elif r['gas'] == 'He' and r['q_Jkg'] > 0.0:
            x['guarda'] = pila('gas', n(alm['m3_almacen_por_t'], 1) + r'~m\textsuperscript{3}')
            x['resto'] = '--'
        else:
            x['guarda'] = '--'
            x['instante'] = '0'                              # ventear no devuelve nada de golpe (aplica y vale 0)
            x['resto'] = pila('ground', 'resupply')
            if r['gas'] == 'H2':                             # el H2 venteado se repone por electrólisis en tierra
                x['energia'] = u(n(r['perdido_kg_por_t'] * p.E_H2_electrolisis / 1e3, 1), 'MWh')
        c2.append(x)
    fila2 = lambda nom, k: nom + ' & ' + ' & '.join(c[k] for c in c2) + r' \\'
    L5 = [r'\begin{tabular}{@{} l c c c c c c @{}}', r'\toprule',
          r' & \multicolumn{4}{c}{\textbf{Hydrogen}} & \multicolumn{2}{c}{\textbf{Helium}} \\',
          r'\cmidrule(lr){2-5}\cmidrule(l){6-7}',
          r'\textbf{Per tonne released} & \makecell{\textbf{burned}\\($\phi = ' + f"{p.phi_objetivo:.2f}" + r'$)} & '
          r'\makecell{\textbf{burned, 1~atm}\\($\phi = ' + f"{1.0 - R[1]['perdida_frac']:.2f}" + r'$)} & '
          r'\makecell{\textbf{compressed}\\(' + f"{hc['P_MPa']:.0f}" + r'~MPa)} & \textbf{vented} & '
          r'\makecell{\textbf{compressed}\\(COSH, ' + f"{e['P_ref_MPa']:.0f}" + r'~MPa)} & \textbf{vented} \\',
          r'\midrule', r'\multicolumn{7}{@{}l}{\textit{Release}} \\',
          fila2(r'Heaviness per kg of gas, $L$', 'L'),
          fila2(r'Gas processed', 'gas'),
          fila2(r'Heat to reject', 'calor'),
          fila2(r'Power, at once (' + f"{t_nat:.0f}" + r'~s)', 'golpe'),
          fila2(r'\makecell[l]{Net power to the plant,\\with a hold-down device (' + f"{t_rec / 60:.1f}" + r'~min)}', 'neta'),
          fila2(r'Lift removed, stored as', 'guarda'),
          r'\addlinespace', r'\multicolumn{7}{@{}l}{\textit{Restoring the lift}} \\',
          fila2(r'Gas returned from onboard tanks', 'reserva'),
          fila2(r'Lift restored at once', 'instante'),
          fila2(r'Composite of those tanks, at least', 'tanque'),
          fila2(r'Rest of the lift restored by', 'resto'),
          fila2(r'Time to restore the rest', 't_resto'),
          fila2(r'Energy in this stage', 'energia'),
          fila2(r'Gas lost', 'perdido'),
          r'\bottomrule', r'\end{tabular}']
    guardar('art_helio_rev2', L5)


# Filas de la tabla de chequeos de datos de referencia (tabla S3 del material suplementario desde la revisión
# 2026-10-01; antes tabla 9 del artículo). exportar_numeros.py toma de aquí el desvío máximo sin la densidad del aire.
SEL_VALIDACION = [
    ('Atmosphere', 'densidad del aire en h0', r'Air density, sea-level reference at $h_0$ (kg/m\textsuperscript{3})',
     'ISA', 4),
    ('', 'altura de escala de la densidad en h0', r'Density scale height at $h_0$ (m)', 'ISA', 0),
    ('Gases', 'rho_H2 a 15 °C, 1 atm (gas ideal)', r'H\textsubscript{2} density, 15\,\textdegree C (kg/m\textsuperscript{3})', 'NIST', 5),
    ('Combustion', 'PCI del H2', r'H\textsubscript{2} lower heating value (MJ/kg)', 'CODATA', 2),
    ('', 'calor latente del agua a 25 °C', r'Latent heat of water, 25\,\textdegree C (kJ/kg)', 'NIST', 1),
    ('', 'PCS = PCI + (M_H2O/M_H2) L_v', r'HHV rebuilt as LHV $+ (M_{H_2O}/M_{H_2}) L_v$ (MJ/kg)', 'CODATA', 2),
    ('Saturation', 'p_sat a 40 °C (IF97)', r'$p_\text{sat}$ at 40\,\textdegree C, IF97 (kPa)', 'IAPWS-95', 4),
    ('', 'p_sat a 100 °C (IF97)', r'$p_\text{sat}$ at 100\,\textdegree C, IF97 (kPa)', 'IAPWS-95', 3),
    ('Electrolysis', 'energía por kg de H2', r'Electricity per kg H\textsubscript{2} (kWh/kg)', 'IEA', 1),
]


def validacion():
    v = rutas.leer_datos('val01_componentes')
    f = {x['magnitud']: x for x in v['filas']}
    # Revisión 2026-09-30, R-13: la tabla es un chequeo de los datos de entrada, no una validación. Salen las filas
    # triviales (g, fracción de O2 en el aire, O2 por kg de H2: dato = fuente por construcción) y la densidad del aire se
    # rotula como lo que es, la referencia a nivel del mar usada en h0 (una elección, no un error). ANTES:
    # sel = [('Atmosphere', 'densidad del aire en h0', r'Air density at $h_0$ (kg/m\textsuperscript{3})', 'ISA', 4),
    #        ('', 'g', r'Gravity (m/s\textsuperscript{2})', 'ISA', 5),
    #        ('', 'fracción másica de O2 en el aire', r'O\textsubscript{2} mass fraction in air', 'US 1976', 4),
    #        ('', 'O2 por kg de H2', r'O\textsubscript{2} per kg H\textsubscript{2}', 'IUPAC', 4),
    sel = SEL_VALIDACION
    L = [r'\begin{tabular}{@{} l l r r r l @{}}', r'\toprule',
         r'\textbf{Block} & \textbf{Quantity} & \textbf{Model} & \textbf{Reference} & \textbf{Dev.} & \textbf{Source} \\',
         r'\midrule']
    for bloque, clave, texto, fuente, dec in sel:
        x = f[clave]
        if bloque and len(L) > 4:
            L.append(r'\addlinespace')
        L.append(f"{bloque} & {texto} & {n(x['modelo'], dec)} & {n(x['referencia'], dec)} & "
                 f"{100 * x['desvio_rel']:+.2f}\\,\\% & {fuente}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    guardar('art_validacion', L)


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    for fn in (plataformas, cota_ideal, requisitos, escalones, helio, validacion, potencia_escalon, ambiente):
        fn()
        print(f'  {fn.__name__}: ok')


if __name__ == '__main__':
    main()
