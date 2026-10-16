# -*- coding: utf-8 -*-
# [EN] Credibility assessment with NASA-STD-7009B and the sensitivity table (Supplementary S4-S6, Table S4, Figure S2).
# =============================================================================
# credibilidad.py — CRED-01: evaluación de credibilidad del paper 1 con los factores de NASA-STD-7009B
# (2024), usada como guía (Apéndice D.1 de la norma: un M&S fuera de un programa de la NASA puede usarla
# como buenas prácticas). Revisión del 2026-10-01.
#
# A. Criticidad [M&S 6], con la matriz del Apéndice D (Figura 1): consecuencia de la decisión que informan
#    los resultados (seguir o no con el ciclo y qué caja negra desarrollar: costo y plazo de I+D, II-III) e
#    influencia del M&S (4, "significant": la decisión sale sobre todo del modelo, con apoyo de sistemas
#    parecidos: la pila de Zakhvatkin, el blimp de Yamasaki, la recuperación de agua de los rígidos).
# B. Niveles 0-4 de los 5 factores de CAPACIDAD [M&S 48] y los 6 de RESULTADOS [M&S 31] (Apéndice E,
#    tablas 6-11 y 13-20), con su evidencia y su brecha, y el umbral que fija este estudio [M&S 43e]. Los
#    niveles son JUICIOS de los autores, no salidas del modelo: viven aquí, en un solo lugar, para que el
#    texto (macros \valCred...), la tabla y la figura no se contradigan. El autor debe confirmarlos.
#    Regla de la norma: cada nivel es un escalón (sin crédito parcial) y vale solo si se cumplen TODAS sus
#    condiciones; los factores no se combinan ni se comparan entre sí (E.2).
# C. Robustez [M&S 30]: sensibilidades locales EXACTAS de la ecuación 7 (elasticidades d ln y / d ln x del
#    flujo instalado) y el cambio de los requisitos sobre el rango de cada entrada, con su naturaleza
#    (diseño, epistémica, aleatoria). Mismas fórmulas que VAL-03 (minimax.resolver; potencia con el calor de
#    condensar phi_objetivo, sec. 5.3); el grupo Y de verificacion.py las cruza con VAL-03 y las tablas.
# Salidas: resultados/datos/credibilidad.json, resultados/tablas/art_sensibilidad.tex (tabla S4 del material
#          suplementario), resultados/figuras/credibilidad.png (figura S1). Revisión 2026-10-01, versión corta: la
#          evaluación completa va al material suplementario; el artículo cita los niveles en su sección 8.4.
# =============================================================================

import io
import sys
import math
import parametros as p
import rutas
import minimax as mm
import fuentes
from exp02_quema_centrada import L_ef
from exp05_condensacion_tiempo_real import q_rech, phi_max


# ── A. Criticidad (Apéndice D, Figura 1): filas influencia 5..1, columnas consecuencia I..V ──────────────
MATRIZ = {5: 'GYRRR', 4: 'GYYRR', 3: 'GYYYR', 2: 'GGGYY', 1: 'GGGGY'}
ROMANO = ('I', 'II', 'III', 'IV', 'V')
CONSECUENCIA = ('II', 'III')        # menor a moderada: costo y plazo de I+D; ninguna decisión de vuelo o seguridad
INFLUENCIA = 4                      # significativa: no hay datos del sistema real; sí de sistemas parecidos

# ── B. Niveles (juicio de los autores; ver la justificación de cada uno) ──────────────────────────────────
UMBRAL = 2   # [M&S 43e] el nivel en que cada factor descansa en evidencia formal y documentada, no informal;
             # elegido por los autores para una decisión de diseño conceptual en la zona amarilla
CAPACIDAD = [
    # (clave de macro, nombre en el artículo, nivel, justificación)
    ('DatosPedigri', 'Data pedigree', 2,
     'Todos los datos de desarrollo trazados a documentos formales (registro de fuentes con SHA-256 y pedigrí '
     'A-C, tabla S3 del suplementario). No llega a 3: la masa añadida viene de flujo potencial acotado con un blimp de 18 m3, un '
     'referente no representativo de un casco de 180 000 m3.'),
    ('Verificacion', 'Verification', 3,
     'Solución exacta (0.003 %), convergencia del paso y la tolerancia, forma cerrada reproducida por la '
     'simulación completa y los chequeos automáticos de punta a punta que detienen la exportación si uno falla; '
     'todos los errores importantes cumplen los criterios declarados (sec. 8). No llega a 4: no hay métodos '
     '"reliable" de estimación de error para todo el modelo (p. ej. un segundo código independiente).'),
    ('Validacion', 'Validation', 2,
     'La cota de condensación contra mediciones de un sistema representativo (pila H2-aire, Zakhvatkin): '
     'favorable con el criterio de una cota superior; la masa vertical contra ensayos de vuelo de un blimp. No '
     'llega a 3: no hay mediciones del sistema real (un dirigible soltando carga).'),
    ('RevisionDes', 'Technical review', 1,
     'Revisión interna informal de coautores que no desarrollaron el modelo (el autor lo confirmó el '
     '2026-10-16) y auditorías internas del código y las fuentes. No llega a 2: la revisión no fue formal, '
     'sancionada y con reglas explícitas (E.3.5); las auditorías con IA no son revisión por pares.'),
    ('GestionDes', 'Process/product management', 1,
     'Proceso documentado (LEEME, reglas del proyecto), parámetros en un solo archivo, cada cifra generada, '
     'control de versiones con git y una compuerta de verificación con sello de los datos (gestión de '
     'configuración formal). Queda en 1 (2026-10-16, tabla 11 de la norma): el 2 '
     'pide además un proceso con roles y requisitos formalmente documentados y aprobados, y aquí son informales.'),
]
RESULTADOS = [
    ('Uso', 'Use assessment', 1,
     'Las corridas están dentro de los límites del M&S y de su dominio de verificación, pero es el PRIMER uso '
     'del modelo: la tabla 15 de la norma lo pone en 1 por definición.'),
    ('PedigriEntrada', 'Input pedigree', 1,
     'Entradas de las naves de documentos formales (CAA 1938, NACA TM 286, TCDS de EASA, Goodyear 1975) salvo '
     'el LCA60T (sitio del fabricante, discrepancia 180 000 / 200 000 m3 declarada); incertidumbre estimada en '
     'las que dimensionan (volumen del casco, inercia, agua retenida). Queda en 1 (2026-10-16, tabla 16 de la norma): el 2 '
     'pide la incertidumbre de todas las entradas al menos estimada, y la de las naves históricas no lo está.'),
    ('Incertidumbre', 'Uncertainty characterization', 2,
     'Fuentes identificadas, cuantificadas y clasificadas (epistémica/aleatoria); propagación por cotas en los '
     'extremos medidos (las cantidades son monótonas en cada entrada). No llega a 3: no hay propagación '
     'probabilística de todas las fuentes (con dos ensayos y seis experimentos no hay distribución que muestrear, '
     'hallazgo S4).'),
    ('Robustez', 'Results robustness', 3,
     'Sensibilidades cuantitativas de todas las entradas de importancia alta de la PIRT (inercia, agua retenida, '
     'temperatura exterior, potencia de rechazo, masa de las otras cajas) y de varias de importancia media y '
     'baja, elasticidades exactas de la ecuación 7 y combinaciones (inercia x agua retenida en las esquinas; '
     'f x potencia específica de la batería). No llega a 4: no se cuantificaron muchas combinaciones.'),
    ('RevisionUso', 'Technical review', 1,
     'Revisión interna informal de los coautores (ver la de desarrollo).'),
    ('GestionUso', 'Process/product management', 1,
     'La misma gestión de configuración: cada corrida del artículo sale de correr_todo.py con su registro y su '
     'versión. Queda en 1, como la de desarrollo (tabla 20 de la norma).'),
]


def zona(consecuencia, influencia):
    return MATRIZ[influencia][ROMANO.index(consecuencia)]


# ── C. Sensibilidades de la forma cerrada ─────────────────────────────────────────────────────────────────
def requisitos(m=None, D=None, phi=1.0, fk=1.0, V_casco=None, nave=None):
    """Flujo instalado, potencia de rechazo y H2 por tonelada en la cota ideal, con las mismas fórmulas que
    VAL-03: inercia (1 + k) rho V_casco por fk, L = L(phi), calor de condensar min(phi, phi_objetivo)."""
    d = dict(p.NAVE_BASE if nave is None else nave)
    if V_casco is not None:
        d['V_casco'] = V_casco
    m = p.m_carga if m is None else m
    D = p.D_lim_oper if D is None else D
    L = L_ef(phi)
    r = mm.resolver(d, m, p.h_max_precision, D, rho_gas=mm.rho_H2, L_gas=L, k_ma=(1.0 + d['k']) * fk - 1.0)
    q = q_rech(0.0, min(phi, p.phi_objetivo))
    return {'dm_star': r['dm_star'], 'P_MW': r['dm_star'] * q / 1e6, 'H2_t': 1000.0 / L, 's': r['s'],
            'M_I': r['M_inercial']}


def elasticidades(M_I, m, D, L):
    """d ln(ṁ*)/d ln x de la ecuación 7, analíticas y numéricas (diferencias centrales en ln x)."""
    s = math.sqrt(1.0 - m / M_I)
    c = 3.0 * (1.0 - s) / (4.0 * s)
    ana = {'m': 1.5 + c, 'M_I': -0.5 - c, 'L': -1.0, 'D': -0.5}
    h = 1e-5
    f = lambda **k: math.log(mm.limpio(k.get('M_I', M_I), k.get('m', m), k.get('D', D), k.get('L', L))['dm_star'])
    base = dict(M_I=M_I, m=m, D=D, L=L)
    num = {x: (f(**{x: base[x] * math.exp(h)}) - f(**{x: base[x] * math.exp(-h)})) / (2.0 * h) for x in base}
    return {'s': s, 'analitica': ana, 'numerica': num}


def fphi_medido():
    """Recolección MEDIDA / cota de equilibrio en Zakhvatkin (VAL-02): mín. y máx. (como VAL-03)."""
    r = [x['medido_sobre_modelo'] for x in rutas.leer_datos('val02_condensacion')['experimentos']]
    return min(r), max(r)


def robustez():
    base = requisitos()
    pm = phi_max(p.T_salida_cond)
    fphi = fphi_medido()
    s = base['s']
    D_max = p.h_max_precision / (1.0 + s)                       # el mayor D_lim que cumple el umbral de 1 m
    V_alt = p.V_casco_FW_alternativo                             # la otra cifra del fabricante (FAQ, FWH26)
    filas = [
        # clave, entrada, naturaleza, rango (texto), extremos de la entrada -> requisitos
        ('m', r'Release size $m_\text{load}$', 'design', r'1--20~t',
         [requisitos(m=1000.0), requisitos(m=20000.0)]),
        ('D', r'Allowed sink $D_\text{lim}$', 'design', f'{p.D_lim_oper:.2f}--{D_max:.2f}~m',
         [requisitos(), requisitos(D=D_max)]),
        ('M_I', r'Vertical inertia $M_I$', 'epistemic', r'0.77--1.20 of potential flow',
         [requisitos(fk=1.20), requisitos(fk=0.77)]),
        ('V', r'Hull volume, LCA60T', 'epistemic',
         f'{p.NAVE_BASE["V_casco"]:,.0f}--{V_alt:,.0f}~m\\textsuperscript{{3}}'.replace(',', '{,}'),
         [requisitos(), requisitos(V_casco=V_alt)]),
        ('phi', r'Water retained $\phi$', 'epistemic',
         f'1 to {pm * fphi[0]:.2f}--{pm * fphi[1]:.2f}',
         [requisitos(), requisitos(phi=pm * fphi[1]), requisitos(phi=pm * fphi[0])]),
    ]
    return {'base': base, 'phi_max': pm, 'fphi': fphi, 'D_max': D_max, 'V_alt': V_alt, 'filas': filas}


def n(x, dec):
    s = f'{x:,.{dec}f}' if abs(x) >= 10000 else f'{x:.{dec}f}'
    return s.replace(',', '{,}')


def rango_txt(vals, dec=None):
    """Rango mín.--máx.; dec=None: 3 cifras significativas para el flujo (de 1.42 a 288 kg/s)."""
    d = (lambda x: 2 if abs(x) < 10 else (1 if abs(x) < 100 else 0)) if dec is None else (lambda x: dec)
    a, b = n(min(vals), d(min(vals))), n(max(vals), d(max(vals)))
    return a if a == b else f'{a}--{b}'


def tabla_sensibilidad(R, E):
    """Tabla S4 del suplementario (inglés): elasticidad del flujo y requisitos sobre el rango de cada entrada."""
    el = {'m': E['analitica']['m'], 'D': E['analitica']['D'], 'M_I': E['analitica']['M_I'],
          'V': E['analitica']['M_I'], 'phi': -(1.0 + p.ratio_O2_H2) / L_ef(1.0)}
    L = [r'\begin{tabular}{@{} l l l r r r r @{}}', r'\toprule',
         r'\textbf{Input} & \textbf{Nature} & \textbf{Range} & \textbf{Elasticity} & '
         r'$\boldsymbol{\dot m^*_{H_2}}$ \textbf{(kg/s)} & $\boldsymbol{P_\text{rej}}$ \textbf{(MW)} & '
         r'\textbf{H\textsubscript{2} (kg/t)} \\', r'\midrule']
    for clave, nom, nat, rng, rs in R['filas']:
        e = el[clave]
        es = (r'$-$' if e < 0 else '') + f'{abs(e):.2f}'
        L.append(f"{nom} & {nat} & {rng} & {es} & {rango_txt([r['dm_star'] for r in rs])} & "
                 f"{rango_txt([r['P_MW'] for r in rs], 0)} & {rango_txt([r['H2_t'] for r in rs], 1)}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('art_sensibilidad'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/credibilidad.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


def graficar():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401  (estilo de journal)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    AZUL, TINTA, GRIS = '#2a78d6', '#0b0b0b', '#52514e'
    grupos = (('M&S capability (development)', CAPACIDAD), ('M&S results (use)', RESULTADOS))
    fig, ejes = plt.subplots(2, 1, figsize=(6.6, 4.6), sharex=True,
                             gridspec_kw={'height_ratios': [len(CAPACIDAD), len(RESULTADOS)], 'hspace': 0.35})
    for ax, (titulo, filas) in zip(ejes, grupos):
        y = list(range(len(filas)))[::-1]
        for yi, (_, nom, nivel, _) in zip(y, filas):
            bajo = nivel < UMBRAL
            ax.barh(yi, nivel, height=0.55, color=AZUL if not bajo else 'white', edgecolor=AZUL, linewidth=1.2,
                    hatch='////' if bajo else None, zorder=2)
            ax.plot([UMBRAL, UMBRAL], [yi - 0.38, yi + 0.38], color=TINTA, lw=2.0, zorder=3)
            ax.text(nivel + 0.1, yi, f'{nivel}', va='center', ha='left', color=TINTA, fontsize=10, zorder=4)
        ax.set_yticks(y, [f[1] for f in filas])
        ax.set_title(titulo, loc='left', fontsize=11, color=TINTA)
        ax.set_xlim(0, 4.3)
        ax.set_ylim(-0.7, len(filas) - 0.3)
        ax.grid(axis='y', visible=False)
        ax.tick_params(axis='y', length=0, right=False)
    ejes[-1].set_xticks(range(5))
    ejes[-1].set_xlabel('Assessed level (NASA-STD-7009B, Appendix E)', color=GRIS)
    ley = [Patch(facecolor=AZUL, edgecolor=AZUL, label='Achieved level'),
           Patch(facecolor='white', edgecolor=AZUL, hatch='////', label='Achieved, below threshold'),
           Line2D([0], [0], color=TINTA, lw=2.0, label=f'Threshold of this study ({UMBRAL})')]
    ejes[0].legend(handles=ley, loc='lower right', fontsize=8.5)
    plot_style.guardar_paper(fig, rutas.figura('credibilidad'))
    plt.close(fig)


def main():
    # aquí y no al importar: verificacion.py importa este módulo y ya envolvió su stdout
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    zonas = [zona(c, INFLUENCIA) for c in CONSECUENCIA]
    d = p.NAVE_BASE
    base = requisitos()
    E = elasticidades(base['M_I'], p.m_carga, p.D_lim_oper, L_ef(1.0))
    E_naves = {}
    for nom, nv in p.PLATAFORMAS.items():
        r = requisitos(nave=nv)
        E_naves[nv['corto']] = elasticidades(r['M_I'], p.m_carga, p.D_lim_oper, L_ef(1.0))['analitica']
    R = robustez()
    pm = R['phi_max']
    el_phi = {'ideal': -(1.0 + p.ratio_O2_H2) * 1.0 / L_ef(1.0), 'atm': -(1.0 + p.ratio_O2_H2) * pm / L_ef(pm)}
    fu = [f for f in fuentes.FUENTES.values() if f.get('ambito') == 'paper 1']
    cuenta = {g: sum(1 for f in fu if f.get('pedigri') == g) for g in ('A', 'B', 'C')}
    datos = {
        'criticidad': {'consecuencia': CONSECUENCIA, 'influencia': INFLUENCIA, 'zonas': zonas},
        'umbral': UMBRAL,
        'capacidad': [{'clave': k, 'factor': f, 'nivel': nv, 'justificacion': j} for k, f, nv, j in CAPACIDAD],
        'resultados': [{'clave': k, 'factor': f, 'nivel': nv, 'justificacion': j} for k, f, nv, j in RESULTADOS],
        'elasticidades_base': E, 'elasticidades_naves': E_naves, 'elasticidad_phi': el_phi,
        'robustez': {'base': R['base'], 'phi_max': pm, 'fphi': R['fphi'], 'D_max': R['D_max'], 'V_alt': R['V_alt'],
                     'filas': [{'clave': c, 'entrada': nom, 'naturaleza': nat, 'rango': rng, 'requisitos': rs}
                               for c, nom, nat, rng, rs in R['filas']]},
        'fuentes_paper1': {'total': len(fu), **cuenta},
    }
    rutas.guardar_datos('credibilidad', datos)
    tabla_sensibilidad(R, E)
    graficar()
    print(f"A. criticidad: consecuencia {'-'.join(CONSECUENCIA)}, influencia {INFLUENCIA} -> zona {'/'.join(zonas)}")
    print(f"B. umbral {UMBRAL}; capacidad " + ', '.join(f'{f} {nv}' for _, f, nv, _ in CAPACIDAD))
    print('   resultados ' + ', '.join(f'{f} {nv}' for _, f, nv, _ in RESULTADOS))
    print(f"   fuentes del paper 1: {len(fu)} ({cuenta})")
    print(f"C. s = {E['s']:.5f}; elasticidades analíticas {E['analitica']}; numéricas "
          f"{ {k: round(v, 6) for k, v in E['numerica'].items()} }")
    print(f"   máx. entre naves: m {max(e['m'] for e in E_naves.values()):.4f}, "
          f"M_I {min(e['M_I'] for e in E_naves.values()):.4f}; phi: ideal {el_phi['ideal']:.3f}, 1 atm {el_phi['atm']:.3f}")
    for c, nom, nat, rng, rs in R['filas']:
        print(f"   {c:4s} {nat:9s} dm* {[round(r['dm_star'], 3) for r in rs]}  P {[round(r['P_MW']) for r in rs]}  "
              f"H2/t {[round(r['H2_t'], 2) for r in rs]}")


if __name__ == '__main__':
    main()
