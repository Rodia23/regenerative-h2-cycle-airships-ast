# [EN] Hulls, altitude responses and design chart of the four airships (Figures 2-4).
# =============================================================================
# exp18_por_nave.py — EXP-18: las figuras y tablas del paper 1 para las cuatro naves (pedido del autor, 2026-10-11).
# EN EL PAPER (autor, 2026-10-11): exp18_cascos_3d.png reemplaza a la figura de los esferoides (sección 3) y
# exp18_alturas_por_nave.png reemplaza a la figura 3 de exp09 (sección 5.1; macro AlturasSube; el agua, al texto). Lo demás (la 4 x 2 con el
# agua, las tablas por nave) queda aparte. Está en correr_todo.py.
#
# (1) Los cascos en 3D (propuesta para el panel (a) de la figura de los esferoides): esferoides de su largo y
#     diámetro publicados (fig_esferoides.datos), a la misma escala, en proyección ortográfica; solo la malla (sin
#     polígonos, pedido del autor), plena adelante y tenue atrás, con su contorno. El LCA60T es el esferoide de la
#     misma esbeltez con su V_casco (el helio publicado), y en discontinua su envolvente publicada de 200 x 50 m (en la
#     figura actual es al revés: la envolvente llena y el esferoide en discontinua). El panel (b), k(l/D) de Lamb.
# (2) La figura 3 (exp09_articulo: altitud y agua retenida con la dinámica del condensador, 1 t, umbral de 1 m) para
#     las cuatro naves, 4 filas x 2 columnas, con los mismos casos: ideal, sin límite, 100 MW y 10 MW (ventear).
# (3) Lo que cambia con la nave en las tablas 6 (freno) y 7 (rutas):
#     - tabla 7: solo "Power, at once" y su tiempo. Todas las rutas operan en el tiempo natural de la nave,
#       t_nat = m / (L ṁ*) (sección 8), así que esa potencia es el calor por tonelada / t_nat; el resto es por
#       tonelada (y el tiempo de recuperación con 10 MW no depende de la nave);
#     - tabla 6: el escalón que absorbe la quema (m_e, de la potencia de tiempo real de la nave), la fuerza y la
#       potencia de las hélices; la grilla termina donde los haces de ACC pesan la carga útil de CADA nave, y solo
#       entran los escalones que caben en su carga útil. El tiempo de recuperación no cambia.
#     Funciones de freno.py (reparto, tiempo_real) y de helices_vectoriales.py (respuesta), sin copiarlas.
# Salidas: resultados/figuras/exp18_cascos_3d.png, exp18_agua_por_nave.png, exp18_alturas_por_nave.png (solo la
#          altura, 0-350 s, la misma ventana en las cuatro naves); resultados/datos/exp18_por_nave.json;
#          resultados/tablas/exp18_tablas_por_nave.txt; exp18_diseno_carga_tiempo.png (parte 4: diagrama de diseño,
#          dos paneles; en el paper, sección 5.2, con la curva de 1 atm como comparación en el panel b)
# =============================================================================
import io
import math
import sys

import numpy as np

import parametros as p
import rutas
import condensador
import drop_flujo_optimo_perdrop as sim
import fig_esferoides as fe
import freno
from helices_vectoriales import respuesta

ORDEN = fe.ORDEN                       # LZ 129, LCA60T, ZR-3, Zeppelin NT (como la figura actual)
CORTO = fe.CORTO
NOMBRE = {d['corto']: n for n, d in p.PLATAFORMAS.items()}
CASOS = (('ideal', None, 'Instantaneous condensation (ideal bound)', 'black', '-'),
         ('inf', None, 'Unlimited power', '#0072B2', '-'),
         ('100', 100e6, '100 MW', '#E69F00', '--'),
         ('10', 10e6, '10 MW', '#D55E00', ':'))           # los de exp09_articulo
ESCALONES_T = (1.0, 10.0, 30.0)                            # los de la tabla 6, si caben en la carga útil
T_VENTANA = 150.0                                          # [s] ventana de la figura del agua (la suelta, ~30 s)

# ── (1) cascos en 3D ─────────────────────────────────────────────────────────
ELEV, AZIM = 22.0, -58.0                   # vista de la figura 3D [grados]
N_MERIDIANOS, PASO_ANILLO_M = 28, 15.0     # malla (pedido del autor: solo la malla, sin polígonos); los anillos a la
                                           # misma distancia en todos los cascos, así la malla también da la escala


def _ojo():
    """Dirección hacia el observador de matplotlib para ELEV y AZIM (ejes a la misma escala)."""
    e, a = math.radians(ELEV), math.radians(AZIM)
    return np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])


def _malla(ax, L, D, x_c, z_c, color, n=160):
    """Malla del esferoide (meridianos y anillos), sin polígonos: los tramos que miran al observador, plenos; los de
    atrás, tenues. Visibilidad de cada punto por el signo de n·ojo, con n la normal del esferoide."""
    d = _ojo()
    a, b = L / 2.0, D / 2.0

    def linea(u, v):
        x, y, z = x_c - a * np.cos(v), b * np.sin(v) * np.cos(u), z_c + b * np.sin(v) * np.sin(u)
        ve = -np.cos(v) / a * d[0] + np.sin(v) * (np.cos(u) * d[1] + np.sin(u) * d[2]) / b >= 0.0
        for m, kw in ((~ve, dict(alpha=0.16, lw=0.45, zorder=1)), (ve, dict(alpha=0.9, lw=0.6, zorder=2))):
            ax.plot(np.where(m, x, np.nan), np.where(m, y, np.nan), np.where(m, z, np.nan), color=color, **kw)

    s = np.linspace(0.0, 1.0, n)
    for u in np.linspace(0.0, 2.0 * np.pi, N_MERIDIANOS, endpoint=False):
        linea(np.full(n, u), np.pi * s)
    for xa in np.arange(PASO_ANILLO_M, L, PASO_ANILLO_M):          # anillos cada PASO_ANILLO_M desde la cola
        linea(2.0 * np.pi * s, np.full(n, np.arccos(1.0 - xa / a)))


def _silueta(L, D, x_c, z_c, n=600):
    """Contorno aparente del esferoide (normal perpendicular a la vista): R cos(u - u0) = (b/a) d_x cot v."""
    d = _ojo()
    a, b = L / 2.0, D / 2.0
    R, u0 = math.hypot(d[1], d[2]), math.atan2(d[2], d[1])
    # el contorno existe donde |c| <= 1: entre v1 y pi - v1, con cot v1 = R a / (b |d_x|); denso cerca de las puntas
    v1 = math.atan2(b * abs(d[0]), R * a) * (1.0 + 1e-9)
    v = v1 + (math.pi - 2.0 * v1) * 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, n)))
    c = np.clip((b / a) * d[0] / np.tan(v) / R, -1.0, 1.0)
    vv = np.concatenate([v, v[::-1], v[:1]])
    u = np.concatenate([u0 + np.arccos(c), (u0 - np.arccos(c))[::-1], u0 + np.arccos(c[:1])])
    return x_c - a * np.cos(vv), b * np.sin(vv) * np.cos(u), z_c + b * np.sin(vv) * np.sin(u)


def figura_cascos(r):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401

    col = {d['corto']: d['color'] for d in p.PLATAFORMAS.values()}
    fig = plt.figure(figsize=(7.4, 3.3))
    ax = fig.add_axes([-0.03, 0.0, 0.68, 1.0], projection='3d')
    bx = fig.add_axes([0.685, 0.17, 0.30, 0.76])
    ax.computed_zorder = False              # los cascos no se tapan entre sí: se dibujan en el orden dado
    GAP = 16.0
    z, centros = 0.0, {}
    for c in ORDEN:
        z -= r['naves'][c]['D_m'] / 2.0
        centros[c] = z
        z -= r['naves'][c]['D_m'] / 2.0 + GAP
    for c in ORDEN:
        n, zc = r['naves'][c], centros[c]
        xc = n['L_m'] / 2.0                 # todos empiezan en x = 0 (se lee el largo en la regla)
        L, D = (n['L_equiv_m'], n['D_equiv_m']) if c == 'FW' else (n['L_m'], n['D_m'])
        _malla(ax, L, D, xc, zc, col[c])
        xs, ys, zs = _silueta(L, D, xc, zc)                         # el contorno, más grueso: la forma se lee limpia
        ax.plot(xs, ys, zs, color=col[c], lw=1.2, zorder=3)
        if c == 'FW':                       # la envolvente publicada, en discontinua
            xs, ys, zs = _silueta(n['L_m'], n['D_m'], xc, zc)
            ax.plot(xs, ys, zs, color=col[c], lw=1.1, ls=(0, (4, 3)), zorder=3)
        ax.text(n['L_m'] + 8.0, 0.0, zc, CORTO[c], fontsize=9, color='#222222', va='center', ha='left', zorder=4)
    # regla de largo bajo los cascos (sin cajas ni grillas 3D)
    z_r = z + GAP * 0.5
    kw = dict(color='#333333', lw=0.9, zorder=4)
    ax.plot([0, 250], [0, 0], [z_r, z_r], **kw)
    for xt in range(0, 251, 50):
        ax.plot([xt, xt], [0, 0], [z_r, z_r + 3.5], **kw)
        ax.text(xt, 0.0, z_r - 8.0, f'{xt}', fontsize=8.5, ha='center', va='top', color='#222222', zorder=4)
    ax.text(270.0, 0.0, z_r, 'm', fontsize=9, ha='left', va='center', color='#222222', zorder=4)
    ax.set_proj_type('ortho')
    x_lim, y_lim, z_lim = (-8.0, 290.0), (-28.0, 28.0), (z_r - 14.0, 2.0)   # margen: no recortar las puntas
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_zlim(*z_lim)
    ax.set_box_aspect((np.ptp(x_lim), np.ptp(y_lim), np.ptp(z_lim)), zoom=1.45)   # misma escala en los 3 ejes
    ax.view_init(elev=ELEV, azim=AZIM)
    ax.set_axis_off()

    # (b) k(l/D) de Lamb, como en fig_esferoides
    f = np.linspace(1.0005, 10.0, 400)
    # (autor, 2026-10-12) la curva punteada, sin la línea de k = 1 ni el punto de la esfera
    bx.plot(f, [p.k_lamb(x) for x in f], color='#555555', lw=1.4, ls=(0, (1, 1.6)), zorder=1)
    for c in ORDEN:
        n = r['naves'][c]
        bx.plot([n['esbeltez']], [n['k']], 'o', ms=7, color=col[c], mec='white', mew=1.2, zorder=4, label=CORTO[c])
    bx.set_xlim(0.7, 10.3)
    bx.set_ylim(0.45, 1.05)
    bx.set_xlabel(r'Fineness ratio $\ell/D$')
    bx.set_ylabel(r'Added-mass coefficient $k$')
    bx.legend(loc='lower right', fontsize=8.5)
    # sin (a) y (b): el texto dice izquierda y derecha (autor, 2026-10-12)
    plot_style.normalizar_paper(fig)
    plot_style.guardar_paper(fig, rutas.figura('exp18_cascos_3d'))
    plt.close(fig)


# ── (2) la figura 3 para las cuatro naves ───────────────────────────────────
def correr_agua():
    series, resumen = {}, {}
    for c in ORDEN:
        nom = NOMBRE[c]
        series[c], resumen[c] = {}, {}
        for clave, P, *_ in CASOS:
            cond = None if clave == 'ideal' else condensador.config(P, tau_s=p.tau_condensador, modo='ventear')
            sim.usar_nave(nom, 1000.0)
            r = sim.correr(n_entregas=1, toneladas_por_drop=1.0, graficar=False, verbose=False, t_post=600.0,
                           cond=cond, devolver_serie=True)
            s = r['serie']
            series[c][clave] = {'t': np.asarray(s['t']), 'dh': np.asarray(s['h']) - p.h0, 'agua': np.asarray(s['agua']),
                                't_rel': s['t_drops'][0], 't_anticipo': r['t_anticipo']}
            resumen[c][clave] = {'t_rel_s': s['t_drops'][0], 't_anticipo_s': r['t_anticipo'],
                                 'dh_arriba_m': r['dh_arriba_m'], 'dh_abajo_m': r['dh_abajo_m'],
                                 'agua_retenida_kg': r['m_agua_fin'], 'agua_venteada_kg': r['m_agua_venteada'],
                                 'cumple_h_max': bool(max(r['dh_arriba_m'], -r['dh_abajo_m']) <= r['h_max_m'])}
    sim.usar_nave(p.nombre_dirigible)
    return series, resumen


def figura_agua(series):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401

    # misma ventana de tiempo y misma escala por columna: las naves se comparan a simple vista
    fig, axes = plt.subplots(4, 2, figsize=(7.4, 9.0), sharex=True, sharey='col')
    for i, c in enumerate(ORDEN):
        s0 = series[c]['ideal']
        t_rel = s0['t_rel']
        t_arr = t_rel - s0['t_anticipo'] - p.n_tau_arranque * p.tau_condensador
        for clave, _, etiqueta, colr, ls in CASOS:
            s = series[c][clave]
            axes[i, 0].plot(s['t'], s['dh'], color=colr, ls=ls, lw=1.4, label=etiqueta)
            axes[i, 1].plot(s['t'], s['agua'], color=colr, ls=ls, lw=1.4)
        for ax in axes[i]:
            ax.axvline(t_arr, color='gray', ls=':', lw=0.9)
            ax.axvline(t_rel, color='gray', ls='--', lw=0.9)
            ax.grid(True, alpha=0.2)
        axes[i, 0].axhspan(-p.h_max_precision, p.h_max_precision, color='gray', alpha=0.12, lw=0)
        axes[i, 0].set_ylabel('Altitude change (m)')
        axes[i, 1].set_ylabel('Water retained (kg)')
        axes[i, 0].text(0.98, 0.92, CORTO[c], transform=axes[i, 0].transAxes, fontsize=10, va='top', ha='right',
                        color='#222222', bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.5))
    axes[0, 0].set_xlim(0.0, T_VENTANA)
    for ax in axes[-1]:
        ax.set_xlabel('Time (s)')
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc='upper center', ncol=2, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(rutas.figura('exp18_agua_por_nave'), dpi=300, bbox_inches='tight')
    plt.close(fig)


T_ASIENTO = 350.0      # [s] ventana de la figura de alturas: el caso más lento (LZ 129, 10 MW) queda bajo 0.1 m a ~300 s


def figura_alturas(series):
    """Solo la altura (autor, 2026-10-11), una nave por fila, con una ventana larga y la misma para todas, para ver
    cuándo se asienta cada caso. (2026-10-12) Los colores de CASOS (el autor prefirió estos a las curvas en negro), la
    subida máxima con 10 MW marcada, espacio arriba para su número y la leyenda abajo."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401

    fig, axes = plt.subplots(len(ORDEN), 1, figsize=(7.4, 8.2), sharex=True, sharey=True)
    for ax, c in zip(axes, ORDEN):
        s0 = series[c]['ideal']
        t_arr = s0['t_rel'] - s0['t_anticipo'] - p.n_tau_arranque * p.tau_condensador
        for clave, _, etiqueta, colr, ls in CASOS:
            s = series[c][clave]
            ax.plot(s['t'], s['dh'], color=colr, ls=ls, lw=1.5, label=etiqueta, zorder=3)
        s10 = series[c]['10']
        i = int(np.argmax(s10['dh']))
        ax.plot(s10['t'][i], s10['dh'][i], 'o', ms=5.5, color=CASOS[-1][3], mec='white', mew=1.1, zorder=5)
        ax.annotate(f"{s10['dh'][i]:.1f} m", xy=(s10['t'][i], s10['dh'][i]), xytext=(0, 7), textcoords='offset points',
                    fontsize=8, fontweight='bold', va='bottom', ha='center', color='#111111', zorder=6,
                    bbox=dict(facecolor='white', edgecolor='none', pad=1.0, alpha=0.9))
        # (autor, 2026-10-12) sin la línea del arranque de la unidad (t_arr; el texto lo dice): solo la suelta
        ax.axvline(s0['t_rel'], color='gray', ls='--', lw=0.9)
        ax.axhspan(-p.h_max_precision, p.h_max_precision, color='gray', alpha=0.12, lw=0, zorder=0)
        ax.set_ylabel('Altitude change (m)')
        ax.grid(True, alpha=0.2)
        ax.set_title(CORTO[c], fontsize=9.5, loc='left', pad=3, color='#111111')      # sin negrita (autor)
    axes[0].set_xlim(0.0, T_ASIENTO)
    # aire arriba para el número de la subida máxima (autor, 2026-10-12)
    pico = max(float(np.max(series[c][k]['dh'])) for c in ORDEN for k in series[c])
    axes[0].set_ylim(-1.25 * p.h_max_precision, 1.3 * pico)
    axes[-1].set_xlabel('Time (s)')
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc='lower center', ncol=2, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 0.0))
    plot_style.normalizar_paper(fig)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    plot_style.guardar_paper(fig, rutas.figura('exp18_alturas_por_nave'))
    plt.close(fig)


# ── (3) tablas 6 y 7 por nave ────────────────────────────────────────────────
def tablas():
    cand = {c['clave']: c for c in rutas.leer_datos('req_masa_rechazo')['candidatos']}
    kwkg = cand['ACC_haces']['kWkg'][0]
    e12 = rutas.leer_datos('exp12_helio')
    R = e12['rutas']
    if not (R[0]['ruta'].startswith('H2 quemado') and '1 atm' in R[1]['ruta'] and 'COSH' in R[3]['ruta']):
        sys.exit('rutas de exp12: el orden cambió')
    # las cuatro columnas de la tabla 7 con calor (las venteadas no tienen)
    calor = {'H2 burned (phi 0.99)': R[0]['calor_GJ_por_t'], 'H2 burned, 1 atm': R[1]['calor_GJ_por_t'],
             'H2 compressed (35 MPa)': e12['h2_comprimido']['calor_GJ_por_t'],
             'He compressed (COSH)': R[3]['calor_GJ_por_t']}
    naves = {}
    for c in ORDEN:
        nom = NOMBRE[c]
        d = p.PLATAFORMAS[nom]
        rt = freno.tiempo_real(nom, 1000.0)
        P_lim = d['carga_max'] / 1000.0 * kwkg * 1e6                  # [W] haces de ACC = carga útil de la nave
        grilla = sorted({x for x in p.freno_P_grilla if x <= P_lim * (1.0 + 1e-9)} | {P_lim})
        filas = [freno.reparto(nom, m * 1000.0, P, cand) for m in ESCALONES_T if m * 1000.0 <= d['carga_max'] + 1e-6
                 for P in grilla]
        resp = respuesta(nom, 1000.0)
        naves[c] = {'nombre': nom, 'carga_t': d['carga_max'] / 1e3, 'P_lim_MW': P_lim / 1e6,
                    't_nat_1t_s': rt['t_nat_s'], 'dm_star_1t_kgps': rt['dm_star'], 'P_rt_1t_MW': rt['P_rt_W'] / 1e6,
                    'tau_techo_1t_s': resp['tau_techo_s'],
                    # calor / t_nat, salvo la columna de phi 0.99, que en la tabla 7 es ṁ* q (la tasa de la ec. 7 se
                    # dimensiona con phi = 1; tablas_articulo, 'golpe'): 446 y no 448 MW en el LCA60T
                    'potencia_golpe_MW': dict({k: v * 1e9 / rt['t_nat_s'] / 1e6 for k, v in calor.items()},
                                              **{'H2 burned (phi 0.99)': rt['P_rt_W'] / 1e6}),
                    'filas': [{'m_t': f['m_kg'] / 1e3, 'P_MW': f['P_W'] / 1e6, 'm_e_kg': f['m_e_kg'],
                               't_rec_s': f['t_rec_s'], 'F_kN': f['F_N'] / 1e3,
                               'P_hel_MW': [min(f['P_hel_W']) / 1e6, max(f['P_hel_W']) / 1e6]} for f in filas]}
    return {'kWkg_haces': kwkg, 'calor_GJ_por_t': calor, 'naves': naves}


def _t(t_s):
    return f'{t_s / 60:.1f} min' if t_s < 2 * 3600 else f'{t_s / 3600:.1f} h'


def escribir_txt(tb, agua):
    L = ['EXP-18 — LAS TABLAS 6 Y 7 DEL PAPER 1 PARA LAS CUATRO NAVES (generado por modelo/exp18_por_nave.py)',
         'Aparte: no está en el paper. 1 t liberada, umbral de 1 m, D_lim 0.10 m, como el caso base.', '',
         'TABLA 7 (rutas por tonelada): solo cambia "Power, at once" y su tiempo. Todas las rutas operan en el tiempo',
         'natural de la nave (t_nat = m / (L ṁ*), sección 8): potencia = calor por tonelada / t_nat. Las demás filas',
         '(L, gas procesado, calor, potencia neta con 10 MW, tanques, energía, gas perdido) son por tonelada: iguales.', '']
    rutas_ = list(tb['calor_GJ_por_t'])
    L.append(f"{'':<14}{'t_nat (s)':>10}" + ''.join(f'{k[:24]:>26}' for k in rutas_))
    for c in ORDEN:
        n = tb['naves'][c]
        L.append(f"{CORTO[c]:<14}{n['t_nat_1t_s']:>10.1f}" +
                 ''.join(f"{n['potencia_golpe_MW'][k]:>23.0f} MW" for k in rutas_))
    L += ['', '(calor por tonelada, GJ: ' + ', '.join(f'{k} {v:.2f}' for k, v in tb['calor_GJ_por_t'].items()) + ')',
          '', 'TABLA 6 (freno): cambian el escalón que absorbe la quema (m_e), la fuerza y la potencia de las hélices (disco',
          f"{p.freno_A_discos[0]:.0f}-{p.freno_A_discos[1]:.0f} m2, FM {p.freno_FM}). El tiempo de recuperación NO cambia. "
          'La grilla termina donde los',
          f"haces de ACC ({tb['kWkg_haces']:.3f} kW/kg) pesan la carga útil de cada nave; solo escalones que caben en ella.", '']
    for c in ORDEN:
        n = tb['naves'][c]
        L.append(f"{n['nombre']} (carga útil {n['carga_t']:.1f} t; haces = carga útil con {n['P_lim_MW']:.2f} MW; "
                 f"tiempo real de 1 t: {n['P_rt_1t_MW']:.0f} MW en {n['t_nat_1t_s']:.1f} s; "
                 f"techo de respuesta del freno, 1 t: {n['tau_techo_1t_s']:.1f} s)")
        L.append(f"  {'escalón':>8}{'P_cap (MW)':>12}{'m_e (kg)':>10}{'recuperación':>14}{'fuerza (kN)':>13}{'hélices (MW)':>15}")
        for f in n['filas']:
            L.append(f"  {f['m_t']:>6.0f} t{f['P_MW']:>12.2f}{f['m_e_kg']:>10.0f}{_t(f['t_rec_s']):>14}"
                     f"{f['F_kN']:>13.1f}{f['P_hel_MW'][0]:>9.2f}-{f['P_hel_MW'][1]:.2f}")
        L.append('')
    L += ['FIGURA 3 POR NAVE (1 t, tau = %g s): subida y bajada máximas, agua retenida al final' % p.tau_condensador]
    for c in ORDEN:
        L.append(f'  {CORTO[c]}:')
        for clave, _, etiqueta, *_ in CASOS:
            x = agua[c][clave]
            L.append(f"    {etiqueta:<42} sube {x['dh_arriba_m']:6.2f} m, baja {x['dh_abajo_m']:6.2f} m, "
                     f"retiene {x['agua_retenida_kg']:6.0f} kg, ventea {x['agua_venteada_kg']:6.0f} kg, "
                     f"cumple 1 m: {'sí' if x['cumple_h_max'] else 'no'}")
    with open(rutas.tabla('exp18_tablas_por_nave').replace('.tex', '.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    return '\n'.join(L)


# ── (4) diagrama de diseño: la carga más grande según el tiempo de recuperación (autor, 2026-10-11) ──────────
# Con el freno se eligen la carga m y el tiempo t; el rechazo instalado es P = m Q / t (ecuación de recuperación,
# Q = calor por tonelada); lo que escala con P pesa a t por MW (haces, motor, compresor del escape: los mismos de
# fig_presupuesto_masa), las cajas fijas b t por tonelada liberada, y la carga m. Con la carga útil C:
#     m_max(t) = C / (1 + b + a Q / t),
# que tiende a C / (1 + b) con t largo (manda lo que va por tonelada; desde 2026-10-13 con el compresor de H2 en
# unidades de catálogo, fpm.suelta_unica_max) y a P = C / a con t corto (manda el rechazo).
# Clase más liviana y compresor del escape (14 atm), como la figura del presupuesto. Las fijas, lineales en m.
# EN EL PAPER (autor, 2026-10-11): dos paneles, (a) m_max por nave en toneladas y (b) m_max / C, donde las cuatro
# curvas se vuelven una y cada nave entra solo por la fracción de su carga útil que es una tonelada. La base queda a
# 14 atm; 1 atm, sin el compresor del escape, va como comparación en (b) (el plan de pasar la base a 1 atm se descartó).
T_DISENO_S = (60.0, 2.0 * 86400.0)          # [s] de 1 min a 48 h
P_ISO_MW = (0.1, 1.0, 10.0, 100.0)          # líneas de rechazo constante
ORDEN_LEYENDA = ('ZNT', 'ZR', 'FW', 'Hind')  # de la nave más chica a la más grande (autor, 2026-10-12)


def datos_diseno():
    # (2026-10-13) el compresor de H2 en unidades de catálogo (fpm.compresor_t): m (1 + b' + a Q / t) + compresor(m)
    # <= C. Con 1 t es una unidad, así que b (las fijas de 1 t) sigue valiendo en el punto de 1 t; el techo ya no es
    # C / (1 + b) y la fracción de la carga útil deja de ser la misma en toda nave (la unidad pesa más en las chicas)
    import fig_presupuesto_masa as fpm
    d = fpm.datos_apilado()
    a, b, Q = d['t_por_MW'], d['fijas_t'], d['calor_t_J'] / 1e6           # [t/MW], [t/t], [MJ/t]
    bp, cu = d['fijas_sin_comp_t'], d['compresor']
    # la misma carta a 1 atm, sin el compresor del escape (comparación, autor 2026-10-11): se retiene phi_max a la
    # salida del condensador; los haces pesan lo mismo por MW, el motor va a la tasa de 1 MW con q(0, phi_max) y Q es
    # el calor de la tonelada con L(phi_max), el mismo de la ruta de 1 atm de la tabla 7 (lo exige la verificación)
    # (2026-10-12) la cuenta vive en fig_presupuesto_masa.datos_apilado, que también da la barra HDM de 1 atm
    phi1 = d['una_atm']['phi']
    a1 = d['una_atm']['t_por_MW']                                                     # [t/MW]
    Q1 = d['una_atm']['calor_t_J'] / 1e6                                              # [MJ/t]
    t = np.geomspace(*T_DISENO_S, 300)
    naves = {}
    for nv in d['naves']:
        C = nv['carga_t']
        m = np.array([fpm.suelta_unica_max(C, bp, a * Q / ti, cu) for ti in t])
        m1 = np.array([fpm.suelta_unica_max(C, bp, a1 * Q1 / ti, cu) for ti in t])
        naves[nv['corto']] = {'carga_t': C, 'm_max_t': m.tolist(), 'm_max_una_atm_t': m1.tolist(),
                              'P_MW': (m * Q / t).tolist(), 'limite_t': fpm.suelta_unica_max(C, bp, 0.0, cu),
                              'P_corto_MW': C / a,
                              't_1t_s': a * Q / (C - 1.0 - b) if C > 1.0 + b else None,
                              't_1t_s_una_atm': a1 * Q1 / (C - 1.0 - b) if C > 1.0 + b else None,
                              'tabla': {f'{s:g}': fpm.suelta_unica_max(C, bp, a * Q / s, cu)
                                        for s in (600.0, 3600.0, 10800.0, 86400.0)}}
    # con la misma carga, el tiempo escala con a Q (t = a Q m / (C - compresor(m) - (1 + b') m)): la razón es la misma
    # en toda nave
    return {'a_t_por_MW': a, 'b_t_por_t': b, 'b_sin_comp_t_por_t': bp, 'compresor': cu, 'Q_MJ_por_t': Q,
            'aQ_s': a * Q, 't_s': t.tolist(), 'naves': naves,
            'una_atm': {'phi': phi1, 'a_t_por_MW': a1, 'Q_MJ_por_t': Q1, 'aQ_s': a1 * Q1,
                        'razon_t': a1 * Q1 / (a * Q)}}


def figura_diseno(dd, vertical=False, x_lineal=False, salida=None):
    """Carta de diseño (Fig. 4): la carga más grande que cabe por nave. (2026-10-14, autor) un solo panel: el de la
    derecha (dividido entre la carga útil) salió del paper; la curva normalizada de la suelta única está en la figura
    S6 y el caso de 1 atm, en la figura 5. x_lineal: variante de prueba; vertical: ya no se usa (era la disposición de
    los dos paneles); salida: ruta del png."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401
    from matplotlib.ticker import FixedLocator, NullFormatter
    from matplotlib.lines import Line2D

    col = {d['corto']: d['color'] for d in p.PLATAFORMAS.values()}
    t = np.array(dd['t_s'])
    a, b, Q = dd['a_t_por_MW'], dd['b_t_por_t'], dd['Q_MJ_por_t']
    u = dd['una_atm']
    esc = 3600.0 if x_lineal else 60.0                     # eje x en horas (lineal) o en minutos (log)
    x = t / esc
    fig, ax = plt.subplots(figsize=(4.8, 3.2))      # (2026-10-14) un panel; antes 1 x 2 de 7.4 x 3.9 in
    if not x_lineal:                                 # la escala antes de los rótulos fijos (si no, vuelven las potencias)
        ax.set_xscale('log')
    ax.set_yscale('log')
    x_max = T_DISENO_S[1] / esc
    x_tec = (0.93 if x_lineal else 0.62) * x_max     # rótulos de la derecha, un poco adentro del borde (autor)

    # izquierda: por nave, en toneladas; rechazo constante: m = P t / Q
    lineas_P = {}
    for P in P_ISO_MW:
        lineas_P[P] = ax.plot(x, P * t / Q, color='#c4c2bc', lw=0.7, ls=(0, (3, 3)), zorder=1)[0]
    for c in ORDEN:
        n = dd['naves'][c]
        ax.plot(x, n['m_max_t'], color=col[c], lw=1.6, zorder=3)
        ax.annotate(f"{CORTO[c]}, {fe_fmt(n['limite_t'])} t", xy=(x_tec, n['m_max_t'][-1]), xytext=(0, 3),
                    textcoords='offset points', fontsize=7, color='#222222', va='bottom', ha='right')
        if n['t_1t_s']:
            ax.plot(n['t_1t_s'] / esc, 1.0, 'o', ms=5, color=col[c], mec='white', mew=0.9, zorder=4)
    ax.set_ylim(0.01, 60.0)
    ax.yaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10]))
    ax.set_yticklabels(['0.01', '0.1', '1', '10'])
    ax.set_ylabel('Largest release that fits (t)')
    ax.legend(handles=[Line2D([], [], ls='none', marker='o', ms=5, color='#555555', mec='white', label='1-t release'),
                       Line2D([], [], color='#c4c2bc', lw=0.7, ls=(0, (3, 3)), label='constant heat rejection')],
              loc='lower right', fontsize=6.5, frameon=False)

    # (2026-10-14) el panel derecho (dividido entre la carga útil) salió del paper


    for e in (ax,):
        if x_lineal:
            e.set_xlim(0.0, x_max)
            e.xaxis.set_major_locator(FixedLocator(list(range(0, int(x_max) + 1, 6))))
            e.set_xlabel('Recovery time chosen (h)', fontsize=9)
        else:
            e.set_xlim(T_DISENO_S[0] / esc, x_max)
            e.xaxis.set_major_locator(FixedLocator([1, 10, 60, 600, 2880]))
            e.set_xticklabels(['1 min', '10 min', '1 h', '10 h', '48 h'])
            e.xaxis.set_minor_formatter(NullFormatter())
            e.set_xlabel('Recovery time chosen', fontsize=9)
        e.yaxis.set_minor_formatter(NullFormatter())
        e.yaxis.label.set_size(9)
        e.tick_params(labelsize=8)
        e.grid(True, which='major', alpha=0.25)
    # (2026-10-14) sin la leyenda de abajo: cada curva lleva el nombre de su nave
    plot_style.normalizar_paper(fig)
    fig.tight_layout()
    # rótulos de potencia sobre la diagonal del panel (de arriba a la izquierda a abajo a la derecha), girados con su
    # línea: se ubican después del layout, en coordenadas de pantalla
    fig.canvas.draw()
    a2d = ax.transAxes.inverted()
    for P, ln in lineas_P.items():
        xs, ys = ln.get_xdata(), ln.get_ydata()
        pts = ax.transData.transform(np.column_stack([xs, ys]))
        uv = a2d.transform(pts)
        f = uv[:, 0] + uv[:, 1] - 1.0
        k = np.where((np.sign(f[:-1]) != np.sign(f[1:])) & (uv[:-1, 0] > 0.02) & (uv[:-1, 0] < 0.98)
                     & (uv[:-1, 1] > 0.02) & (uv[:-1, 1] < 0.98))[0]
        if not len(k):
            continue
        i = int(k[0])
        d = pts[i + 1] - pts[i]
        ang = np.degrees(np.arctan2(d[1], d[0]))
        nrm = np.array([-d[1], d[0]]) / np.hypot(*d)                 # normal a la línea, hacia arriba-izquierda
        # el lado de la línea con más espacio libre respecto de las curvas de las naves
        libre = []
        for sgn in (1.0, -1.0):
            centro = pts[i] + sgn * 10.0 * nrm
            dmin = min(np.min(np.hypot(*(ax.transData.transform(np.column_stack([x, dd['naves'][c]['m_max_t']]))
                                         - centro).T)) for c in ORDEN)
            libre.append(dmin)
        sgn = 1.0 if libre[0] >= libre[1] else -1.0
        off = sgn * 3.0 * nrm * 72.0 / fig.dpi
        ax.annotate(f'{P:g} MW', xy=(xs[i], ys[i]), xytext=tuple(off), textcoords='offset points', fontsize=7.5,
                    color='#8a8882', ha='center', va='bottom' if sgn > 0 else 'top', rotation=ang,
                    rotation_mode='anchor')
    if salida:
        fig.savefig(salida, dpi=300, bbox_inches='tight')
    else:
        plot_style.guardar_paper(fig, rutas.figura('exp18_diseno_carga_tiempo'))
    plt.close(fig)


def fe_fmt(x):
    return f'{x:.0f}' if x >= 10 else (f'{x:.1f}' if x >= 1 else f'{x:.2f}')


def main():
    if (sys.stdout.encoding or '').lower() != 'utf-8':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    r = fe.datos()
    figura_cascos(r)
    series, agua = correr_agua()
    figura_agua(series)
    figura_alturas(series)
    tb = tablas()
    dd = datos_diseno()
    figura_diseno(dd)
    rutas.guardar_datos('exp18_por_nave', {'tablas': tb, 'agua': agua,
                                           'diseno': {k: v for k, v in dd.items() if k != 't_s'}})
    print(escribir_txt(tb, agua))
    L = ['', 'DIAGRAMA DE DISEÑO (con freno): la carga más grande que cabe según el tiempo de recuperación elegido',
         "m_max(t): m (1 + b' + a Q / t) + compresor(m) = C; P = m Q / t. Clase más liviana, compresor del escape",
         "(14 atm); tanques y PEM lineales en la carga, el compresor de H2 en unidades de catálogo (al menos una).",
         f"a = {dd['a_t_por_MW']:.3f} t por MW, b' = {dd['b_sin_comp_t_por_t']:.3f} t por tonelada liberada, una unidad "
         f"de {dd['compresor']['kg']:.0f} kg por {dd['compresor']['cap_kgph'] / dd['compresor']['kgph_por_t']:.2f} t "
         f"(b = {dd['b_t_por_t']:.2f} con 1 t), Q = "
         f"{dd['Q_MJ_por_t']:.0f} MJ por tonelada, a Q = {dd['aQ_s']:.0f} s"]
    for c in ORDEN:
        n = dd['naves'][c]
        L.append(f"  {CORTO[c]:<12} carga útil {n['carga_t']:.2f} t; límite (t largo) {n['limite_t']:.2f} t; rechazo con "
                 f"t corto {n['P_corto_MW']:.1f} MW; 1 t en " + ('— (no cabe)' if n['t_1t_s'] is None else
                                                                f"{n['t_1t_s'] / 60:.1f} min") + '; '
                 + ', '.join((f"{float(s) / 60:g} min" if float(s) < 3600 else f"{float(s) / 3600:g} h") + f": {m:.2f} t"
                            for s, m in n['tabla'].items()))
    u = dd['una_atm']
    L.append(f"A 1 atm, sin compresor del escape (comparación, panel b): phi = {u['phi']:.3f}, a = {u['a_t_por_MW']:.3f} "
             f"t por MW, Q = {u['Q_MJ_por_t']:.0f} MJ por tonelada; con la misma carga, el tiempo x {u['razon_t']:.3f} "
             'en toda nave; 1 t en ' + ', '.join(f"{CORTO[c]} {dd['naves'][c]['t_1t_s_una_atm'] / 60:.1f} min"
                                                 for c in ORDEN if dd['naves'][c]['t_1t_s_una_atm']) + '; el techo no cambia')
    with open(rutas.tabla('exp18_tablas_por_nave').replace('.tex', '.txt'), 'a', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
