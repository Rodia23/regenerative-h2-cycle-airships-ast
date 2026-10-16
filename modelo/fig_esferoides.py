# [EN] Hulls as prolate spheroids and the added-mass coefficient (Section 3.1, Figure 2).
# =============================================================================
# fig_esferoides.py — Los cascos como esferoides prolatos y el k(l/D) de Lamb (paper 1, punto E, 2026-10-08)
#
# (a) Los cuatro cascos como esferoides de su largo y diámetro publicados, a la misma escala. En el LCA60T, en
#     discontinua, el esferoide de la misma esbeltez que encierra el V_casco del modelo (su volumen de helio
#     publicado, que está dentro de una envolvente mayor).
# (b) El coeficiente de masa añadida transversal k = beta0/(2 - beta0) contra la esbeltez (parametros.k_lamb),
#     con la esfera (0.5) y el cilindro infinito (1) como límites, y las cuatro naves.
# Poco texto en la figura (pedido del autor): solo el nombre corto de cada nave; el resto va al pie.
# Salidas: resultados/figuras/esferoides.png y resultados/datos/fig_esferoides.json (cifras del pie).
# =============================================================================
import math

import numpy as np

import parametros as p
import rutas

ORDEN = ['Hind', 'FW', 'ZR', 'ZNT']
CORTO = {'Hind': 'LZ 129', 'FW': 'LCA60T', 'ZR': 'ZR-3', 'ZNT': 'Zeppelin NT'}


def datos():
    naves = {}
    for nom, d in p.PLATAFORMAS.items():
        V_esf = math.pi / 6.0 * d['L'] * d['D'] ** 2
        naves[d['corto']] = {'nombre': nom, 'L_m': d['L'], 'D_m': d['D'], 'esbeltez': d['esbeltez'],
                             'k': p.k_lamb(d['esbeltez']), 'k_modelo': d['k'], 'V_casco_m3': d['V_casco'],
                             'V_esferoide_m3': V_esf, 'razon_V': V_esf / d['V_casco']}
    # LCA60T: el esferoide de la misma esbeltez con el V_casco del modelo
    fw = naves['FW']
    D_eq = (6.0 * fw['V_casco_m3'] / (math.pi * fw['esbeltez'])) ** (1.0 / 3.0)
    fw['D_equiv_m'], fw['L_equiv_m'] = D_eq, fw['esbeltez'] * D_eq
    otras = [naves[c]['razon_V'] for c in ORDEN if c != 'FW']
    return {'naves': naves, 'razon_V_otras': [min(otras), max(otras)]}


def figura(r):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401
    from matplotlib.patches import Ellipse

    col = {d['corto']: d['color'] for d in p.PLATAFORMAS.values()}
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 3.5), gridspec_kw={'width_ratios': [1.25, 1.0]})

    # (a) vista lateral, misma escala; cada esferoide empieza en x = 0 para leer el largo en el eje. Primero los ejes y
    # el ajuste de márgenes; los cascos se dibujan después, cuando ya se conoce el recuadro del panel
    X_MIN, X_MAX, MARGEN = -5.0, 300.0, 5.0
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_xlabel('Length (m)')
    ax.set_yticks([])
    ax.tick_params(axis='y', which='both', left=False, right=False)
    ax.grid(axis='y', visible=False)

    # (b) k(l/D) de Lamb
    f = np.linspace(1.0005, 10.0, 400)
    bx.plot(f, [p.k_lamb(x) for x in f], color='#555555', lw=1.6, label='Prolate spheroid', zorder=1)
    bx.axhline(1.0, color='#555555', lw=0.9, ls=(0, (4, 3)), zorder=0)
    bx.plot([1.0], [0.5], 'o', ms=6, mfc='white', mec='#555555', mew=1.2, zorder=3)
    for c in ORDEN:
        n = r['naves'][c]
        bx.plot([n['esbeltez']], [n['k']], 'o', ms=7, color=col[c], mec='white', mew=1.2, zorder=4, label=CORTO[c])
    bx.set_xlim(0.7, 10.3)
    bx.set_ylim(0.45, 1.05)
    bx.set_xlabel(r'Fineness ratio $\ell/D$')
    bx.set_ylabel(r'Added-mass coefficient $k$')
    bx.legend(loc='lower right', fontsize=8.5)

    fig.tight_layout()
    # (a) y (b) de la misma altura (autor, 2026-10-08): con aspecto 1 el alto de (a) lo fijan los cascos, así que la
    # separación entre ellos se elige para que el dibujo llene el recuadro de (a), que es tan alto como el de (b)
    caja = ax.get_position()
    ancho_in, alto_in = fig.get_size_inches()
    alto_datos = (X_MAX - X_MIN) * (caja.height * alto_in) / (caja.width * ancho_in)
    suma_D = sum(r['naves'][c]['D_m'] for c in ORDEN)
    sep = (alto_datos - suma_D - 2.0 * MARGEN) / (len(ORDEN) - 1)
    y = -MARGEN
    for c in ORDEN:
        n = r['naves'][c]
        yc = y - n['D_m'] / 2.0
        ax.add_patch(Ellipse((n['L_m'] / 2.0, yc), n['L_m'], n['D_m'], facecolor=col[c], alpha=0.15,
                             edgecolor='none'))
        ax.add_patch(Ellipse((n['L_m'] / 2.0, yc), n['L_m'], n['D_m'], facecolor='none', edgecolor=col[c], lw=1.6))
        if c == 'FW':
            ax.add_patch(Ellipse((n['L_m'] / 2.0, yc), n['L_equiv_m'], n['D_equiv_m'], facecolor='none',
                                 edgecolor=col[c], lw=1.1, ls=(0, (4, 3))))
        ax.text(n['L_m'] + 6.0, yc, CORTO[c], va='center', ha='left', fontsize=9, color='#222222')
        y -= n['D_m'] + sep
    ax.set_ylim(-alto_datos, 0.0)
    ax.set_aspect('equal')       # a escala; el recuadro ya tiene esta proporción, así que no se encoge
    # rótulos (a) y (b) a la misma altura, sobre el borde izquierdo de cada panel
    alto = max(ax.get_position().y1, bx.get_position().y1) + 0.02
    for eje, letra in ((ax, '(a)'), (bx, '(b)')):
        fig.text(eje.get_position().x0, alto, letra, fontsize=11, va='bottom', ha='left')
    fig.savefig(rutas.figura('esferoides'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    r = datos()
    rutas.guardar_datos('fig_esferoides', r)
    figura(r)
    fw = r['naves']['FW']
    print(f"esferoide / V_casco: {r['razon_V_otras'][0]:.2f}-{r['razon_V_otras'][1]:.2f} (otras), "
          f"{fw['razon_V']:.2f} (LCA60T); LCA60T con su V_casco: {fw['L_equiv_m']:.0f} x {fw['D_equiv_m']:.0f} m")


if __name__ == '__main__':
    main()
