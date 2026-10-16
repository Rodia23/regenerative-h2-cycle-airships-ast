# [EN] Plot style of the figures.
# =============================================================================
# plot_style.py — estilo visual "de journal" para todas las figuras del paper
#
# Uso: agregar, justo despues de "import matplotlib.pyplot as plt":
#
#     import plot_style
#
# No cambia colores ni datos de ningun script -- solo tipografia, grilla,
# grosor de linea y estilo general, via rcParams globales de matplotlib.
# No requiere paquetes externos (SciencePlots, etc.), solo matplotlib base.
# =============================================================================

import matplotlib.pyplot as plt

plt.rcParams.update({
    # Tipografia serif (look "de journal" en vez del sans-serif por defecto)
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Nimbus Roman', 'Liberation Serif', 'DejaVu Serif', 'serif'],
    'mathtext.fontset': 'stix',   # simbolos matematicos estilo Times, no el azul/sans por defecto

    'font.size': 11,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'axes.labelweight': 'bold',      # títulos de los ejes en negrita (autor, 2026-10-12)
    'pdf.fonttype': 42,              # PDF con las fuentes TrueType incrustadas (Elsevier: Type 1 o TrueType)
    'ps.fonttype': 42,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,

    # Lineas y ejes: mas deliberado, menos "default"
    'lines.linewidth': 1.6,
    'axes.linewidth': 0.8,
    'axes.edgecolor': '#333333',
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,

    # Grilla sutil (ayuda a leer valores sin dominar la figura)
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.5,
    'grid.linestyle': '--',
    'axes.axisbelow': True,

    # Leyenda: marco fino, fondo solido (no el gris/alpha por defecto)
    'legend.frameon': True,
    'legend.framealpha': 0.95,
    'legend.edgecolor': '#333333',
    'legend.fancybox': False,

    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})


# ── Figuras del paper (autor, 2026-10-12) ────────────────────────────────────────────────────────────────────────
# Una sola escala tipográfica para las Figs. 2 a 6, a ancho completo (7.4 in, ~188 mm; 190 mm en la revista). Elsevier
# pide letras de 7 pt impresos como mínimo (6 pt en subíndices); los rótulos dentro de la gráfica van a 7.5 como mínimo.
T_EJE, T_MARCA, T_PANEL, T_LEYENDA, T_ROTULO = 10, 9, 10, 8, 7.5


def normalizar_paper(fig):
    """Aplica la escala del paper: títulos de ejes, números de los ejes, títulos de panel, leyendas y un mínimo para los
    rótulos dentro de la gráfica (los más grandes quedan como están). Llamar antes del layout."""
    for ax in fig.axes:
        for lab in (ax.xaxis.label, ax.yaxis.label):
            lab.set_fontsize(T_EJE)
            lab.set_fontweight('bold')
        ax.tick_params(labelsize=T_MARCA)
        for tit in (ax.title, ax._left_title, ax._right_title):
            tit.set_fontsize(T_PANEL)
            tit.set_fontweight('normal')
        if ax.get_legend():
            for tx in ax.get_legend().get_texts():
                tx.set_fontsize(T_LEYENDA)
        for tx in ax.texts:
            if tx.get_fontsize() < T_ROTULO:
                tx.set_fontsize(T_ROTULO)
    for leg in fig.legends:
        for tx in leg.get_texts():
            tx.set_fontsize(T_LEYENDA)


def guardar_paper(fig, ruta_png, dpi=300):
    """Guarda la figura del paper o del suplementario en PNG (los documentos de trabajo). Con FIGURAS_PDF=1 en el
    entorno, también en PDF vectorial (lo pide exportar_figuras_envio.py para la versión de envío; Elsevier no lista
    el PNG). Desde 2026-10-13 también la usan las seis figuras del suplementario."""
    import os
    fig.savefig(ruta_png, dpi=dpi, bbox_inches='tight')
    if os.environ.get('FIGURAS_PDF') == '1':
        fig.savefig(ruta_png[:-4] + '.pdf', bbox_inches='tight')
