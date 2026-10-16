# -*- coding: utf-8 -*-
# [EN] Reactive control versus clean burn (Section 4, Figure S4, Table S6).
# =============================================================================
# reactivo_vs_limpia.py — por qué se opera con quema previa (revisión 2026-10-01, sec. 4 y tabla 3 del paper 1)
#
# Sin propulsión, la única forma de ganar tiempo es empezar a quemar antes de soltar. Para cada nave (1 t,
# D_lim operacional) se compara la quema limpia del artículo con el control REACTIVO puro: sin quema previa, el
# controlador del apéndice B responde solo al soltar, saturado a la tasa instalada, con el agua retenida al
# instante (el mismo L que la quema limpia; t3_barrido_reactivo.py acreditaba solo la pérdida de flotación).
# Se busca, por bisección sobre la tasa instalada, la menor que mantiene la subida bajo cada umbral:
#   - la excursión de la quema limpia, (1+s) D_lim (la misma precisión);
#   - h_max de precisión (1 m) y de seguridad (5 m).
# Cota analítica: una quema reactiva IDEAL, a tasa plena desde el instante de soltar, frena la subida que deja la
# carga con la fuerza g (m - L ṁ t) sobre M_I - m: sube h = (2/3) g m^3 / ((M_I - m) L^2 ṁ^2), así que
#   ṁ_r(h) = sqrt(2 g m^3 / (3 (M_I - m) L^2 h)),   ṁ_r / ṁ* = sqrt(4 (1+s)^3 D_lim / (s^2 h)).
# El controlador real reacciona con retardo y necesita más.
# Barrido (autor, 2026-10-12): la subida reactiva simulada en una grilla de tasas instaladas y la excursión de la
# quema limpia con la forma cerrada de la ec. (5), (1+s) D_lim = g m^3 / (6 M_I L^2 (1+s)^2 dm^2), en los mismos ejes
# (figura del suplementario S7; reemplaza a barrido_reactivo.png, del modelo viejo sin el agua retenida, y a
# dlim_flujo.png). Con el MISMO desplazamiento h, la quema limpia necesita dm_c(h) y el reactivo su tasa de umbral:
# la razón reactivo / limpia a 1 m y a 5 m (razon_mismo).
# Salidas: resultados/datos/reactivo_vs_limpia.json (lo leen tablas_articulo.py, exportar_numeros.py y
#          verificacion.py, grupo Z); resultados/figuras/limpia_vs_reactivo.png
# =============================================================================

import io
import sys
import math
import parametros as p
import minimax as mm
import rutas
import drop_flujo_optimo_perdrop as sim

TOL = 1.0005          # bisección: hasta que hi/lo < TOL (0.05 % en la tasa)
TASAS_BARRIDO = [0.3 * 1000.0 ** (i / 35.0) for i in range(36)]      # 0.3 a 300 kg/s, log [kg/s]
ORDEN_LEYENDA = ('ZNT', 'ZR', 'FW', 'Hind')                         # de la nave más chica a la más grande


def subida(dm):
    """Subida máxima sobre h0 con control reactivo puro saturado a dm [kg/s]."""
    return sim.correr(graficar=False, verbose=False, dm_opt=dm, precond=False)['dh_arriba_m']


def tasa_minima(objetivo, lo, hi):
    """Menor tasa instalada (dentro de TOL) cuya subida no pasa de objetivo [m]; la subida baja con la tasa."""
    assert subida(lo) > objetivo >= subida(hi), (lo, hi, objetivo)
    while hi / lo >= TOL:
        mid = math.sqrt(lo * hi)
        if subida(mid) > objetivo:
            lo = mid
        else:
            hi = mid
    return hi


def excursion_limpia(M_I, m, L, s, dm):
    """Excursión (1+s) D_lim de la quema limpia con la tasa dm (ec. (5) despejada) [m]."""
    return p.g * m**3 / (6.0 * M_I * L**2 * (1.0 + s)**2 * dm**2)


def tasa_limpia(M_I, m, L, s, h):
    """Tasa de la quema limpia cuya excursión es h (la inversa de excursion_limpia) [kg/s]."""
    return math.sqrt(p.g * m**3 / (6.0 * M_I * L**2 * (1.0 + s)**2 * h))


def figura(filas, tasas):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401
    from matplotlib.ticker import FixedLocator, NullFormatter
    from matplotlib.lines import Line2D
    import fig_esferoides as fe

    por_corto = {f['corto']: (nom, f) for nom, f in filas.items()}
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.set_xscale('log')
    ax.set_yscale('log')
    for h, lab, ls in ((p.h_max_precision, f'precision {p.h_max_precision:g} m', ':'),
                       (p.h_max_seguridad, f'safety {p.h_max_seguridad:g} m', (0, (5, 3)))):
        ax.axhline(h, color='#555555', lw=0.9, ls=ls, zorder=1)
        ax.annotate(lab, xy=(0.93 * tasas[-1], h), xytext=(0, 2.5), textcoords='offset points', fontsize=7.5,
                    color='#555555', va='bottom', ha='right')              # a la derecha, donde no pasan las curvas
    for c in ORDEN_LEYENDA:
        nom, f = por_corto[c]
        col = p.PLATAFORMAS[nom]['color']
        ax.plot(tasas, f['barrido']['excursion_limpia_m'], color=col, lw=1.8, zorder=3)
        ax.plot(tasas, f['barrido']['subida_reactivo_m'], color=col, lw=1.4, ls=(0, (4, 2)), zorder=3)
        ax.plot(f['dm_limpia'], f['excursion_limpia_m'], 'o', ms=6, color=col, mec='white', mew=1.0, zorder=5)
        for k in ('precision', 'seguridad'):
            ax.plot(f['reactivo_sim'][k], f['objetivos_m'][k], 'o', ms=5.5, mfc='white', mec=col, mew=1.4, zorder=5)
    ax.set_xlim(tasas[0], tasas[-1])
    ax.set_ylim(0.01, 100)
    ax.xaxis.set_major_locator(FixedLocator([0.3, 1, 3, 10, 30, 100, 300]))
    ax.set_xticklabels(['0.3', '1', '3', '10', '30', '100', '300'])
    ax.yaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10, 100]))
    ax.set_yticklabels(['0.01', '0.1', '1', '10', '100'])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel('Installed burn rate (kg/s)', fontsize=9)
    ax.set_ylabel('Maximum vertical excursion (m)', fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(True, which='major', alpha=0.25)
    ax.legend(handles=[Line2D([], [], color='#333333', lw=1.8, label='clean burn: sink before the release, no rise'),
                       Line2D([], [], color='#333333', lw=1.4, ls=(0, (4, 2)),
                              label='reactive control: rise after the release'),
                       Line2D([], [], ls='none', marker='o', ms=6, color='#333333', mec='white',
                              label=f'clean burn at $D_{{lim}}$ = {p.D_lim_oper:g} m'),
                       Line2D([], [], ls='none', marker='o', ms=5.5, mfc='white', mec='#333333', mew=1.4,
                              label='reactive control at the thresholds')],
              loc='lower left', fontsize=7, frameon=False)
    fig.legend(handles=[Line2D([], [], color=p.PLATAFORMAS[por_corto[c][0]]['color'], lw=1.8, label=fe.CORTO[c])
                        for c in ORDEN_LEYENDA], loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=8,
               frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    plot_style.guardar_paper(fig, rutas.figura('limpia_vs_reactivo'))
    plt.close(fig)


def ideal(M_I, m, L, h):
    """Tasa de la quema reactiva ideal (tasa plena desde el instante de soltar) que sube h [m]."""
    return math.sqrt(2.0 * p.g * m**3 / (3.0 * (M_I - m) * L**2 * h))


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    m = float(p.m_carga)
    filas = {}
    for nom, d in p.PLATAFORMAS.items():
        sim.usar_nave(nom, m)
        r0 = sim.correr(graficar=False, verbose=False)            # la quema limpia del artículo
        dm_s = r0['dm_opt']
        M_I = r0['M_inercial']
        s = math.sqrt(1.0 - m / M_I)
        exc = abs(r0['dh_abajo_m'])                               # (1+s) D_lim, simulado
        objetivos = {'igual': exc, 'precision': p.h_max_precision, 'seguridad': p.h_max_seguridad}
        sim_r = {k: tasa_minima(h, 0.2 * dm_s, 400.0) for k, h in objetivos.items()}
        ide_r = {k: ideal(M_I, m, mm.L_eff, h) for k, h in objetivos.items()}
        filas[nom] = {'corto': d['corto'], 'dm_limpia': dm_s, 'subida_limpia_m': r0['dh_arriba_m'],
                      'excursion_limpia_m': exc, 's': s, 'M_I': M_I, 'objetivos_m': objetivos,
                      'reactivo_sim': sim_r, 'reactivo_ideal': ide_r,
                      'razon_sim': {k: v / dm_s for k, v in sim_r.items()},
                      'razon_ideal': {k: v / dm_s for k, v in ide_r.items()},
                      'razon_ideal_formula': {k: math.sqrt(4.0 * (1.0 + s)**3 * p.D_lim_oper / (s**2 * h))
                                              for k, h in objetivos.items()}}
        f = filas[nom]
        # barrido y mismo desplazamiento (2026-10-12): la nave sigue cargada en sim
        f['barrido'] = {'subida_reactivo_m': [subida(x) for x in TASAS_BARRIDO],
                        'excursion_limpia_m': [excursion_limpia(M_I, m, mm.L_eff, s, x) for x in TASAS_BARRIDO]}
        f['excursion_limpia_formula_m'] = excursion_limpia(M_I, m, mm.L_eff, s, dm_s)
        f['razon_mismo'] = {k: sim_r[k] / tasa_limpia(M_I, m, mm.L_eff, s, objetivos[k])
                            for k in ('precision', 'seguridad')}
        print(f"{d['corto']:4s} limpia {dm_s:.2f} kg/s (excursión {exc:.3f} m, sube {r0['dh_arriba_m']:.4f} m)")
        for k in objetivos:
            print(f"     reactivo {k:10s} ({objetivos[k]:.3f} m): sim {sim_r[k]:7.2f} kg/s ({f['razon_sim'][k]:5.2f}x)"
                  f"  ideal {ide_r[k]:7.2f} kg/s ({f['razon_ideal'][k]:5.2f}x)")
    rutas.guardar_datos('reactivo_vs_limpia', {'D_lim_m': p.D_lim_oper, 'm_carga_kg': m, 'L': mm.L_eff,
                                               'tasas_barrido_kgps': TASAS_BARRIDO, 'naves': filas})
    figura(filas, TASAS_BARRIDO)
    for nom, f in filas.items():
        print(f"{f['corto']:4s} mismo desplazamiento: reactivo / limpia {f['razon_mismo']['precision']:.2f}x (1 m), "
              f"{f['razon_mismo']['seguridad']:.2f}x (5 m), {f['razon_sim']['igual']:.1f}x (la excursión de la limpia)")


if __name__ == '__main__':
    main()
