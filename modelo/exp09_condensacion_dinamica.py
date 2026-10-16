# -*- coding: utf-8 -*-
# [EN] Simulations with the dynamics of the heat-rejection unit (Section 5.1, Figure 3).
# =============================================================================
# exp09_condensacion_dinamica.py — EXP-09: las simulaciones del ARTÍCULO con la
# dinámica del condensador (condensador.py), en vez de la condensación
# instantánea (secciones 5.1 y 9 del paper).
#
# Maniobra: la quema limpia del artículo tal cual (caso base, 1 t, h_max = 1 m,
# D_lim = 0.10 m, su controlador PD). Casos:
#   - ideal: condensación instantánea, 9 kg de agua por kg de H2 (el artículo);
#   - condensador con capacidad P_cap (∞ = solo el límite de equilibrio phi,
#     1000, 100 y 10 MW), retardo térmico tau, y el vapor que no alcanza:
#       'ventear' (sale con el N2) o 'pulmon' (gas dentro de la envolvente que
#       empuja hacia arriba hasta que se condensa). OJO (hallazgo M3): vapor puro a
#       T_pulmon y 1 atm no es un estado físico (se condensaría); como gas iría con su
#       N2, que también sustenta: el pulmón SUBESTIMA la subida. Fuera del paper 1
#       desde 2026-09-29 (solo como contraste interno, con esta salvedad);
#   - sensibilidad al retardo (0, 5, 30 s) con capacidad sobrada.
# Operación (desde 2026-09-27, nota N22 del autor): la unidad se PRE-ENCIENDE n_tau_arranque·tau antes
# de la quema previa (condensador.config, prearranque=True), como describe el paper (sec. 2.3); antes
# arrancaba desde cero con la quema. Ese caso queda como una sola sensibilidad: 'sin_prearranque_inf'.
# Además: las cuatro plataformas con un condensador de 100 MW (ventear).
# =============================================================================

import io
import sys
import numpy as np
import parametros as p
import rutas
import condensador
import drop_flujo_optimo_perdrop as sim

P_LISTA_MW = (None, 1000.0, 100.0, 10.0)        # None = sin límite de capacidad
TAU_LISTA = (0.0, 5.0, 30.0)
T_POST = 600.0                                  # más largo que el artículo: el pulmón tarda en vaciarse


def correr_pf(cond, t_post=T_POST):
    sim.usar_nave(p.nombre_dirigible, 1000.0)      # caso base del artículo
    return sim.correr(n_entregas=1, toneladas_por_drop=1.0, graficar=False, verbose=False, t_post=t_post,
                      cond=cond, devolver_serie=True)


def resumen(nombre, r, cond):
    # H2 quemado desde los ESTADOS del integrador: el agua solo se produce al quemar, así que
    # quemado = (retenida + venteada + retenida como vapor) / (M_H2O/M_H2). Antes se derivaba m_H2 en la
    # malla de salida y se sumaba lo positivo: con quema y reinyección alternando dentro de un
    # intervalo eso veía el neto (retención > 100 % en el Zeppelin NT) y dependía de la malla (VER-02).
    ret = r['m_agua_fin']
    producida = ret + r['m_agua_venteada'] + r['m_vapor_pulmon_fin']
    quemado = producida / (1.0 + p.ratio_O2_H2)
    return {'caso': nombre, 'P_cap_MW': None if cond is None or not np.isfinite(cond['P_cap']) else cond['P_cap'] / 1e6,
            'tau_s': None if cond is None else cond['tau'], 'modo': 'ideal' if cond is None else cond['modo'],
            'prearranque': None if cond is None else cond['prearranque'],
            'H2_quemado_kg': quemado, 'H2_neto_kg': r['m_H2_neto'], 'agua_producida_kg': producida,
            'agua_retenida_kg': ret, 'agua_venteada_kg': r['m_agua_venteada'],
            'vapor_pulmon_fin_kg': r['m_vapor_pulmon_fin'],
            'frac_retenida': ret / producida if producida > 0 else float('nan'),
            'dh_arriba_m': r['dh_arriba_m'], 'dh_abajo_m': r['dh_abajo_m'], 'h_max_m': r['h_max_m'],
            'cumple_h_max': bool(max(r['dh_arriba_m'], -r['dh_abajo_m']) <= r['h_max_m']),
            'dm_opt_kgps': r['dm_opt'], 'P_calor_pico_MW': r['dm_opt'] * (cond['q'] if cond else 0) / 1e6}


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    casos, series = [], {}
    r = correr_pf(None)
    casos.append(resumen('ideal', r, None)); series['ideal'] = r['serie']
    t_anticipo = r['t_anticipo']          # la misma quema en todas las corridas del caso base
    for modo in ('ventear', 'pulmon'):
        for P in P_LISTA_MW:
            c = condensador.config(None if P is None else P * 1e6, tau_s=p.tau_condensador, modo=modo)
            nombre = f'{modo}_{"inf" if P is None else int(P)}MW'
            r = correr_pf(c)
            casos.append(resumen(nombre, r, c)); series[nombre] = r['serie']
    for tau in TAU_LISTA:
        c = condensador.config(1000e6, tau_s=tau, modo='ventear')
        nombre = f'tau_{tau:g}s'
        r = correr_pf(c)
        casos.append(resumen(nombre, r, c)); series[nombre] = r['serie']
    # sensibilidad: la unidad NO se pre-enciende; su potencia sigue a la quema desde cero con retardo
    c = condensador.config(None, tau_s=p.tau_condensador, modo='ventear', prearranque=False)
    r = correr_pf(c)
    casos.append(resumen('sin_prearranque_inf', r, c)); series['sin_prearranque_inf'] = r['serie']

    plataformas = []
    for n, d in p.PLATAFORMAS.items():
        sim.usar_nave(n, 1000.0)
        for cond in (None, condensador.config(100e6, modo='ventear')):
            r = sim.correr(n_entregas=1, toneladas_por_drop=1.0, graficar=False, verbose=False, t_post=T_POST,
                           cond=cond, devolver_serie=True)
            x = resumen('ideal' if cond is None else 'ventear_100MW', r, cond)
            x.update({'plataforma': n, 'corto': d['corto']})
            plataformas.append(x)

    c0 = condensador.config(None)
    rutas.guardar_datos('exp09_condensacion_dinamica', {
        'phi_max': c0['phi'], 'q_MJkg': c0['q'] / 1e6, 'h_agua_MJkg': c0['h_agua'] / 1e6,
        'tau_base_s': p.tau_condensador, 'T_pulmon_C': p.T_pulmon, 'T_post_s': T_POST,
        'casos_base': casos, 'plataformas': plataformas})
    graficar(series, casos)
    graficar_articulo(series, t_anticipo)

    print(f'phi_max = {c0["phi"]:.3f}, q = {c0["q"]/1e6:.1f} MJ/kg, h_agua = {c0["h_agua"]/1e6:.2f} MJ/kg de agua')
    print(f'{"caso":<18}{"H2 kg":>8}{"retenida":>10}{"venteada":>10}{"pulmón":>8}{"% ret":>7}'
          f'{"arriba m":>10}{"abajo m":>9}{"cumple":>8}')
    for x in casos:
        print(f'{x["caso"]:<18}{x["H2_quemado_kg"]:>8.1f}{x["agua_retenida_kg"]:>10.1f}{x["agua_venteada_kg"]:>10.1f}'
              f'{x["vapor_pulmon_fin_kg"]:>8.1f}{100*x["frac_retenida"]:>7.0f}{x["dh_arriba_m"]:>10.2f}'
              f'{x["dh_abajo_m"]:>9.2f}{str(x["cumple_h_max"]):>8}')
    print('\nCuatro plataformas, 1 t (ideal vs condensador de 100 MW, ventear):')
    for x in plataformas:
        print(f'  {x["plataforma"]:<22}{x["caso"]:<15} H2 {x["H2_quemado_kg"]:6.1f} kg  retenida {100*x["frac_retenida"]:4.0f} %'
              f'  arriba {x["dh_arriba_m"]:5.2f} m  abajo {x["dh_abajo_m"]:5.2f} m  cumple {x["cumple_h_max"]}')


def graficar(series, casos):
    import matplotlib.pyplot as plt
    import plot_style  # noqa
    fig, axes = plt.subplots(4, 2, figsize=(14, 14), sharex=True)
    estilos = {'ideal': ('black', '-'), 'inf': ('#0072B2', '-'), '1000': ('#009E73', '--'),
               '100': ('#E69F00', '-.'), '10': ('#D55E00', ':')}
    for col, modo in enumerate(('ventear', 'pulmon')):
        for nombre in ['ideal'] + [f'{modo}_{"inf" if P is None else int(P)}MW' for P in P_LISTA_MW]:
            s = series[nombre]
            clave = 'ideal' if nombre == 'ideal' else nombre.split('_')[1].replace('MW', '')
            col_, ls = estilos[clave]
            etiqueta = ('Instantaneous condensation (article)' if nombre == 'ideal'
                        else f'Condenser {"unlimited" if clave == "inf" else clave + " MW"}')
            t = s['t']
            axes[0, col].plot(t, s['h'], color=col_, ls=ls, lw=1.8, label=etiqueta)
            axes[1, col].plot(t, np.clip(-np.gradient(s['m_H2'], t), 0, None), color=col_, ls=ls, lw=1.6)
            axes[2, col].plot(t, s['agua'], color=col_, ls=ls, lw=1.8)
            if nombre != 'ideal':
                axes[3, col].plot(t, s['vent'] if modo == 'ventear' else s['pulmon'], color=col_, ls=ls, lw=1.8)
        axes[0, col].axhline(p.h0, color='gray', ls=':', lw=0.8)
        axes[0, col].set_title('Uncondensed vapour VENTED with the N$_2$' if modo == 'ventear'
                               else 'Uncondensed vapour HELD as gas in the envelope (buffer)', fontsize=10)
        axes[0, col].set_ylabel('Altitude [m]')
        axes[1, col].set_ylabel('H$_2$ burn rate [kg/s]')
        axes[2, col].set_ylabel('Water kept aboard [kg]')
        axes[3, col].set_ylabel('Vented water [kg]' if modo == 'ventear' else 'Vapour held as gas [kg]')
        axes[3, col].set_xlabel('Time [s]')
        axes[0, col].legend(fontsize=8, loc='lower right')
        for ax in axes[:, col]:
            ax.grid(True, alpha=0.2)
            ax.set_xlim(0, 200)
    fig.suptitle(f'Article base case ({p.nombre_dirigible}, 1 t, sized burn): condensation dynamics '
                 f'(thermal lag {p.tau_condensador:g} s)', fontsize=11)
    fig.tight_layout()
    fig.savefig(rutas.figura('exp09_condensacion_dinamica'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def graficar_articulo(series, t_anticipo):
    """Versión liviana para el paper 1: altitud y agua retenida, vapor venteado. El pulmón sale del paper (hallazgo
    M3, autor 2026-09-29: vapor puro a 70 °C y 1 atm no es un estado físico).
    Con el pre-encendido (nota N22) el retardo casi no cuenta: las etiquetas no lo nombran (tau va en el pie)
    y dos líneas verticales marcan la orden de arranque de la unidad (punteada) y la liberación (guiones)."""
    import matplotlib.pyplot as plt
    import plot_style  # noqa
    curvas = (('ideal', 'Instantaneous condensation (ideal bound)', 'black', '-'),
              ('ventear_infMW', 'Unlimited power', '#0072B2', '-'),
              ('ventear_100MW', '100 MW', '#E69F00', '--'),
              ('ventear_10MW', '10 MW', '#D55E00', ':'))
    t_rel = series['ideal']['t_drops'][0]
    t_arr = t_rel - t_anticipo - p.n_tau_arranque * p.tau_condensador
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for clave, etiqueta, col, ls in curvas:
        s = series[clave]
        axes[0].plot(s['t'], s['h'] - p.h0, color=col, ls=ls, lw=1.7, label=etiqueta)
        axes[1].plot(s['t'], s['agua'], color=col, ls=ls, lw=1.7)
    for ax in axes:
        ax.axvline(t_arr, color='gray', ls=':', lw=1.0)
        ax.axvline(t_rel, color='gray', ls='--', lw=1.0)
    axes[0].axhspan(-1.0, 1.0, color='gray', alpha=0.12, lw=0)
    axes[0].set_ylabel('Altitude change [m]')
    axes[1].set_ylabel('Water retained on board [kg]')
    for ax in axes:
        ax.set_xlabel('Time [s]')
        ax.set_xlim(0, 200)
        ax.grid(True, alpha=0.2)
    axes[0].legend(fontsize=8, loc='upper right')
    fig.tight_layout()
    fig.savefig(rutas.figura('exp09_articulo'), dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
