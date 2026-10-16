# -*- coding: utf-8 -*-
# [EN] Small releases sized by the installed heat rejection, with its lag (supporting step).
# =============================================================================
# exp10_escalones_condensador.py — EXP-10: ¿aguanta la liberación por escalones
# (EXP-06) cuando el condensador tiene retardo térmico (EXP-09)?
#
# EXP-06 dimensiona el escalón máximo con un condensador que responde al
# instante: el pico de quema del minimax por q no pasa de la potencia instalada.
# EXP-09 mostró que lo que rompe el minimax de 1 t es el retardo. Aquí se simula
# un escalón de ese tamaño (h_max = 1 m, D_lim = 0.10 m, el minimax y el PD del
# artículo) en las cuatro plataformas, con 10 y 100 MW instalados y el vapor que
# no alcanza a condensar venteado (retardo 0, 5 y 30 s) o retenido como gas
# (pulmón, 5 s). Sale: si cumple h_max, cuánta agua queda y cuánto H2 cuesta
# por tonelada frente a ventear H2 directamente.
# =============================================================================

import numpy as np
import parametros as p
import rutas
import condensador
import drop_flujo_optimo_perdrop as sim
from exp09_condensacion_dinamica import resumen

P_LISTA_MW = (10.0, 100.0)
CASOS = (('ventear', 0.0), ('ventear', 5.0), ('ventear', 30.0), ('pulmon', 5.0))
H_ESCALON = 1.0                                  # umbral con el que EXP-06 dimensiona el escalón
T_POST = 600.0


def escalones_exp06():
    return {(x['plataforma'], x['P_instalada_MW']): x['paso_max_kg']
            for x in rutas.leer_datos('exp06_requisitos_minimos')['pasos'] if x['h_m'] == H_ESCALON}


def duracion_quema(s):
    q = np.clip(-np.gradient(s['m_H2'], s['t']), 0.0, None)
    return float(np.sum(np.diff(s['t'])[q[:-1] > 0.05 * q.max()])) if q.max() > 0 else 0.0


def asentado(s):
    """La nave terminó quieta en h0 (si no, el lazo perdió estabilidad y no hay balance final)."""
    v = np.gradient(s['h'], s['t'])
    return bool(abs(s['h'][-1] - p.h0) < 0.05 and abs(v[-1]) < 0.01)


def main():
    pasos = escalones_exp06()
    L_ideal = p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1)) + p.ratio_O2_H2
    L_venteo = L_ideal - (1.0 + p.ratio_O2_H2)
    filas, series = [], {}
    for n, d in p.PLATAFORMAS.items():
        for P in P_LISTA_MW:
            m = pasos[(n, P)]
            sim.usar_nave(n, m)
            for modo, tau in (('ideal', None),) + CASOS:
                c = None if modo == 'ideal' else condensador.config(P * 1e6, tau_s=tau, modo=modo)
                nombre = 'ideal' if c is None else f'{modo}_tau{tau:g}s'
                r = sim.correr(n_entregas=1, toneladas_por_drop=m / 1000.0, graficar=False, verbose=False,
                               t_post=T_POST, cond=c, devolver_serie=True)
                x = resumen(nombre, r, c)
                x.update({'plataforma': n, 'corto': d['corto'], 'P_instalada_MW': P, 'escalon_kg': m,
                          'H2_por_t_kg': 1000.0 * x['H2_quemado_kg'] / m,
                          'dur_quema_s': duracion_quema(r['serie']), 'asentado': asentado(r['serie'])})
                filas.append(x)
                series[(d['corto'], P, nombre)] = r['serie']

    c0 = condensador.config(None)
    rutas.guardar_datos('exp10_escalones_condensador', {
        'h_escalon_m': H_ESCALON, 'P_instalada_MW': list(P_LISTA_MW), 'phi_max': c0['phi'],
        'H2_por_t_ideal_kg': 1000.0 / L_ideal, 'H2_por_t_venteo_kg': 1000.0 / L_venteo,
        'T_post_s': T_POST, 'casos': filas})
    graficar(series, filas, 1000.0 / L_venteo)
    tabla(filas)

    print(f'{"plat":<6}{"P MW":>6}{"escalón":>9}{"caso":>18}{"H2 kg":>8}{"kg/t":>7}{"% ret":>7}'
          f'{"arriba":>8}{"abajo":>8}{"cumple":>8}{"quema s":>9}')
    for x in filas:
        print(f'{x["corto"]:<6}{x["P_instalada_MW"]:>6.0f}{x["escalon_kg"]:>9.0f}{x["caso"]:>18}'
              f'{x["H2_quemado_kg"]:>8.2f}{x["H2_por_t_kg"]:>7.1f}{100*x["frac_retenida"]:>7.0f}'
              f'{x["dh_arriba_m"]:>8.2f}{x["dh_abajo_m"]:>8.2f}{str(x["cumple_h_max"]):>8}{x["dur_quema_s"]:>9.1f}')
    print(f'\nReferencias: ideal {1000/L_ideal:.1f} kg/t; ventear H2 directamente {1000/L_venteo:.1f} kg/t')


def tabla(filas):
    """Por plataforma e instalación: H2 por tonelada [kg/t] y agua retenida [%]; † = no cumple h_max."""
    casos = list(dict.fromkeys(x['caso'] for x in filas))
    lin = [r'\begin{tabular}{lrrr' + 'r' * len(casos) + '}', r'\toprule',
           r'Plataforma & $P_\text{cap}$ [MW] & Escalón [kg] & Quema [s] & Ideal & Retardo 0~s & Retardo 5~s'
           r' & Retardo 30~s & Pulmón 5~s \\', r'\midrule']
    claves = list(dict.fromkeys((x['corto'], x['P_instalada_MW']) for x in filas))
    for corto, P in claves:
        fila = {x['caso']: x for x in filas if x['corto'] == corto and x['P_instalada_MW'] == P}
        x0 = fila['ideal']
        celdas = [f"{fila[c]['H2_por_t_kg']:.0f} ({100 * min(fila[c]['frac_retenida'], 1.0):.0f})"
                  + ('' if fila[c]['cumple_h_max'] else r'$^\dagger$') for c in casos]
        lin.append(f"{corto} & {P:.0f} & {x0['escalon_kg']:.0f} & {x0['dur_quema_s']:.0f} & "
                   + ' & '.join(celdas) + r' \\')
    lin += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('exp10_escalones_condensador'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lin) + '\n')


def graficar(series, filas, kg_t_venteo):
    import matplotlib.pyplot as plt
    import plot_style  # noqa
    estilos = {'ideal': ('black', '-'), 'ventear_tau0s': ('#0072B2', '-'), 'ventear_tau5s': ('#E69F00', '--'),
               'ventear_tau30s': ('#D55E00', ':'), 'pulmon_tau5s': ('#009E73', '-.')}
    etiquetas = {'ideal': 'Instantaneous condensation (ideal)', 'ventear_tau0s': 'Vented, lag 0 s',
                 'ventear_tau5s': 'Vented, lag 5 s', 'ventear_tau30s': 'Vented, lag 30 s',
                 'pulmon_tau5s': 'Held as gas, lag 5 s'}
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, (corto, P) in zip(axes[:2], ((p.NAVE_BASE['corto'], 100.0), ('ZNT', 10.0))):
        esc = next(x['escalon_kg'] for x in filas if x['corto'] == corto and x['P_instalada_MW'] == P)
        for nombre, (col, ls) in estilos.items():
            s = series[(corto, P, nombre)]
            ax.plot(s['t'], s['h'] - p.h0, color=col, ls=ls, lw=1.7, label=etiquetas[nombre])
        ax.axhspan(-H_ESCALON, H_ESCALON, color='gray', alpha=0.12, lw=0)
        ax.set_xlim(0, 300)
        ax.set_title(f'{corto}, {P:g} MW installed, step {esc:.0f} kg', fontsize=10)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Altitude change [m]')
        ax.grid(True, alpha=0.2)
    axes[0].legend(fontsize=8, loc='lower right')
    ax = axes[2]
    cortos = list(dict.fromkeys(x['corto'] for x in filas))
    nombres = list(estilos)
    ancho = 0.8 / (len(nombres) * len(P_LISTA_MW))
    for i, corto in enumerate(cortos):
        k = 0
        for P in P_LISTA_MW:
            for nombre in nombres:
                x = next(y for y in filas if y['corto'] == corto and y['P_instalada_MW'] == P and y['caso'] == nombre)
                col, _ = estilos[nombre]
                ax.bar(i - 0.4 + (k + 0.5) * ancho, x['H2_por_t_kg'], ancho * 0.9, color=col,
                       alpha=1.0 if P == P_LISTA_MW[0] else 0.55, hatch='' if x['cumple_h_max'] else '//',
                       edgecolor='white' if x['cumple_h_max'] else 'black', lw=0.4)
                k += 1
    ax.axhline(kg_t_venteo, color='gray', ls='--', lw=1)
    ax.text(len(cortos) - 0.5, kg_t_venteo, ' venting H$_2$', va='bottom', ha='right', fontsize=8, color='gray')
    ax.set_xticks(range(len(cortos)))
    ax.set_xticklabels(cortos)
    ax.set_ylabel('H$_2$ burned per tonne released [kg/t]')
    ax.set_title('Cost per tonne (solid: 10 MW, light: 100 MW; hatched: misses $h_{max}$)', fontsize=9)
    ax.set_yscale('log')
    ax.grid(True, axis='y', alpha=0.2)
    fig.suptitle('Stepwise release sized by installed condensation power, with condenser dynamics', fontsize=11)
    fig.tight_layout()
    fig.savefig(rutas.figura('exp10_escalones_condensador'), dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
