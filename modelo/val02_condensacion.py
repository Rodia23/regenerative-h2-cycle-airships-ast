# -*- coding: utf-8 -*-
# [EN] Condensation bound against published hydrogen-air measurements (Supplementary S4.3, Figure S1).
# =============================================================================
# val02_condensacion.py — VAL-02: la ecuación de condensación del paper 1 (Ec. 5)
# contra un sistema H2-aire publicado con mediciones:
#   Zakhvatkin, Schechter, Avrahami, "Water Recuperation from Hydrogen Fuel Cell
#   during Aerial Mission", Energies 15 (2022) 6848.
# El escape de una pila PEM de cátodo abierto es H2 + aire con exceso lambda, como
# nuestra combustión: el agua que se puede recuperar al enfriar a T sale del mismo
# balance de Dalton que la Ec. 5, con N = moles de gas seco que salen por mol de agua.
# A. Su figura 4 (fracción recuperable vs T para varios lambda, aire seco) con
#    nuestro código (IAPWS-IF97) frente a sus ecuaciones (Tetens).
# B. Sus seis experimentos (tabla 1 y apéndice A): agua esperada con nuestro código y
#    sus flujos medidos, frente a su valor esperado y al agua realmente recogida.
# =============================================================================

import math
import numpy as np
import parametros as p
import rutas
from exp05_condensacion_tiempo_real import p_sat

X_O2_ZAK = 0.21                      # fracción molar de O2 que usa el paper
T_MIN = 60.0                         # [min] duración de cada experimento (apéndice A)
AH_IN = 0.005                        # [g/L] humedad absoluta del aire de ventilación (apéndice A)
RHO_H2 = 0.1227                      # [g/L] H2 a 1.5 bar y 25 °C (apéndice A)
RHO_AIRE = 1.2                       # [g/L] (apéndice A)
M_AIRE = 28.97                       # [g/mol] (apéndice A)
LAMBDAS = (1.0, 2.0, 5.0, 10.0)
# Tabla 1 de Zakhvatkin et al. 2022: aire entrada, aire salida (Ec. 15), H2 [L/min], lambda, T_amb [°C],
# AH de saturación [g/m3], agua esperada [g], agua recogida [g], k esperado [%], k logrado [%]
EXPERIMENTOS = (
    (8.0, 7.6, 0.60, 2.4, 26.9, 25.4, 30.7, 22.0, 72.75, 52.13),
    (7.0, 6.5, 0.60, 2.1, 27.4, 26.1, 29.6, 25.2, 70.64, 60.14),
    (10.0, 9.4, 0.78, 2.4, 25.0, 23.1, 38.3, 27.7, 70.53, 51.01),
    (13.0, 12.2, 0.99, 2.4, 24.7, 22.5, 49.1, 40.8, 70.65, 58.71),
    (15.0, 14.2, 1.12, 2.5, 25.7, 23.8, 65.0, 48.5, 72.38, 54.01),
    (17.0, 16.1, 1.30, 2.4, 26.7, 22.1, 64.8, 52.9, 71.05, 58.00),
)


def e_tetens(T_C):
    """Presión de saturación de Tetens [Pa], la que usa el paper (su Ec. 16)."""
    return 610.78 * math.exp(17.27 * T_C / (T_C + 237.3))


def k_seco(T_C, lam, psat):
    """Fracción recuperable del agua producida, aire seco, gas saturado a T y 1 atm (Ec. 5 con lambda)."""
    N = 0.5 * lam / X_O2_ZAK - 0.5          # gas seco que sale por mol de agua: aire - O2 consumido
    ps = psat(T_C)
    return max(0.0, 1.0 - N * ps / (p.P_atmos - ps))


def experimento(fila):
    V_in, V_out_zak, V_H2, lam, T, AH_sat_zak, esp_zak, recog, k_esp_zak, k_log = fila
    n_H2 = V_H2 * RHO_H2 / (p.M_H2 * 1e3)                 # mol/min
    n_aire = V_in * RHO_AIRE / M_AIRE                      # mol/min
    producida = n_H2 * p.M_H2O * 1e3 * T_MIN               # g
    entra = AH_IN * V_in * T_MIN                           # g
    seco_sale = n_aire - 0.5 * n_H2                        # mol/min de gas seco que sale
    ps = p_sat(T)
    sale = seco_sale * ps / (p.P_atmos - ps) * p.M_H2O * 1e3 * T_MIN   # g de vapor que sale saturado
    esperada = producida + entra - sale
    # su apéndice A con SUS números (aire de salida y AH de saturación de la tabla 1): coherencia de la tabla
    esperada_zak_recalc = producida + entra - V_out_zak * AH_sat_zak / 1e3 * T_MIN
    return {'lambda': lam, 'T_amb_C': T, 'producida_g': producida, 'entra_g': entra, 'sale_g': sale,
            'esperada_zak_recalc_g': esperada_zak_recalc,
            'esperada_g': esperada, 'k_esperado': esperada / (producida + entra),
            'esperada_zak_g': esp_zak, 'k_esperado_zak': k_esp_zak / 100.0,
            'recogida_g': recog, 'k_logrado': k_log / 100.0,
            'lambda_calc': n_aire * X_O2_ZAK / (0.5 * n_H2),
            'desvio_vs_zak': esperada / esp_zak - 1.0, 'medido_sobre_modelo': recog / esperada}


def main():
    T = np.linspace(-20.0, 40.0, 121)
    curvas = {lam: {'if97': [k_seco(t, lam, p_sat) for t in T], 'tetens': [k_seco(t, lam, e_tetens) for t in T]}
              for lam in LAMBDAS}
    dif_max = max(abs(a - b) for lam in LAMBDAS for a, b in zip(curvas[lam]['if97'], curvas[lam]['tetens']))
    exps = [experimento(f) for f in EXPERIMENTOS]
    rutas.guardar_datos('val02_condensacion', {
        'fuente': 'Zakhvatkin, Schechter, Avrahami, Energies 15 (2022) 6848',
        'T_C': list(T), 'curvas': {str(k): v for k, v in curvas.items()}, 'dif_max_curvas': dif_max,
        'experimentos': exps})
    graficar(T, curvas, exps)
    print(f'A. figura 4: diferencia máxima IF97 vs Tetens = {100*dif_max:.2f} puntos de %')
    print(f'B. {"lam":>4}{"T":>6}{"esp. nuestro":>13}{"esp. Zak":>9}{"desvío":>8}{"recogida":>9}{"medido/modelo":>14}')
    for x in exps:
        print(f'   {x["lambda"]:>4.1f}{x["T_amb_C"]:>6.1f}{x["esperada_g"]:>13.1f}{x["esperada_zak_g"]:>9.1f}'
              f'{100*x["desvio_vs_zak"]:>7.1f}%{x["recogida_g"]:>9.1f}{100*x["medido_sobre_modelo"]:>13.0f}%')


def graficar(T, curvas, exps):
    import matplotlib.pyplot as plt
    import plot_style  # noqa
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    cols = {1.0: '#0072B2', 2.0: '#009E73', 5.0: '#E69F00', 10.0: '#D55E00'}
    # su modelo: banda gruesa y clara; el nuestro: puntos espaciados encima (nota N41 del autor: con
    # línea discontinua fina las dos curvas se confundían en una sola)
    paso = max(1, len(T) // 14)
    for lam, c in curvas.items():
        axes[0].plot(T, 100 * np.array(c['tetens']), color=cols[lam], lw=6, alpha=0.28, solid_capstyle='round')
        axes[0].plot(T, 100 * np.array(c['if97']), color=cols[lam], ls='none', marker='o', markersize=4,
                     markevery=paso, label=f'$\\lambda$ = {lam:g}')
    axes[0].set_xlabel('Temperature after cooling [°C]')
    axes[0].set_ylabel('Recoverable product water [%]')
    # ANTES (revisión 2026-09-30, P-08: la guía pide figuras sin título interno; el caption ya lo dice):
    # axes[0].set_title('Zakhvatkin et al. 2022, Fig. 4 (thick bands) vs this model (dots)', fontsize=9)
    axes[0].legend(fontsize=8, loc='lower left')
    # margen bajo 0 y sobre 100: por encima del punto de rocío la fracción es 0 (tramo plano) y con el eje
    # cortado en 0 la banda y los puntos quedaban partidos por el borde (observación del autor, 2026-09-27)
    axes[0].set_ylim(-4, 104)
    axes[0].axhline(0.0, color='gray', lw=0.6, zorder=0)
    i = np.arange(1, len(exps) + 1)
    axes[1].bar(i - 0.27, [x['esperada_g'] for x in exps], 0.26, color='#0072B2', label='This model')
    # ANTES (revisión 2026-09-30, I-01: el texto compara con su apéndice A recalculado, la figura mostraba su tabla 1):
    # axes[1].bar(i, [x['esperada_zak_g'] for x in exps], 0.26, color='#56B4E9', label='Zakhvatkin et al. (expected)')
    # AHORA: barra = su apéndice A recalculado con sus flujos (la referencia del texto); cruz = su tabla 1, que en los
    # experimentos 5 y 6 no coincide con su propio apéndice (fuentes.py, ZAK22, hallazgos)
    axes[1].bar(i, [x['esperada_zak_recalc_g'] for x in exps], 0.26, color='#56B4E9',
                label='Their Appendix A, recomputed')
    axes[1].bar(i + 0.27, [x['recogida_g'] for x in exps], 0.26, color='#D55E00', label='Measured')
    axes[1].plot(i, [x['esperada_zak_g'] for x in exps], 'x', color='black', ms=7, mew=1.6, zorder=5,
                 label='Their Table 1 (expected)')
    axes[1].set_xticks(i)
    axes[1].set_xlabel('Experiment')
    axes[1].set_ylabel('Water collected in 60 min [g]')
    axes[1].legend(fontsize=8, loc='upper left')
    for ax in axes:
        ax.grid(True, alpha=0.2)
    fig.tight_layout()
    plot_style.guardar_paper(fig, rutas.figura('val02_condensacion'), dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    main()
