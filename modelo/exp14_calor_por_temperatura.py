# -*- coding: utf-8 -*-
# [EN] Heat by temperature band; cold band at 1 and 14 atm (Section 9, Appendix C).
# =============================================================================
# exp14_calor_por_temperatura.py — EXP-14: no todo el calor es igual de difícil.
# (Surgió de la idea de la "tina" del autor, 2026-09-25.)
#
# El calor a rechazar por kg de H2 (q, EXP-05) se reparte en dos bandas:
#   - CALIENTE: todo lo que sale por encima del punto de rocío (gas caliente);
#   - FRÍA: lo que sale entre el punto de rocío y el fin de la condensación
#     (casi todo latente), cerca de la temperatura del aire.
# Cada banda se dimensiona con su propio kW/kg:
#   - caliente: intercambiador del ciclo de turbina de H2, FlyZero 20 kW/kg (2026),
#     solo "> 20" en 2035 -> se usa 20 en los dos años (decisión del autor, 2026-09-25: el 30
#     no está impreso en la fuente; ver fuentes.FZ22);
#   - fría: radiador de aire cuyo kW/kg es proporcional a la diferencia de
#     temperatura con el aire, calibrado con FlyZero: 10 kW/kg rechazando a 100 °C
#     (2026) y 20 kW/kg a 160 °C (2035), con aire a la temperatura de diseño.
#     SUPUESTO declarado: kW/kg = k (T - T_aire), integrado a lo largo de la banda
#     con el aire de enfriamiento en contracorriente (T_aire_exterior -> +dT_aire).
# Presión: condensar la misma fracción phi a mayor presión sube la banda fría.
# Tina (condensador de contacto directo): cuánta agua haría falta para guardar el
# calor frío de una tonelada, y cuánto aporta el propio lastre acumulado.
# =============================================================================

import numpy as np
import parametros as p
import rutas
from exp05_condensacion_tiempo_real import p_sat, N_INERTE, AIRE_ESTEQ, w_comp_escape

PRESIONES_ATM = (1.0, 2.0, 5.0, 10.0, 20.0)
ETAS = (0.0, p.rend_gen)
KWKG_CALIENTE = (20.0, 20.0)            # FlyZero p. 7: intercambiadores del ciclo de turbina de H2, 2026 / 2035 (">20")
T_RECHAZO_FLYZERO = (100.0, 160.0)      # °C a los que FlyZero da 10 y 20 kW/kg (2026 / 2035)
N_PASOS = 4000


def L_v(T_C):
    """Calor latente del agua [J/kg], ajuste lineal 0-100 °C (IAPWS: 2501 kJ/kg a 0 °C, 2257 a 100 °C)."""
    return 2.5009e6 - 2.4436e3 * T_C


def T_sat(pv):
    lo, hi = -20.0, 370.0
    for _ in range(80):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if p_sat(m) < pv else (lo, m)
    return 0.5 * (lo + hi)


def phi_T(T_C, P):
    ps = p_sat(T_C)
    return max(0.0, 1.0 - N_INERTE * ps / (P - ps))


def banda_fria(P, phi_obj):
    """Enfriar el escape del rocío al fin de la condensación (phi_obj) a presión P.
    Devuelve temperaturas, la curva (T, calor acumulado) y el calor de la banda [J/kg H2]."""
    r = 1.0 + p.ratio_O2_H2
    inerte = AIRE_ESTEQ - p.ratio_O2_H2
    x_w = 1.0 / (1.0 + N_INERTE)
    T_r = T_sat(x_w * P)
    T_f = T_sat((1.0 - phi_obj) / (1.0 - phi_obj + N_INERTE) * P)
    T = np.linspace(T_r, T_f, N_PASOS)
    phi = np.array([phi_T(t, P) for t in T])
    # calor liberado desde el rocío hasta T: sensible del gas + sensible del agua + latente de lo condensado
    Q = (inerte * p.cp_N2 + r * p.cp_vapor) * (T_r - T) + r * phi * L_v(T)
    return {'T_rocio_C': T_r, 'T_fin_C': T_f, 'T': T, 'Q': Q, 'Q_banda': float(Q[-1])}


def masa_banda_fria(b, k, T_aire0=None):
    """Masa de radiador por W de calor de la banda [kg/W]: suma de (fracción de calor) / (k (T - T_aire)),
    con el aire en contracorriente: entra a T_aire0 (por defecto T_aire_exterior) por la punta fría y sale
    dT_aire_enfr más caliente."""
    T_aire0 = p.T_aire_exterior if T_aire0 is None else T_aire0
    Qb = b['Q_banda']
    T = b['T'][::-1]                                    # de la punta fría (fin de la condensación) al rocío
    Qc = Qb - b['Q'][::-1]                              # calor acumulado desde la punta fría: 0 -> Qb
    T_aire = T_aire0 + p.dT_aire_enfr * Qc / Qb
    frac = np.diff(Qc) / Qb
    dT = 0.5 * ((T[1:] - T_aire[1:]) + (T[:-1] - T_aire[:-1]))
    if np.any(dT <= 0.0):
        return float('inf')
    return float(np.sum(frac / (k * 1e3 * dT)))


def q_total(eta, phi, T_out_C, T_in_C=None):
    """Mismo balance que q_rech (EXP-05) con la temperatura de salida (y, opcional, la de entrada del H2 y
    el aire; por defecto T_aire_exterior) como parámetros."""
    r = 1.0 + p.ratio_O2_H2
    T_in_C = p.T_aire_exterior if T_in_C is None else T_in_C
    entra = (p.cp_H2 + AIRE_ESTEQ * p.cp_aire) * (T_in_C - p.T_ref_comb)
    sale = ((AIRE_ESTEQ - p.ratio_O2_H2) * p.cp_N2
            + r * (phi * p.cp_agua_liq + (1.0 - phi) * p.cp_vapor)) * (T_out_C - p.T_ref_comb)
    return (1.0 - eta) * p.Q_comb_H2 + phi * r * p.L_v_agua + entra - sale


def main():
    phi_obj = phi_T(p.T_salida_cond, p.P_atmos)          # la fracción del caso base (40 °C, 1 atm)
    k_fria = tuple(kw / (T - p.T_aire_exterior) for kw, T in zip((10.0, 20.0), T_RECHAZO_FLYZERO))
    k_rho = p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1)) - 1.0
    B_t = 1000.0 / (k_rho + (1.0 + p.ratio_O2_H2) * phi_obj)
    filas = []
    for P_atm in PRESIONES_ATM:
        P = P_atm * p.P_atmos
        b = banda_fria(P, phi_obj)
        W = w_comp_escape(P, b['T_fin_C'])          # mínimo: solo la fase gaseosa (M-1, auditoría 3)
        for eta in ETAS:
            qt = q_total(eta, phi_obj, b['T_fin_C'])
            q_fria = b['Q_banda'] + W                    # la compresión también sale como calor (conservador: en la banda fría)
            q_cal = qt - b['Q_banda']
            for i, ext in enumerate(('2026', '2035')):
                m_f = masa_banda_fria(b, k_fria[i])            # kg por W de la banda fría
                m_c = 1.0 / (KWKG_CALIENTE[i] * 1e3)           # kg por W de la banda caliente
                kg_por_kW = (q_cal * m_c + q_fria * m_f) / (qt + W) * 1e3
                filas.append({'P_atm': P_atm, 'eta': eta, 'tec': ext, 'T_rocio_C': b['T_rocio_C'],
                              'T_fin_C': b['T_fin_C'], 'q_total_MJkg': (qt + W) / 1e6,
                              'q_caliente_MJkg': q_cal / 1e6, 'q_fria_MJkg': q_fria / 1e6,
                              'frac_fria': q_fria / (qt + W), 'W_comp_MJkg': W / 1e6,
                              'kWkg_fria_efectivo': 1.0 / (m_f * 1e3),
                              'kg_por_kW': kg_por_kW,
                              'kg_por_kW_articulo': 1.0 / p.kWkg_radiador[i],
                              'razon_vs_articulo': kg_por_kW * p.kWkg_radiador[i]})
    # Tina: guardar el calor frío de una tonelada en agua (1 atm, eta = 0)
    b1 = banda_fria(p.P_atmos, phi_obj)
    Q_frio_t = B_t * b1['Q_banda']
    tina = {'Q_frio_por_t_GJ': Q_frio_t / 1e9,
            'agua_t_por_t': {dT: Q_frio_t / (p.cp_agua_liq * dT) / 1e3 for dT in (5.0, 10.0)},
            'lastre_por_t_t': B_t * (1.0 + p.ratio_O2_H2) * phi_obj / 1e3,
            'frac_lastre_10K': B_t * (1.0 + p.ratio_O2_H2) * phi_obj * p.cp_agua_liq * 10.0 / Q_frio_t,
            'phi_tina': {T: phi_T(T, p.P_atmos) for T in (40.0, 45.0, 60.0, 73.0)}}
    # (2026-10-13) la banda fría de los dos casos de retención del paper (auditoría, objeción 5.1): 1 atm (phi_max a
    # 40 °C) y la presión que retiene phi objetivo (0.99) a 40 °C (la de EXP-06, ~14 atm). Fracción = banda / q(0, phi),
    # el q de la ec. (4) del paper, sin el trabajo de compresión (va aparte, W_comp)
    P_obj = rutas.leer_datos('exp06_requisitos_minimos')['P_req_atm'] * p.P_atmos
    bandas = {}
    for clave, P_b, phi_b in (('1atm', p.P_atmos, phi_obj), ('objetivo', P_obj, p.phi_objetivo)):
        bb = banda_fria(P_b, phi_b)
        qb = q_total(0.0, phi_b, bb['T_fin_C'])
        bandas[clave] = {'P_atm': P_b / p.P_atmos, 'phi': phi_b, 'T_rocio_C': bb['T_rocio_C'], 'T_fin_C': bb['T_fin_C'],
                         'Q_banda_MJkg': bb['Q_banda'] / 1e6, 'q_MJkg': qb / 1e6, 'frac': bb['Q_banda'] / qb,
                         'W_comp_MJkg': w_comp_escape(P_b, bb['T_fin_C']) / 1e6}
    rutas.guardar_datos('exp14_calor_por_temperatura', {
        'phi_obj': phi_obj, 'k_fria_kWkg_por_K': k_fria, 'T_rechazo_flyzero_C': T_RECHAZO_FLYZERO,
        'kWkg_caliente': KWKG_CALIENTE, 'H2_por_t_kg': B_t, 'filas': filas, 'tina': tina, 'bandas': bandas})
    graficar(b1, filas)
    print(f'phi objetivo {phi_obj:.3f}; k banda fria {k_fria[0]:.3f} / {k_fria[1]:.3f} kW/kg por K')
    print(f'{"P atm":>6}{"eta":>6}{"tec":>6}{"rocio":>7}{"fin":>6}{"q":>7}{"fria":>7}{"%fria":>7}{"kW/kg fria":>11}'
          f'{"kg/kW":>8}{"art":>7}{"x art":>7}')
    for f in filas:
        print(f'{f["P_atm"]:>6.0f}{f["eta"]:>6.2f}{f["tec"]:>6}{f["T_rocio_C"]:>7.0f}{f["T_fin_C"]:>6.0f}'
              f'{f["q_total_MJkg"]:>7.1f}{f["q_fria_MJkg"]:>7.1f}{100*f["frac_fria"]:>7.0f}{f["kWkg_fria_efectivo"]:>11.2f}'
              f'{f["kg_por_kW"]:>8.3f}{f["kg_por_kW_articulo"]:>7.3f}{f["razon_vs_articulo"]:>7.2f}')
    print('tina:', {k: (round(v, 3) if isinstance(v, float) else v) for k, v in tina.items()})


def graficar(b1, filas):
    import matplotlib.pyplot as plt
    import plot_style  # noqa
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    ax = axes[0]
    ax.plot(b1['Q'] / 1e6, b1['T'], color='#0072B2', lw=2, label='Exhaust (condensing band)')
    T_aire = p.T_aire_exterior + p.dT_aire_enfr * b1['Q'] / b1['Q'][-1]
    ax.plot(b1['Q'] / 1e6, T_aire[::-1], color='#E69F00', lw=1.5, ls='--', label='Cooling air (counterflow)')
    ax.set_xlabel('Heat released below the dew point [MJ per kg H$_2$]')
    ax.set_ylabel('Temperature [°C]')
    ax.set_title('1 atm: the cold band is a few degrees above the air', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)
    ax = axes[1]
    for eta, col in ((0.0, '#0072B2'), (p.rend_gen, '#D55E00')):
        for tec, ls in (('2026', '-'), ('2035', '--')):
            xs = [f for f in filas if f['eta'] == eta and f['tec'] == tec]
            ax.plot([f['P_atm'] for f in xs], [f['kg_por_kW'] for f in xs], color=col, ls=ls, marker='o', ms=4,
                    label=f'$\\eta$ = {eta:.2f}, {tec}')
    for i, tec in enumerate(('2026', '2035')):
        ax.axhline(1.0 / p.kWkg_radiador[i], color='gray', ls=':' if i else '-.', lw=1)
    ax.text(20, 1.0 / p.kWkg_radiador[0], ' article, 2026', fontsize=8, color='gray', va='bottom', ha='right')
    ax.text(20, 1.0 / p.kWkg_radiador[1], ' article, 2035', fontsize=8, color='gray', va='bottom', ha='right')
    ax.set_xscale('log')
    ax.set_xlabel('Exhaust pressure during condensation [atm]')
    ax.set_ylabel('Heat-rejection mass [kg per kW rejected]')
    ax.set_title('Two-band sizing vs a single radiator figure', fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(rutas.figura('exp14_calor_por_temperatura'), dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
