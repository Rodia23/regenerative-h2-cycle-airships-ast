# -*- coding: utf-8 -*-
# [EN] Heat to reject and condensable fraction of the exhaust (Equations 4-5, Section 3.3).
# =============================================================================
# exp05_condensacion_tiempo_real.py — EXP-05 (análisis, SIN simular): la
# condensación en tiempo real como BRECHA DE DISEÑO. Pone números y márgenes
# para trabajos futuros; no busca soluciones.
#
# Supuesto acordado (EXP-02): el agua queda a bordo como líquido al instante.
# Para eso hay que enfriar TODO el escape (N2 + vapor) por debajo del rocío:
#   q_rech = (1 − η)·PCI + φ·9·L_v      [J por kg de H2 quemado]
#   η = fracción del PCI que sale como trabajo de eje (0 = sin elegir equipo)
#   φ = fracción del agua que condensa a T_salida (el resto sale con el N2)
# Calor de UNA liberación (fijo, no depende de la maniobra): Q = (m/L(φ))·q_rech.
# Lo que la maniobra decide es en cuánto tiempo hay que rechazarlo:
#   - de golpe (quema centrada, EXP-02):   P = ṁ·q_rech durante T
#   - con una capacidad P_cap:             t_min = Q / P_cap
#   - una sola liberación limitada por P:  H = g m³ q² / (24 M_e L² P²)
#   - liberación fraccionada en n pasos iguales con la misma bajada total H:
#         ṁ_n = ṁ_1/n,  P_n = P_1/n,  duración total n·T_1  (mismo Q)
# Referencia de capacidad (orden de magnitud, no diseño): el casco completo
# como intercambiador pasivo, P_casco = A_total·h·ΔT.
#
#
# Palancas registradas (trabajo futuro): presión (sube el rocío y φ), membranas
# selectivas al agua y separación de gotas. Límite común: el N2 que sale no puede
# llevarse el latente (ΔT = 9 L_v / (m_N2 cp_N2)); siempre hace falta un sumidero.
#
# Se evalúan fórmulas (instantáneo): no es una simulación.
# =============================================================================

import io
import sys
import csv
import math
import parametros as p
import rutas
import minimax as mm
from exp02_quema_centrada import L_ef, quema_centrada, cargas_eval, H_LISTA

PASO_FRACCIONADO_KG = 500.0      # ejemplo del autor: soltar "de a poco"
P_LISTA_ATM = (1.0, 3.0, 10.0)   # presión total del escape, en atmósferas (palanca de presión)
AIRE_ESTEQ = p.ratio_O2_H2 / p.w_O2_aire          # kg de aire por kg de H2 (≈ 34.3)
N_INERTE = (AIRE_ESTEQ - p.ratio_O2_H2) / p.M_inerte_aire * p.M_H2   # mol N2+Ar por mol H2O


_IF97 = (0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
         -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
         -0.23855557567849, 0.65017534844798e3)


def p_sat(T_C):
    """Presión de saturación del agua [Pa], IAPWS-IF97 (ec. 30), válida de 0 a 374 °C."""
    n1, n2, n3, n4, n5, n6, n7, n8, n9, n10 = _IF97
    T = T_C + 273.15
    th = T + n9 / (T - n10)
    A, B, C = th * th + n1 * th + n2, n3 * th * th + n4 * th + n5, n6 * th * th + n7 * th + n8
    return (2.0 * C / (-B + math.sqrt(B * B - 4.0 * A * C))) ** 4 * 1e6


# Densidad del agua líquida saturada, IAPWS: Wagner & Pruss, J. Phys. Chem. Ref. Data 22 (1993) 783, ec. 2.6
# (la de la release de propiedades de saturación); con la IAPWS-95 del NIST WebBook coincide en < 10 ppm de 0 a
# 100 °C (verificación, grupo L). La usa la caja del tanque de agua (balance_masa)
_SAT_RHO = (1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5)
_SAT_EXP = (1.0 / 3.0, 2.0 / 3.0, 5.0 / 3.0, 16.0 / 3.0, 43.0 / 3.0, 110.0 / 3.0)


def rho_liq_sat(T_C):
    """Densidad del agua líquida saturada [kg/m3] a T_C [°C] (IAPWS, Wagner & Pruss 1993)."""
    tau = 1.0 - (T_C + 273.15) / 647.096
    return 322.0 * (1.0 + sum(b * tau ** e for b, e in zip(_SAT_RHO, _SAT_EXP)))


def phi_max(T_C, P=p.P_atmos):
    """Fracción del agua que condensa si la mezcla sale saturada a T_C."""
    ps = p_sat(T_C)
    return max(0.0, min(1.0, 1.0 - N_INERTE * ps / (P - ps)))


def T_rocio(P=p.P_atmos):
    x = 1.0 / (1.0 + N_INERTE)                     # fracción molar de H2O en el escape
    lo, hi = 0.0, 300.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if p_sat(mid) < x * P else (lo, mid)
    return 0.5 * (lo + hi), x


def w_comp_escape(P, T_C, P0=p.P_atmos):
    """Trabajo MÍNIMO [J por kg de H2] para comprimir el escape de P0 a P, isotérmico y reversible a T_C.
    Solo se comprime la fase gaseosa (hallazgo M-1, auditoría 3; antes (1 + N) R T ln(P/P0), toda el agua como
    gas, 1.5 veces más a 14 atm y 40 °C): mientras el escape no está saturado van 1 + N mol de gas por mol de H2O;
    desde la presión de rocío P_r = p_s (1 + N) el vapor queda en p_s y lo que condensa sale líquido, con
    n_gas = N P / (P - p_s). Integrando n_gas R T dP / P:
        W = (1 + N) R T ln(P_r / P0) [si P_r > P0]  +  N R T ln((P - p_s) / (max(P0, P_r) - p_s))."""
    if P <= P0:
        return 0.0
    ps = p_sat(T_C)
    RT = p.R_univ * (T_C + 273.15) / p.M_H2          # por kg de H2 (1 mol de H2O por mol de H2)
    P_r = ps * (1.0 + N_INERTE)
    if P <= P_r:                                     # no llega a saturar: todo gas
        return (1.0 + N_INERTE) * RT * math.log(P / P0)
    seco = (1.0 + N_INERTE) * RT * math.log(P_r / P0) if P_r > P0 else 0.0
    return seco + N_INERTE * RT * math.log((P - ps) / (max(P0, P_r) - ps))


def q_rech(eta, phi, lam=1.0):
    """Calor a rechazar por kg de H2 [J/kg]: balance de energía con referencia a T_ref_comb (25 °C).
    Entran 1 kg de H2 y lam x aire estequiométrico a T_aire_exterior; sale el escape a T_salida_cond
    con la fracción phi del agua condensada (líquida) y el resto como vapor. eta = fracción del PCI
    que sale como trabajo. ÚNICA fórmula de q del modelo (la usan EXP-05 a EXP-11)."""
    r_agua = 1.0 + p.ratio_O2_H2
    aire = lam * AIRE_ESTEQ
    entra = (p.cp_H2 + aire * p.cp_aire) * (p.T_aire_exterior - p.T_ref_comb)
    sale = ((aire - p.ratio_O2_H2) * p.cp_N2
            + r_agua * (phi * p.cp_agua_liq + (1.0 - phi) * p.cp_vapor)) * (p.T_salida_cond - p.T_ref_comb)
    return (1.0 - eta) * p.Q_comb_H2 + phi * r_agua * p.L_v_agua + entra - sale


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    k0 = p.k_masa_anadida
    Tdp, x_h2o = T_rocio()
    phi = phi_max(p.T_salida_cond)
    L = L_ef(phi)
    L0 = L_ef(0.0)
    q = q_rech(0.0, phi)                 # sin elegir equipo (todo el calor va al escape)
    q_gen = q_rech(p.rend_gen, phi)      # si una fracción rend_gen sale como trabajo de eje
    q_lat = (1.0 + p.ratio_O2_H2) * p.L_v_agua   # lo que cuenta hoy el artículo (solo latente, φ = 1)
    c_art = rutas.leer_datos('condensacion')
    dl = rutas.leer_datos('dlim_flujo')

    filas = []
    for n, d in p.PLATAFORMAS.items():
        Pc_min, Pc_max = (d['A_total'] * h * p.dT_casco for h in p.h_conv_casco)
        c = 0.5 * p.rho_aire_std * p.Cd * d['Area']
        k0 = d['k']                                          # masa añadida de esta nave
        for m_t in cargas_eval(n):
            m = m_t * 1000.0
            M = p.rho_aire_std * d['V_casco']                # aire del casco: inercia (1 + k)·M (inercia.py)
            Me = (1.0 + k0) * M
            Q = m / L * q
            fila = {'plataforma': n, 'corto': d['corto'], 'carga_t': m_t, 'M_total_kg': M,
                    'H2_kg': m / L, 'Q_GJ': Q / 1e9, 'P_casco_min_MW': Pc_min / 1e6, 'P_casco_max_MW': Pc_max / 1e6,
                    't_min_casco_max_min': Q / Pc_max / 60.0, 't_min_casco_min_min': Q / Pc_min / 60.0,
                    'P_10min_MW': Q / (p.t_cond_objetivo * 60.0) / 1e6,
                    'H_una_liberacion_casco_max_m': p.g * m**3 * q**2 / (24.0 * Me * L**2 * Pc_max**2)}
            for H in H_LISTA:
                r = quema_centrada(M, m, H, k0, phi, c)
                fila[f'P_rt_H{H:g}_MW'] = r['dm_kgps'] * q / 1e6
                fila[f'T_H{H:g}_s'] = r['T_s']
            # liberación fraccionada, misma bajada total H = h_max_seguridad (masa añadida, sin arrastre)
            H5 = p.h_max_seguridad
            r1 = quema_centrada(M, m, H5, k0, phi)
            P1 = r1['dm_kgps'] * q
            n_min = max(1, math.ceil(P1 / Pc_max))
            n_500 = max(1, round(m / PASO_FRACCIONADO_KG))
            fila.update({'P1_H5_MW': P1 / 1e6, 'T1_H5_s': r1['T_s'], 'n_min_casco_max': n_min,
                         'paso_min_kg': m / n_min, 'duracion_n_min_min': n_min * r1['T_s'] / 60.0,
                         'n_500kg': n_500, 'P_500kg_MW': P1 / n_500 / 1e6, 'dm_500kg_kgps': r1['dm_kgps'] / n_500,
                         'duracion_500kg_min': n_500 * r1['T_s'] / 60.0})
            filas.append(fila)

    with open(rutas.dato('exp05_condensacion_tiempo_real', 'csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        for r in filas:
            w.writerow({kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in r.items()})

    def buscar(n, m_t):
        return next(f for f in filas if f['plataforma'] == n and f['carga_t'] == m_t)
    pf, fw = buscar(p.nombre_dirigible, 1.0), buscar('Flying Whales LCA60T', 20.0)
    dm_art = dl['plataformas'][p.nombre_dirigible]['dm_star_kgps']
    rutas.guardar_datos('exp05_condensacion_tiempo_real', {
        'T_rocio_C': Tdp, 'x_H2O': x_h2o, 'T_salida_C': p.T_salida_cond, 'phi_max_salida': phi,
        'phi_max_T1': phi_max(p.T1 - 273.15), 'L_phi': L, 'L_phi0': L0, 'aire_esteq': AIRE_ESTEQ,
        'q_rech_MJkg': q / 1e6, 'q_rech_gen_MJkg': q_gen / 1e6, 'q_latente_MJkg': q_lat / 1e6,
        'eta_gen': p.rend_gen, 'factor_vs_latente': q / q_lat,
        'articulo_Qcond_min_kW': c_art['Q_cond_min_kW'],
        'articulo_base_dm_star_kgps': dm_art, 'articulo_base_P_rt_MW': dm_art * q / 1e6,
        'presion': [{'P_atm': P, 'T_rocio_C': T_rocio(P * p.P_atmos)[0],
                     'phi_40C': phi_max(p.T_salida_cond, P * p.P_atmos),
                     'W_iso_MJkg': (1.0 + N_INERTE) / p.M_H2 * p.R_univ * (p.T_salida_cond + 273.15)
                                   * math.log(P) / 1e6} for P in P_LISTA_ATM],
        'dT_N2_latente_K': (1.0 + p.ratio_O2_H2) * p.L_v_agua / ((AIRE_ESTEQ - p.ratio_O2_H2) * p.cp_N2),
        'base_1t': pf, 'FW_20t': fw, 'paso_fraccionado_kg': PASO_FRACCIONADO_KG,
    })

    print(f'Escape estequiométrico: {AIRE_ESTEQ:.1f} kg aire/kg H2, x_H2O = {x_h2o:.3f}, rocío {Tdp:.1f} °C')
    print(f'A {p.T_salida_cond:g} °C condensa φ = {phi:.3f} (a {p.T1 - 273.15:g} °C: {phi_max(p.T1 - 273.15):.3f});'
          f'  L(φ) = {L:.2f}, L0 = {L0:.2f}')
    print(f'Calor a rechazar: {q/1e6:.1f} MJ/kg H2 (η = 0), {q_gen/1e6:.1f} (η = {p.rend_gen:g});'
          f' solo latente (artículo): {q_lat/1e6:.1f}  -> x{q/q_lat:.1f}')
    print(f'Artículo, caso base 1 t: ṁ* = {dm_art:.2f} kg/s -> {dm_art*q/1e6:.0f} MW en tiempo real'
          f' (artículo: Q_cond,min = {c_art["Q_cond_min_kW"]:.0f} kW en 10 min, solo latente)\n')
    print(f'{"plataforma":<22}{"carga":>6}{"Q GJ":>7}{"casco MW":>12}{"t_min min":>13}{"P10min MW":>10}'
          f'{"P_rt 5m":>9}{"H 1 lib.":>10}{"n_min":>7}{"n500":>6}{"P500 MW":>9}')
    for f in filas:
        print(f'{f["plataforma"]:<22}{f["carga_t"]:>6g}{f["Q_GJ"]:>7.1f}'
              f'{f["P_casco_min_MW"]:>6.1f}-{f["P_casco_max_MW"]:<5.1f}'
              f'{f["t_min_casco_max_min"]:>6.0f}-{f["t_min_casco_min_min"]:<6.0f}{f["P_10min_MW"]:>10.1f}'
              f'{f["P_rt_H5_MW"]:>9.0f}{f["H_una_liberacion_casco_max_m"]:>10.0f}{f["n_min_casco_max"]:>7d}'
              f'{f["n_500kg"]:>6d}{f["P_500kg_MW"]:>9.1f}')
    graficar(filas, q)


def graficar(filas, q):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import plot_style  # noqa
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.8))
    P = np.logspace(0, 4, 200)                        # MW
    for n, d in p.PLATAFORMAS.items():
        cap = cargas_eval(n)[-1]
        f = next(f for f in filas if f['plataforma'] == n and f['carga_t'] == cap)
        a1.plot(H_LISTA, [f[f'P_rt_H{H:g}_MW'] for H in H_LISTA], color=d['color'], lw=2.2, marker='o', ms=4,
                label=f'{n} ({cap:g} t)')
        a1.axhline(f['P_casco_max_MW'], color=d['color'], ls=':', lw=1.4)
        a2.plot(P, f['Q_GJ'] * 1e3 / P / 60.0, color=d['color'], lw=2.2, label=f'{n} ({cap:g} t)')
        a2.plot([f['P_casco_min_MW'], f['P_casco_max_MW']],
                [f['t_min_casco_min_min'], f['t_min_casco_max_min']], color=d['color'], lw=7,
                alpha=0.45, solid_capstyle='round')
    x = rutas.leer_datos('exp05_condensacion_tiempo_real')
    a1.plot([p.h_max_precision], [x['articulo_base_P_rt_MW']], marker='*', ms=13, color='black', ls='none',
            label=f'Article design point ({p.nombre_dirigible}, 1 t, minimax)')
    a1.set_title('(a) Real-time heat rejection, single release (centered burn)', fontsize=10)
    a1.set_xscale('log'); a1.set_yscale('log')
    a1.set_xticks(H_LISTA); a1.set_xticklabels([f'{h:g}' for h in H_LISTA])
    a1.set_xlabel(r'Allowed total descent  $H$  [m]')
    a1.set_ylabel('Heat to reject during the burn  [MW]')
    a2.axhline(p.t_cond_objetivo, color='black', ls='--', lw=1.1)
    a2.text(1.2, p.t_cond_objetivo * 1.12, f'{p.t_cond_objetivo:g} min (article window)', fontsize=8)
    a2.set_title('(b) Minimum time to reject the heat of one release', fontsize=10)
    a2.set_xscale('log'); a2.set_yscale('log')
    a2.set_xlabel('Installed heat-rejection power  [MW]')
    a2.set_ylabel('Time  [min]')
    for ax in (a1, a2):
        ax.grid(True, which='both', alpha=0.2)
    # una sola leyenda abajo: los colores de plataforma son los mismos en ambos paneles
    h, l = a1.get_legend_handles_labels()
    h += [Line2D([], [], color='gray', ls=':', lw=1.4), Line2D([], [], color='gray', lw=7, alpha=0.45)]
    l += ['(a) Whole hull as passive exchanger, forced air (upper ref.)',
          '(b) Whole-hull range, calm to forced air']
    fig.legend(h, l, loc='lower center', ncol=3, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(rutas.figura('exp05_condensacion_tiempo_real'), dpi=200, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
