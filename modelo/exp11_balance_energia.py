# -*- coding: utf-8 -*-
# [EN] Energy balance of the black boxes and self-powered cooling (Section 5.1, Appendix C).
# =============================================================================
# exp11_balance_energia.py — EXP-11: el ciclo con cajas negras (paradigma del
# autor: grupo electrógeno + baterías que alimentan el enfriamiento).
#
# Cajas negras (ninguna es un equipo elegido):
#   unidad de conversión de energía: H2 + aire -> trabajo eléctrico eta·PCI + escape
#   buffer de energía:              electricidad <-> electricidad (Wh/kg, kW/kg)
#   unidad de rechazo de calor:     escape + aire ambiente -> agua líquida + N2;
#                                   mover el aire pide la fracción f del calor rechazado
#   unidad de regeneración:         agua + sol -> H2
#
# A. Autoabasto: eta·PCI >= f·Q + W_comp, con Q = (1-eta) PCI + 9 phi L_v + W_comp
#    (la compresión isotérmica también se rechaza como calor). Casos: 1 atm (phi de
#    equilibrio a 40 °C) y presurizado a la presión de EXP-06 (phi = 0.99).
# B. Buffer: deja el condensador andando n_tau·tau antes de la quema (EXP-10 pide
#    tiempo de respuesta corto) y, como cota conservadora, mueve solo el aire hasta
#    el fin de la quema. ¿Lo recarga la quema del mismo escalón?
# C. Regeneración con la retención real (EXP-10): el agua venteada es H2 perdido.
# =============================================================================

import numpy as np
import parametros as p
import rutas
from exp05_condensacion_tiempo_real import phi_max, q_rech

ETA_BARRIDO = np.linspace(0.0, 0.45, 91)
CASOS_REGEN = ('ideal', 'ventear_tau0s', 'ventear_tau5s')


def q_total(eta, phi, W_comp):
    """Calor a rechazar por kg de H2 [J/kg]: balance de EXP-05 + compresión (que también se rechaza)."""
    return q_rech(eta, phi) + W_comp


def eta_minima(phi, f, W_comp):
    """Menor eta con la que el trabajo alcanza para mover el aire y comprimir."""
    return (f * (q_rech(0.0, phi) + W_comp) + W_comp) / ((1.0 + f) * p.Q_comb_H2)


def razon_autoabasto(eta, phi, f, W_comp):
    """Trabajo generado / trabajo que pide el enfriamiento (por kg de H2)."""
    return eta * p.Q_comb_H2 / (f * q_total(eta, phi, W_comp) + W_comp)


def main():
    rho_H2 = p.P_atmos / (p.R_H2 * p.T1)
    k_H2 = p.rho_aire_std / rho_H2 - 1.0                       # rho_aire/rho_H2 - 1
    e6 = rutas.leer_datos('exp06_requisitos_minimos')
    e10 = rutas.leer_datos('exp10_escalones_condensador')
    reg = rutas.leer_datos('fuente_O2')['regeneracion_solar']

    # ── A: autoabasto ──────────────────────────────────────────────────────
    presiones = {'1atm': (phi_max(p.T_salida_cond), 0.0, 1.0),
                 'presurizado': (e6['phi_objetivo'], e6['W_comp_MJkg'] * 1e6, e6['P_req_atm'])}
    autoabasto = []
    for clave, (phi, W, P_atm) in presiones.items():
        B_t = 1000.0 / (k_H2 + (1.0 + p.ratio_O2_H2) * phi)                       # H2 por tonelada soltada [kg/t]
        for f, ext in zip(p.perdida_rad, ('conservador', 'optimista')):
            eta_min = eta_minima(phi, f, W)
            Q_t = B_t * q_total(p.rend_gen, phi, W)
            E_gen = B_t * p.rend_gen * p.Q_comb_H2
            E_pide = f * Q_t + B_t * W
            autoabasto.append({'caso': clave, 'extremo': ext, 'phi': phi, 'P_atm': P_atm, 'f': f,
                               'W_comp_MJkg': W / 1e6, 'eta_min': eta_min, 'eta_ref': p.rend_gen,
                               'H2_por_t_kg': B_t, 'E_gen_kWh_t': E_gen / 3.6e6, 'E_pide_kWh_t': E_pide / 3.6e6,
                               'razon': E_gen / E_pide,
                               'razon_barrido': [razon_autoabasto(e, phi, f, W) for e in ETA_BARRIDO]})

    # ── B: buffer que arranca el condensador antes de soltar ───────────────
    tau = p.tau_condensador
    buffer = []
    # rechazo: masa por kW del sistema de gestión térmica completo (FlyZero p. 4; caja negra, hallazgos I5/S1
    # de la auditoría, 2026-09-27; antes dos bandas de EXP-14), 2026 = conservador, 2035 = optimista
    kgkw_rad = tuple(1.0 / k for k in p.kWkg_TMS_sistema)
    for x in e10['casos']:
        if x['caso'] != 'ventear_tau5s':
            continue
        P_cap = x['P_instalada_MW'] * 1e6
        t_buf = p.n_tau_arranque * tau + x['dur_quema_s']
        E_esc = p.rend_gen * p.Q_comb_H2 * x['H2_quemado_kg']  # lo que genera la quema del escalón
        carga_max = p.PLATAFORMAS[x['plataforma']]['carga_max']
        for f, kWkg, kgkw, ext in zip(p.perdida_rad, p.kWkg_bateria, kgkw_rad, ('conservador', 'optimista')):
            m_rad = P_cap / 1e3 * kgkw
            P_buf = f * P_cap
            E_buf = P_buf * t_buf
            m_E = E_buf / 3600.0 / p.Whkg_bateria
            m_P = P_buf / 1e3 / kWkg
            buffer.append({'plataforma': x['plataforma'], 'corto': x['corto'], 'P_instalada_MW': x['P_instalada_MW'],
                           'extremo': ext, 'escalon_kg': x['escalon_kg'], 'P_buf_MW': P_buf / 1e6, 't_buf_s': t_buf,
                           'E_buf_kWh': E_buf / 3.6e6, 'tasa_C': 3600.0 / t_buf,   # la de un paquete por ENERGÍA
                           # la del paquete que se usa (el mayor de energía y potencia): P / (m Wh/kg) [1/h]
                           'tasa_C_paquete': P_buf / (max(m_E, m_P) * p.Whkg_bateria),
                           'm_bat_energia_kg': m_E, 'm_bat_potencia_kg': m_P, 'm_bat_kg': max(m_E, m_P),
                           'manda': 'potencia' if m_P >= m_E else 'energía',
                           'm_bat_por_escalon': max(m_E, m_P) / x['escalon_kg'],
                           'E_escalon_kWh': E_esc / 3.6e6, 'recarga': E_esc / E_buf,
                           'm_rad_kg': m_rad, 'kg_por_kW': (max(m_E, m_P) + m_rad) / (P_cap / 1e3),
                           'frac_carga': (max(m_E, m_P) + m_rad) / carga_max})

    # ── C: regeneración con la retención real ──────────────────────────────
    B_ideal = 1000.0 / (k_H2 + 1.0 + p.ratio_O2_H2)
    regen = []
    for n, d in p.PLATAFORMAS.items():
        kg_dia = reg[n]['regen_kg_dia']
        regen.append({'plataforma': n, 'corto': d['corto'], 'P_instalada_MW': None, 'caso': 'forma_cerrada',
                      'H2_por_t_kg': B_ideal, 'retenida': 1.0, 'recuperado_kg_t': B_ideal, 'perdido_kg_t': 0.0,
                      'dias_por_t': B_ideal / kg_dia, 'regen_kg_dia': kg_dia})
        for x in e10['casos']:
            if x['plataforma'] != n or x['caso'] not in CASOS_REGEN:
                continue
            r = min(x['frac_retenida'], 1.0)
            B = x['H2_por_t_kg']
            regen.append({'plataforma': n, 'corto': d['corto'], 'P_instalada_MW': x['P_instalada_MW'],
                          'caso': x['caso'], 'H2_por_t_kg': B, 'retenida': r, 'recuperado_kg_t': r * B,
                          'perdido_kg_t': (1.0 - r) * B, 'dias_por_t': r * B / kg_dia, 'regen_kg_dia': kg_dia})

    # ── D: ¿otro sumidero que la atmósfera? Nitrógeno líquido que se evapora y se ventea ──
    # NIST (SRD 69): N2 a 1 atm, T_sat = 77.355 K, h_fg = 199.18 kJ/kg; el gas sale a la temperatura exterior.
    phi_obj = p.phi_objetivo
    B_obj = 1000.0 / (k_H2 + (1.0 + p.ratio_O2_H2) * phi_obj)          # H2 por tonelada, condensando todo
    absorbe_LN2 = 199.18e3 + p.cp_N2 * (p.T_aire_exterior + 273.15 - 77.355)   # J por kg de LN2
    sumidero = {'H2_por_t_kg': B_obj, 'Q_por_t_GJ': B_obj * q_rech(0.0, phi_obj) / 1e9,
                'Q_latente_por_t_GJ': B_obj * (1.0 + p.ratio_O2_H2) * phi_obj * p.L_v_agua / 1e9,
                'absorbe_LN2_kJkg': absorbe_LN2 / 1e3}
    sumidero['LN2_t_por_t'] = sumidero['Q_por_t_GJ'] * 1e9 / absorbe_LN2 / 1000.0
    sumidero['LN2_latente_t_por_t'] = sumidero['Q_latente_por_t_GJ'] * 1e9 / absorbe_LN2 / 1000.0

    rutas.guardar_datos('exp11_balance_energia', {
        'eta_ref': p.rend_gen, 'eta_barrido': list(ETA_BARRIDO), 'tau_s': tau, 'n_tau_arranque': p.n_tau_arranque,
        'H2_por_t_venteo_kg': e10['H2_por_t_venteo_kg'], 'H2_por_t_ideal_kg': B_ideal,
        'autoabasto': autoabasto, 'buffer': buffer, 'regeneracion': regen, 'sumidero_LN2': sumidero})
    graficar(autoabasto, regen, e10['H2_por_t_venteo_kg'])
    tabla_buffer(buffer)

    print('A. Autoabasto (eta_ref = %.2f)' % p.rend_gen)
    for a in autoabasto:
        print(f"  {a['caso']:<12}{a['extremo']:<12} phi {a['phi']:.2f}  eta_min {a['eta_min']:.3f}  "
              f"genera {a['E_gen_kWh_t']:5.0f} kWh/t  pide {a['E_pide_kWh_t']:5.0f} kWh/t  razón {a['razon']:.2f}")
    print('B. Buffer (tau = %g s, arranque %g tau)' % (tau, p.n_tau_arranque))
    for b in buffer:
        print(f"  {b['corto']:<5}{b['P_instalada_MW']:>5.0f} MW {b['extremo']:<12} P {b['P_buf_MW']:5.1f} MW  "
              f"t {b['t_buf_s']:5.1f} s  E {b['E_buf_kWh']:6.1f} kWh  {b['tasa_C']:4.0f} C  "
              f"masa {b['m_bat_kg']:6.0f} kg ({b['manda']})  recarga {b['recarga']:.2f}  "
              f"rad {b['m_rad_kg']:6.0f} kg  {b['kg_por_kW']:.2f} kg/kW  {100*b['frac_carga']:5.1f} % carga")
    print('C. Regeneración por tonelada')
    for g in regen:
        print(f"  {g['corto']:<5}{str(g['P_instalada_MW']):>6} {g['caso']:<15} quemado {g['H2_por_t_kg']:5.1f}  "
              f"recuperado {g['recuperado_kg_t']:5.1f}  perdido {g['perdido_kg_t']:5.1f} kg/t  "
              f"{g['dias_por_t']:.2f} días/t")


def tabla_buffer(buffer):
    lin = [r'\begin{tabular}{lrrrrrrrr}', r'\toprule',
           r'Plataforma & $P_\text{cap}$ [MW] & $P_\text{buf}$ [MW] & $t_\text{buf}$ [s] & $E_\text{buf}$ [kWh]'
           r' & Tasa [C] & Buffer [t] & Recarga [--] & Rad.\,+\,buffer [\% carga] \\', r'\midrule']
    claves = list(dict.fromkeys((b['corto'], b['P_instalada_MW']) for b in buffer))
    for corto, P in claves:
        c, o = [b for b in buffer if b['corto'] == corto and b['P_instalada_MW'] == P]
        lin.append(f"{corto} & {P:.0f} & {o['P_buf_MW']:g}--{c['P_buf_MW']:g} & {c['t_buf_s']:.0f} & "
                   f"{o['E_buf_kWh']:.0f}--{c['E_buf_kWh']:.0f} & {c['tasa_C']:.0f} & "
                   f"{o['m_bat_kg'] / 1e3:.1f}--{c['m_bat_kg'] / 1e3:.1f} & {c['recarga']:.1f}--{o['recarga']:.1f} & "
                   f"{100 * o['frac_carga']:.0f}--{100 * c['frac_carga']:.0f}" + r' \\')
    lin += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('exp11_buffer'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lin) + '\n')


def graficar(autoabasto, regen, kg_t_venteo):
    import matplotlib.pyplot as plt
    import plot_style  # noqa
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    ax = axes[0]
    estilos = {('1atm', 'conservador'): ('#0072B2', '-'), ('1atm', 'optimista'): ('#0072B2', '--'),
               ('presurizado', 'conservador'): ('#D55E00', '-'), ('presurizado', 'optimista'): ('#D55E00', '--')}
    for a in autoabasto:
        col, ls = estilos[(a['caso'], a['extremo'])]
        etiqueta = f"{a['P_atm']:.0f} atm, $\\phi$ = {a['phi']:.2f}, f = {a['f']:.2f}"
        ax.plot(ETA_BARRIDO, a['razon_barrido'], color=col, ls=ls, lw=1.8, label=etiqueta)
    ax.axhline(1.0, color='gray', lw=0.8)
    ax.axvline(p.rend_gen, color='black', ls=':', lw=1)
    ax.text(p.rend_gen, 0.2, f'  $\\eta$ = {p.rend_gen:.2f} (article)', fontsize=8)
    ax.set_xlabel(r'Work fraction of the LHV exported by the conversion unit, $\eta$')
    ax.set_ylabel('Work generated / work needed to cool')
    ax.set_ylim(0, 4)
    ax.set_xlim(0, ETA_BARRIDO[-1])
    ax.set_title('Self-powered heat rejection (black boxes)', fontsize=10)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.2)

    ax = axes[1]
    cortos = list(dict.fromkeys(g['corto'] for g in regen))
    orden = [('forma_cerrada', None), ('ventear_tau0s', 100.0), ('ventear_tau5s', 100.0), ('ventear_tau5s', 10.0)]
    nombres = {'forma_cerrada': 'Ideal (closed form)', 'ventear_tau0s': 'Lag 0 s', 'ventear_tau5s': 'Lag 5 s'}
    ancho = 0.8 / len(orden)
    for i, corto in enumerate(cortos):
        for k, (caso, P) in enumerate(orden):
            g = next(y for y in regen if y['corto'] == corto and y['caso'] == caso and y['P_instalada_MW'] == P)
            xk = i - 0.4 + (k + 0.5) * ancho
            ax.bar(xk, g['recuperado_kg_t'], ancho * 0.9, color='#0072B2',
                   label='Recovered by electrolysis' if (i, k) == (0, 0) else None)
            ax.bar(xk, g['perdido_kg_t'], ancho * 0.9, bottom=g['recuperado_kg_t'], color='#D55E00',
                   label='Lost with the vented vapour' if (i, k) == (0, 1) else None)
    ax.axhline(kg_t_venteo, color='gray', ls='--', lw=1)
    ax.text(-0.45, kg_t_venteo, 'venting H$_2$ (all lost)', va='bottom', ha='left', fontsize=8, color='gray')
    ax.set_xticks(range(len(cortos)))
    ax.set_xticklabels(cortos)
    ax.set_ylabel('H$_2$ per tonne released [kg/t]')
    ax.set_title('Per platform: ideal, lag 0 s (100 MW), lag 5 s (100 MW), lag 5 s (10 MW)', fontsize=9)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, axis='y', alpha=0.2)
    fig.suptitle('The cycle with black boxes: energy balance and regeneration with real water retention',
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(rutas.figura('exp11_balance_energia'), dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
