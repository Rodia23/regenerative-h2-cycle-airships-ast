# -*- coding: utf-8 -*-
# [EN] Bounds and sensitivities of the phenomena ranking (Supplementary Table S2, Section 9).
# =============================================================================
# val03_cotas.py — VAL-03: cotas y sensibilidades de los fenómenos de la PIRT que no
# tienen validación directa (PIRT: tabla S2 del suplementario).
#
# A. P8, retardo del condensador tau. Un radiador de masa concentrada tiene
#    tau = C/UA = c/k por kg, con k [W/(kg K)] la calibración de FlyZero de EXP-14
#    (kW/kg = k (T - T_aire)) y c el calor específico del aluminio (JANAF). Cota alta:
#    la constante de tiempo medida de un intercambiador enfriado por aire de
#    laboratorio (Ghatos et al. 2023, registro GH23). Se comparan con el tau supuesto.
# B. P5, arrastre vertical. El caso base del artículo y uno con condensador real
#    (100 MW, tau = 5 s, ventear) con Cd = 0 (sin arrastre), 0.5 (modelo), 0.8 (JPL,
#    registro IZ25) y 1.0.
# C. P6, gas isotérmico. Tiempo térmico del H2 a escala real: capacidad del gas sobre
#    h_int A, con h_int de la correlación de convección natural de Churchill & Chu (1975)
#    y propiedades del H2 de NIST, casco elipsoide de esbeltez 4 (la misma de k). Un
#    cambio de sobrecalentamiento de equilibrio dT_eq (5-20 K) durante la maniobra da una
#    deriva de flotabilidad M_total dT(t)/T; se compara con la carga y el escalón.
# D. P3 y P4, propagación de incertidumbre con las formas cerradas del minimax:
#    masa vertical TOTAL (propia + añadida) medida/analítica 0.77-1.20 (Yamasaki & Goto,
#    MLG(2,2)), aplicada a toda la inercia como en VAL-01, y agua retenida
#    la recolección MEDIDA / cota de los experimentos de Zakhvatkin (VAL-02, ~0.68-0.85). Se reportan el MÍNIMO y
#    el MÁXIMO en los extremos medidos (esquinas; las cantidades son monótonas en cada entrada), no percentiles:
#    con dos ensayos y seis experimentos no hay distribución que muestrear (hallazgo S4, autor 2026-09-29;
#    antes Monte Carlo uniforme, p5-p95).
# E. (revisión 2026-10-01, R-07) Humedad del aire de entrada. El modelo usa aire seco. El aire de combustión húmedo
#    trae x mol de vapor por mol de agua producida; como el escape sale saturado a T_salida, el vapor que se va
#    (1 - phi_max por mol de agua producida) no depende de lo que entra: todo el vapor que entra condensa. Suma
#    lastre (r x por kg de H2) y calor (su latente, más su calor sensible de entrada y salida), y el aire húmedo es un
#    poco menos denso (rho_aire/rho_H2 por 1 - (1 - M_H2O/M_aire) p_v/P). A 1 atm y phi_max.
# Salidas: resultados/datos/val03_cotas.json, resultados/tablas/val03_*.tex,
#          resultados/figuras/val03_cotas.png
# =============================================================================

import math
import os
import numpy as np
import parametros as p
import rutas
import fuentes
import minimax as mm
import condensador
import drop_flujo_optimo_perdrop as sim
from exp02_quema_centrada import L_ef
from exp05_condensacion_tiempo_real import q_rech, phi_max, p_sat, AIRE_ESTEQ

AQUI = os.path.dirname(os.path.abspath(__file__))
CD_LISTA = (0.0, 0.5, 0.8, 1.0)
DT_EQ = (5.0, 10.0, 20.0)                # [K] cambio del sobrecalentamiento de equilibrio (rango declarado)
H_FACTOR = (0.5, 2.0)                    # incertidumbre declarada de h_int (recinto, no placa plana)
ESBELTEZ = p.NAVE_BASE['esbeltez']       # la del caso base (su k de masa añadida sale de ella)
FK = (0.77, 1.20)                        # masa vertical total medida / analítica (YG03, MLG(2,2))


def cp_aluminio():
    """Cp del Al(cr) a 298.15 K de la tabla JANAF Al-001 [J/(kg K)]."""
    for l in open(os.path.join(AQUI, 'datos_fuente', 'janaf', 'Al-001.txt'), encoding='utf-8'):
        c = l.strip().split('\t')
        if c[0] == '298.15':
            return float(c[1]) / 26.9815385e-3          # M del Al (IUPAC)
    raise ValueError('Al-001 sin fila 298.15')


def props_H2(T=300.0):
    """rho, cp, mu, k del H2 a 1 atm y T [K] de NIST (isobara guardada en datos_fuente/nist)."""
    for l in open(os.path.join(AQUI, 'datos_fuente', 'nist', 'H2_isobara_0.101325MPa_280-320K.txt'), encoding='utf-8'):
        c = l.strip().split('\t')
        try:
            if abs(float(c[0]) - T) < 1e-6:
                return {'rho': float(c[2]), 'cp': float(c[8]) * 1e3, 'mu': float(c[11]), 'k': float(c[12])}
        except ValueError:
            continue
    raise ValueError(f'sin fila a {T} K')


def parte_A():
    c_al = cp_aluminio()
    k = [x * 1e3 for x in rutas.leer_datos('exp14_calor_por_temperatura')['k_fria_kWkg_por_K']]
    gh = fuentes.FUENTES['GH23']['datos']
    return {'c_Al_J_kgK': c_al, 'k_W_kgK': k, 'tau_FlyZero_s': [c_al / x for x in k],
            'tau_Ghatos_forzada_s': [60.0 * x for x in gh['tau_forzada_min'][0]],
            'tau_Ghatos_natural_s': [60.0 * x for x in gh['tau_natural_min'][0]],
            'tau_supuesto_s': p.tau_condensador, 'tau_cota_s': 30.0}


def corrida(Cd, con_cond):
    Cd0 = p.Cd
    p.Cd = Cd
    try:
        cond = condensador.config(100e6, tau_s=p.tau_condensador, modo='ventear') if con_cond else None
        r = sim.correr(graficar=False, verbose=False, cond=cond, devolver_serie=True)
    finally:
        p.Cd = Cd0
    s = r['serie']
    quemado = s['m_H2'][0] - s['m_H2']
    tot = float(np.max(quemado))
    t1, t99 = (float(s['t'][np.argmax(quemado >= f * tot)]) for f in (0.01, 0.99))
    return {'Cd': Cd, 'excursion_m': r['error_max_m'], 'H2_kg': tot,
            'agua_kg': float(s['agua'][-1]), 'v_max_ms': r['v_max_ms'], 'dur_quema_s': t99 - t1}


def parte_B():
    out = {}
    for caso, cc in (('base', False), ('condensador_100MW_tau5', True)):
        filas = [corrida(Cd, cc) for Cd in CD_LISTA]
        ref = next(f for f in filas if f['Cd'] == p.Cd)
        for f in filas:
            f['d_excursion_rel'] = f['excursion_m'] / ref['excursion_m'] - 1.0
            f['d_H2_rel'] = f['H2_kg'] / ref['H2_kg'] - 1.0
        out[caso] = filas
    return out


def termico(m_h2, pr):
    """Tiempo térmico del gas m_h2 [kg] en un casco elipsoide de esbeltez ESBELTEZ: capacidad del gas sobre h_int A."""
    rho_h2 = p.P_atmos / (p.R_H2 * p.T1)
    V = m_h2 / rho_h2
    b = (V / ((4.0 / 3.0) * math.pi * ESBELTEZ)) ** (1.0 / 3.0)      # semieje menor
    a = ESBELTEZ * b
    e = math.sqrt(1.0 - (b / a) ** 2)
    A = 2.0 * math.pi * b * b * (1.0 + a / (b * e) * math.asin(e))
    # Churchill & Chu (1975), placa vertical, todo Ra; L = diámetro (a Ra alto h no depende de L)
    nu = pr['mu'] / pr['rho']
    alpha = pr['k'] / (pr['rho'] * pr['cp'])
    Pr = nu / alpha
    dT_pared = 10.0
    L = 2.0 * b
    Ra = p.g * (1.0 / 300.0) * dT_pared * L ** 3 / (nu * alpha)
    Nu = (0.825 + 0.387 * Ra ** (1.0 / 6.0) / (1.0 + (0.492 / Pr) ** (9.0 / 16.0)) ** (8.0 / 27.0)) ** 2
    h = Nu * pr['k'] / L
    C_gas = m_h2 * p.cp_H2
    tau_th = [C_gas / (f * h * A) for f in (H_FACTOR[1], 1.0, H_FACTOR[0])]      # rápido, medio, lento
    return {'V_m3': V, 'A_m2': A, 'eje_a_m': a, 'eje_b_m': b, 'Pr': Pr, 'Ra': Ra, 'Nu': Nu, 'h_W_m2K': h,
            'C_gas_J_K': C_gas, 'tau_th_s': tau_th}


def parte_C(base):
    pr = props_H2(300.0)
    eq = mm.equilibrio(p.nombre_dirigible, p.m_carga)          # gas de equilibrio del caso base con 1 t
    m_h2, M_tot = eq['m_gas'], eq['M_total']
    T = termico(m_h2, pr)
    tau_th = T['tau_th_s']
    dur_minimax = base['dur_quema_s']                                           # del 1 al 99 % del H2 quemado
    # J por tonelada, con el phi objetivo como el resto del paper (antes phi = 1: 10.54 min en vez de 10.57,
    # y la discusión decía 10.5 mientras la sec. 5.3 daba 106 min para 10 t)
    Q_t = q_rech(0.0, p.phi_objetivo) * 1000.0 / L_ef(p.phi_objetivo)
    dur_pasos = Q_t / 10e6
    deriva = {}
    for dT in DT_EQ:
        for nom, t in (('minimax', dur_minimax), ('pasos_10MW_1t', dur_pasos)):
            dm = [M_tot * dT * (1.0 - math.exp(-t / tt)) / p.T1 for tt in tau_th]
            deriva[f'{dT:g}K_{nom}'] = {'dT_K': dT, 'duracion_s': t, 'dm_kg': dm}
    # La entrega de referencia de la sección 5.3 (EXP-22, 2026-10-13): sueltas una tras otra durante horas, con la
    # nave cargada al salir (equipo + carga de negocio) y su propio tiempo térmico. En horas el gas llega a la
    # temperatura nueva: la deriva es casi toda M dT/T1, del orden de varias sueltas.
    E = rutas.leer_datos('exp22_energia_potencia')['E']
    eq_e = mm.equilibrio(E['nave'], E['carga0_t'] * 1e3)
    T_e = termico(eq_e['m_gas'], pr)
    dur_entrega = E['n'] * E['t_paso_s']
    for dT in DT_EQ:
        dm = [eq_e['M_total'] * dT * (1.0 - math.exp(-dur_entrega / tt)) / p.T1 for tt in T_e['tau_th_s']]
        deriva[f'{dT:g}K_entrega'] = {'dT_K': dT, 'duracion_s': dur_entrega, 'dm_kg': dm}
    return {'m_H2_kg': m_h2, 'M_total_kg': M_tot, **T, 'props_H2': pr,
            'dur_minimax_s': dur_minimax, 'dur_pasos_s': dur_pasos, 'deriva': deriva,
            'entrega': {'M_total_kg': eq_e['M_total'], 'm_H2_kg': eq_e['m_gas'], 'carga_t': E['carga0_t'],
                        'tau_th_s': T_e['tau_th_s'], 'duracion_s': dur_entrega}}


def qois(fk, fphi, phi_cota):
    phi = phi_cota * fphi
    L = L_ef(phi)
    r = mm.resolver(p.NAVE_BASE, p.m_carga, p.h_max_precision, p.D_lim_oper, rho_gas=mm.rho_H2, L_gas=L,
                    k_ma=(1.0 + p.NAVE_BASE['k']) * fk - 1.0)     # M_I = (1 + k) rho V_casco por fk
                    # (antes k*fk: escalaba solo la parte añadida y el rango salía ~la mitad del de VAL-01)
    # la cota ideal combina la tasa con phi = 1 y el calor de condensar phi_objetivo (sec. 5.3 del artículo, la
    # misma regla de \valReqPcondBaseUno); antes q(phi = 1): base 447 MW frente a los 446 del texto (rev. 2026-10-01)
    q = q_rech(0.0, min(phi, p.phi_objetivo))
    return {'dm_star': r['dm_star'], 'P_MW': r['dm_star'] * q / 1e6, 'H2_t': 1000.0 / L,
            'Q_GJ_t': q * 1000.0 / L / 1e9}


def fphi_medido():
    """Recolección MEDIDA / cota de equilibrio en los experimentos de Zakhvatkin (VAL-02): mín. y máx.
    (antes FPHI = 0.68-1.0 copiado a mano: el 1.0 es la cota, no una medición; auditoría R2-04)."""
    r = [x['medido_sobre_modelo'] for x in rutas.leer_datos('val02_condensacion')['experimentos']]
    return min(r), max(r)


def parte_D():
    """Dos escenarios: (i) cota ideal, phi = 1 por definición: solo la masa añadida es incierta;
    (ii) condensación real: phi = phi_max(T_salida_cond) x eficiencia de recolección medida (VAL-02);
    su valor base, el centro del rango medido. Rango = mínimo y máximo en las esquinas de los extremos medidos."""
    fphi = fphi_medido()
    out = {'phi_max': phi_max(p.T_salida_cond), 'fphi_medido': fphi}
    for esc, cota, usa_phi in (('ideal', 1.0, False), ('real', out['phi_max'], True)):
        fb = 0.5 * (fphi[0] + fphi[1]) if usa_phi else 1.0
        base = qois(1.0, fb, cota)
        esquinas = [qois(a, b, cota) for a in FK for b in (fphi if usa_phi else (1.0,))]
        oat = {'masa_anadida': [qois(x, fb, cota) for x in FK]}
        if usa_phi:
            oat['agua_retenida'] = [qois(1.0, x, cota) for x in fphi]
        out[esc] = {'base': base, 'oat': oat,
                    'rango': {k: {'min': min(e[k] for e in esquinas), 'max': max(e[k] for e in esquinas)}
                              for k in base}}
    return out


def humedo(HR, eta=0.0, T_out=None, P=None):
    """E. Calor, pesadez y requisitos por tonelada con aire de entrada húmedo (HR = 0: aire seco del modelo)."""
    T_out = p.T_salida_cond if T_out is None else T_out
    P = p.P_atmos if P is None else P
    r = 1.0 + p.ratio_O2_H2
    pv = HR * p_sat(p.T_aire_exterior)
    n_aire = AIRE_ESTEQ / p.M_aire * p.M_H2              # mol de aire seco por mol de agua producida
    x = n_aire * pv / (P - pv)                           # mol de vapor que entra por mol de agua producida
    pm = phi_max(T_out, P)
    razon = (p.rho_aire_std / mm.rho_H2) * (1.0 - (1.0 - p.M_H2O / p.M_aire) * pv / P)
    L = razon - 1.0 + r * (pm + x)
    q = (q_rech(eta, pm) + r * x * (p.L_v_agua + p.cp_vapor * (p.T_aire_exterior - p.T_ref_comb)
                                    - p.cp_agua_liq * (T_out - p.T_ref_comb)))
    return {'HR': HR, 'x_mol': x, 'agua_extra_kg_por_kgH2': r * x, 'phi_max': pm, 'L': L, 'q_MJkg': q / 1e6,
            'H2_t': 1000.0 / L, 'Q_GJ_t': q * 1000.0 / L / 1e9}


def parte_E():
    seco = humedo(0.0)
    filas = [humedo(h) for h in p.HR_aire]
    for f in filas:
        f['dH2_t_rel'] = f['H2_t'] / seco['H2_t'] - 1.0
        f['dQ_t_rel'] = f['Q_GJ_t'] / seco['Q_GJ_t'] - 1.0
        f['dq_rel'] = f['q_MJkg'] / seco['q_MJkg'] - 1.0
    return {'seco': seco, 'humedo': filas}


def tablas(A, B, C, D):
    L = [r'\begin{tabular}{l r}', r'\toprule', r'Estimación de $\tau$ & [s] \\', r'\midrule',
         f"Radiador de FlyZero de aluminio concentrado, $c/k$ (2026 / 2035) & "
         f"{A['tau_FlyZero_s'][0]:.1f} / {A['tau_FlyZero_s'][1]:.1f} \\\\",
         f"Supuesto del modelo (EXP-09) & {A['tau_supuesto_s']:.0f} \\\\",
         f"Intercambiador de laboratorio, aire forzado 0.75 m/s (Ghatos 2023) & "
         f"{A['tau_Ghatos_forzada_s'][0]:.0f}--{A['tau_Ghatos_forzada_s'][1]:.0f} \\\\",
         f"Cota conservadora del modelo & {A['tau_cota_s']:.0f} \\\\", r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('val03_tau'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/val03_cotas.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')
    L = [r'\begin{tabular}{l ' + 'r ' * len(CD_LISTA) + '}', r'\toprule',
         '$C_d$ & ' + ' & '.join(f'{c:g}' for c in CD_LISTA) + r' \\', r'\midrule']
    for caso, filas in B.items():
        nom = 'base' if caso == 'base' else r'100 MW, $\tau$ = 5 s'
        L.append(f'Excursión, {nom} [m] & ' + ' & '.join(f"{f['excursion_m']:.3f}" for f in filas) + r' \\')
        L.append(f'H$_2$ quemado, {nom} [kg] & ' + ' & '.join(f"{f['H2_kg']:.2f}" for f in filas) + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('val03_cd'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/val03_cotas.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')
    L = [r'\begin{tabular}{l l r r r}', r'\toprule', r'Escenario & Cantidad & base & mín. & máx. \\',
         r'\midrule']
    for esc in ('ideal', 'real'):
        for k, nom in (('dm_star', r'$\dot m^*$ [kg/s]'), ('P_MW', 'potencia de rechazo [MW]'),
                       ('H2_t', 'H$_2$ por tonelada [kg]'), ('Q_GJ_t', 'calor por tonelada [GJ]')):
            m, b0 = D[esc]['rango'][k], D[esc]['base'][k]
            L.append(f"{esc} & {nom} & {b0:.3g} & {m['min']:.3g} & {m['max']:.3g} \\\\")
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('val03_incertidumbre'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/val03_cotas.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


def graficar(A, B, C):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    et = ['FlyZero Al\n(2026)', 'FlyZero Al\n(2035)', 'assumed', 'Ghatos 2023\n(forced air)', 'bound']
    va = [A['tau_FlyZero_s'][0], A['tau_FlyZero_s'][1], A['tau_supuesto_s'], np.mean(A['tau_Ghatos_forzada_s']),
          A['tau_cota_s']]
    err = [0, 0, 0, (A['tau_Ghatos_forzada_s'][1] - A['tau_Ghatos_forzada_s'][0]) / 2.0, 0]
    ax[0].bar(range(5), va, yerr=err, color=['#0072B2', '#0072B2', 'black', '#E69F00', 'gray'], capsize=4)
    ax[0].set_xticks(range(5), et, fontsize=8)
    ax[0].set_ylabel('Time constant $\\tau$ [s]')
    ax[0].set_title('Condenser lag: estimates vs assumption')
    for caso, ls in (('base', '-o'), ('condensador_100MW_tau5', '--s')):
        f = B[caso]
        ax[1].plot([x['Cd'] for x in f], [100 * x['d_excursion_rel'] for x in f], ls,
                   label='ideal bound' if caso == 'base' else '100 MW, $\\tau$ = 5 s')
    ax[1].axhline(0, color='gray', lw=0.8)
    ax[1].set_xlabel('Vertical drag coefficient $C_d$')
    ax[1].set_ylabel('Change in peak excursion [%]')
    ax[1].set_title('Drag sensitivity')
    ax[1].legend()
    t = np.linspace(0, max(C['dur_pasos_s'], C['dur_minimax_s']) * 1.1, 200)
    for dT, col in zip(DT_EQ, ('#56B4E9', '#0072B2', '#D55E00')):
        dm = C['M_total_kg'] * dT * (1 - np.exp(-t / C['tau_th_s'][1])) / p.T1
        ax[2].plot(t / 60, dm, color=col, label=f'$\\Delta T_{{eq}}$ = {dT:g} K')
    ax[2].axvline(C['dur_minimax_s'] / 60, color='gray', ls=':')
    ax[2].axvline(C['dur_pasos_s'] / 60, color='gray', ls='--')
    ax[2].axhline(p.m_carga, color='black', lw=0.8)
    ax[2].set_xlabel('Time into the maneuver [min]')
    ax[2].set_ylabel('Buoyancy drift [kg]')
    ax[2].set_title('Gas superheat drift (mid $h$)')
    ax[2].legend(fontsize=8)
    for a_ in ax:
        a_.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(rutas.figura('val03_cotas'), dpi=150)
    plt.close(fig)


def main():
    A = parte_A()
    B = parte_B()
    C = parte_C(B['base'][1])
    D = parte_D()
    E = parte_E()
    rutas.guardar_datos('val03_cotas', {'tau': A, 'cd': B, 'termica': C, 'incertidumbre': D, 'humedad': E,
                                        'rangos': {'FK': FK, 'FPHI': D['fphi_medido'], 'DT_EQ': DT_EQ, 'H_FACTOR': H_FACTOR}})
    tablas(A, B, C, D)
    graficar(A, B, C)
    print(f"A. tau: FlyZero c/k {A['tau_FlyZero_s'][0]:.1f} / {A['tau_FlyZero_s'][1]:.1f} s; supuesto "
          f"{A['tau_supuesto_s']:g} s; Ghatos forzada {A['tau_Ghatos_forzada_s'][0]:.0f}-{A['tau_Ghatos_forzada_s'][1]:.0f} s")
    for caso, filas in B.items():
        print(f"B. {caso}: " + '; '.join(f"Cd {f['Cd']:g}: exc {f['excursion_m']:.3f} m ({100*f['d_excursion_rel']:+.1f} %), "
                                         f"H2 {f['H2_kg']:.2f} kg ({100*f['d_H2_rel']:+.2f} %)" for f in filas))
    print(f"C. H2 {C['m_H2_kg']:.0f} kg, V {C['V_m3']:.0f} m3, A {C['A_m2']:.0f} m2, Ra {C['Ra']:.2e}, Nu {C['Nu']:.0f}, "
          f"h {C['h_W_m2K']:.2f} W/m2K, tau_th {[round(x/60,1) for x in C['tau_th_s']]} min; maniobra minimax "
          f"{C['dur_minimax_s']:.0f} s, pasos 1 t a 10 MW {C['dur_pasos_s']/60:.1f} min")
    for k, x in C['deriva'].items():
        print(f"   {k:20s} deriva {[round(v) for v in x['dm_kg']]} kg (h x2, x1, x0.5)")
    for esc in ('ideal', 'real'):
        print(f"D. {esc} (base / mín. / máx.):" + (f" phi_max {D['phi_max']:.3f}" if esc == 'real' else ''))
        for k, m in D[esc]['rango'].items():
            print(f"   {k:8s} {D[esc]['base'][k]:.4g} / {m['min']:.4g} / {m['max']:.4g}")
        for k, (lo, hi) in D[esc]['oat'].items():
            print(f"   OAT {k:14s} dm* {lo['dm_star']:.3f}-{hi['dm_star']:.3f}  P {lo['P_MW']:.0f}-{hi['P_MW']:.0f} MW  "
                  f"H2/t {lo['H2_t']:.1f}-{hi['H2_t']:.1f}")
    for f in E['humedo']:
        print(f"E. HR {100*f['HR']:.0f} %: entra {f['agua_extra_kg_por_kgH2']:.2f} kg de agua por kg de H2; H2/t "
              f"{100*f['dH2_t_rel']:+.1f} %, calor/t {100*f['dQ_t_rel']:+.1f} %, q {100*f['dq_rel']:+.1f} %")


if __name__ == '__main__':
    main()
