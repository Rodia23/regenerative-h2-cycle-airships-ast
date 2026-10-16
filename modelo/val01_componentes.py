# -*- coding: utf-8 -*-
# [EN] Model inputs checked against published reference data (Supplementary S4.2, Table S3).
# =============================================================================
# val01_componentes.py — VAL-01: validación por componente contra datos publicados
# (Paso 2 del plan; punto 4 del checklist del paper 1).
#
# Cada bloque de ecuaciones del modelo se compara con su dato de referencia
# publicado, SIN elegir equipo: atmósfera (ISA), gases (CODATA, NIST), termoquímica
# de la combustión (CODATA), estequiometría (IUPAC), saturación del agua
# (IAPWS-IF97 y NIST/IAPWS-95) y electrólisis (termoneutral e IEA 2025).
#
# Historia: la primera corrida (commit a2f66f0) encontró que el modelo usaba valores
# redondeados (9.81, 8.314, 120 MJ/kg, 1:8:9, L_v a 100 °C, H_s 8500 m, 45 kWh/kg).
# Desde el 2026-09-24 parametros.py usa los valores publicados. Este script deja:
#   1) cada constante del modelo frente a su fuente (tabla val01_componentes.tex);
#   2) los SUPUESTOS de modelado que quedan y su efecto (val01_supuestos.tex);
#   3) las cifras de salida con los valores anteriores y con los actuales
#      (val01_cifras.tex), para que el cambio quede documentado.
#
# VALIDACIÓN (dato externo) != VERIFICACIÓN (fórmula vs simulación, verificacion.py).
# El grupo H de verificacion.py comprueba que las implementaciones usadas aquí
# reproduzcan los valores de verificación publicados (IF97, tabla 35) y que los
# datos de referencia sean coherentes entre sí.
# Salidas: resultados/datos/val01_componentes.json y las tres tablas.
# =============================================================================

import io
import sys
import parametros as p
import rutas
import drop_flujo_optimo_perdrop as sim
import exp05_condensacion_tiempo_real as e5
from exp05_condensacion_tiempo_real import p_sat, N_INERTE, q_rech

# ── Fuentes ──────────────────────────────────────────────────────────────────
FUENTES = {
    'ISA':    'ISO 2533:1975 / U.S. Standard Atmosphere 1976 (NOAA-S/T 76-1562)',
    'CODATA': 'Cox, Wagman & Medvedev, CODATA Key Values for Thermodynamics (1989); '
              'R: CODATA 2018 (Tiesinga et al., Rev. Mod. Phys. 93, 2021)',
    'IUPAC':  'IUPAC 2005, pesos atómicos estándar (Wieser, Pure Appl. Chem. 78, 2006): H 1.00794, O 15.9994',
    'IF97':   'IAPWS R7-97(2012), IAPWS-IF97, ec. 30 y tabla 35',
    'NIST':   'NIST Chemistry WebBook, SRD 69 (Lemmon, McLinden & Friend): agua IAPWS-95, '
              'H2 Leachman et al. 2009, N2 Span et al. 2000; consultado 2026-09-24',
    'IEA':    'IEA, Global Hydrogen Review 2025, Assumptions annex, p. 6 (electrólisis 63 % PCI en 2024)',
    'YG03':   'Yamasaki & Goto, Identification of Blimp Dynamics via Flight Tests, Trans. JSASS 46(153) '
              '195-205, 2003, tabla 4 (dirigible Sky Probe-J, 18 m3, esbeltez 2.2)',
}

# ── Datos de referencia (valores publicados, sin redondear) ──────────────────
G0_ISA   = 9.80665           # m/s2
T0_ISA   = 288.15            # K
P0_ISA   = 101325.0          # Pa
RHO0_ISA = 1.225             # kg/m3
R_AIRE_ISA = 287.05287       # J/kg K
LAPSE_ISA  = 0.0065          # K/m (troposfera)
X_O2_ISA   = 0.209476        # fracción molar de O2 en el aire seco (US 1976)
M_AIRE_ISA = 28.9644e-3      # kg/mol
R_CODATA   = 8.314462618     # J/mol K
M_H2_REF   = 2.01588e-3      # kg/mol
M_O2_REF   = 31.9988e-3      # kg/mol
M_H2O_REF  = 18.01528e-3     # kg/mol
DHF_H2O_L  = 285.830e3       # J/mol, -entalpía de formación del agua líquida (CODATA)
DHF_H2O_G  = 241.826e3       # J/mol, -entalpía de formación del vapor (CODATA)
# NIST WebBook (IAPWS-95), saturación: T [°C] -> (p_sat [kPa], h_l [kJ/kg], h_v [kJ/kg], cp_l, cp_v [J/kg K])
NIST_SAT = {25.0: (3.169929, 104.8292, 2546.505, 4181.600, 1911.823),
            40.0: (7.384938, 167.5330, 2573.510, 4179.646, 1931.421),
            60.0: (19.94643, 251.1804, 2608.835, None, None),
            100.0: (101.4180, 419.1662, 2675.570, None, None)}
NIST_RHO_H2_15C = 0.08520506     # kg/m3 a 15 °C y 101.325 kPa (gas real)
NIST_CP_N2_27C  = 1041.357       # J/kg K a 27 °C y 101.325 kPa
IEA_EFIC_ELECTROLISIS = 0.63     # PCI, promedio global 2024
# Masa vertical (masa propia + masa añadida, MLG(2,2) [kg]) de un dirigible real: estimación analítica
# (teoría potencial, como la del modelo) frente a la identificada en ensayo restringido y en vuelo libre.
YG03_MASA_VERTICAL = {'analitica': 38.44, 'restringido': 29.65, 'vuelo_libre': 46.19}

# ── Valores que usaba el modelo hasta el commit a2f66f0 (para documentar el cambio) ──
ANTERIORES = {'g': 9.81, 'R_univ': 8.314, 'M_H2': 0.002016, 'ratio_O2_H2': 8.0, 'Q_comb_H2': 120e6,
              'L_v_agua': 2257e3, 'H_escala': 8500.0, 'E_H2_electrolisis': 45.0, 'w_O2_aire': 0.2314,
              'M_inerte_aire': 0.02816}


def scale_height_isa(h):
    """Altura de escala LOCAL de la densidad en la troposfera ISA: 1 / (g/(R T) - Γ/T)."""
    T = T0_ISA - LAPSE_ISA * h
    return 1.0 / (G0_ISA / (R_AIRE_ISA * T) - LAPSE_ISA / T)


def rho_isa(h):
    T = T0_ISA - LAPSE_ISA * h
    return RHO0_ISA * (T / T0_ISA) ** (G0_ISA / (R_AIRE_ISA * LAPSE_ISA) - 1.0)


def cifras(rho_a, rho_h2, r_o2, n_inerte, E_el, q):
    """Cifras de salida del trabajo (mismas fórmulas del modelo) para un juego de datos de entrada."""
    L = rho_a / rho_h2 + r_o2                                   # L_eff, agua retenida (phi = 1)

    def phi_max(T_C, P):
        ps = p_sat(T_C)
        return max(0.0, min(1.0, 1.0 - n_inerte * ps / (P - ps)))

    x = 1.0 / (1.0 + n_inerte)                                  # fracción molar de H2O en el escape
    lo, hi = 0.0, 300.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if p_sat(mid) < x * p.P_atmos else (lo, mid)
    P_phi = p_sat(p.T_salida_cond) * (1.0 + n_inerte / (1.0 - p.phi_objetivo))
    return {'L_eff': L, 'H2_por_t_kg': 1000.0 / L, 'agua_por_t_kg': (1.0 + r_o2) * 1000.0 / L,
            'q_MJkg': q / 1e6, 'phi_max_40C_1atm': phi_max(p.T_salida_cond, p.P_atmos),
            'T_rocio_1atm_C': 0.5 * (lo + hi), 'P_phi_objetivo_atm': P_phi / p.P_atmos,
            'E_electrolisis_kWhkg': E_el}


def sensibilidad_escala():
    """Simulación del artículo (PF, 1 t, h_max = 1 m) con la H_s actual y con la anterior (8500 m)."""
    H_modelo = p.H_escala
    out = {}
    for nombre, H in (('actual', H_modelo), ('anterior', ANTERIORES['H_escala'])):
        p.H_escala = H
        try:
            r = sim.correr(graficar=False, verbose=False)
        finally:
            p.H_escala = H_modelo
        out[nombre] = {'H_s_m': H, 'error_max_m': r['error_max_m'], 'dh_arriba_m': r['dh_arriba_m'],
                       'dh_abajo_m': r['dh_abajo_m'], 't_estab_s': r['t_estab_s'], 'm_H2_neto_kg': r['m_H2_neto']}
    return out


def q_con(**cambios):
    """q_rech(0, phi_objetivo) con parámetros cambiados temporalmente."""
    viejos = {k: getattr(p, k) for k in cambios}
    try:
        for k, v in cambios.items():
            setattr(p, k, v)
        return q_rech(0.0, p.phi_objetivo)
    finally:
        for k, v in viejos.items():
            setattr(p, k, v)


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

    # ── Valores de referencia derivados ──────────────────────────────────────
    r_o2_ref = 0.5 * M_O2_REF / M_H2_REF
    r_h2o_ref = M_H2O_REF / M_H2_REF
    lhv_ref = DHF_H2O_G / M_H2_REF
    hhv_ref = DHF_H2O_L / M_H2_REF
    hfg25 = (NIST_SAT[25.0][2] - NIST_SAT[25.0][1]) * 1e3
    w_o2_ref = X_O2_ISA * M_O2_REF / M_AIRE_ISA
    m_inerte_ref = (M_AIRE_ISA - X_O2_ISA * M_O2_REF) / (1.0 - X_O2_ISA)
    n_inerte_ref = 0.5 * (1.0 - X_O2_ISA) / X_O2_ISA           # mol N2+Ar por mol H2O
    Hs_isa = scale_height_isa(p.h0)
    E_el_iea = lhv_ref / 3.6e6 / IEA_EFIC_ELECTROLISIS
    E_el_termoneutral = hhv_ref / 3.6e6
    rho_h2_mod = p.P_atmos / (p.R_H2 * p.T1)
    sens = sensibilidad_escala()

    # ── 1) Cada constante frente a su fuente ─────────────────────────────────
    filas = []

    def fila(bloque, magnitud, unidad, modelo, referencia, fuente, nota, tol=0.01, dec=4):
        desv = (modelo - referencia) / referencia
        filas.append({'bloque': bloque, 'magnitud': magnitud, 'unidad': unidad, 'modelo': modelo,
                      'referencia': referencia, 'desvio_rel': desv, 'tol_rel': tol,
                      'dentro': abs(desv) <= tol, 'fuente': fuente, 'nota': nota, 'dec': dec})

    fila('Atmósfera', 'T de referencia', 'K', p.T1, T0_ISA, 'ISA', 'estado de referencia', dec=2)
    fila('Atmósfera', 'P de referencia', 'Pa', p.P_atmos, P0_ISA, 'ISA', 'estado de referencia', dec=0)
    fila('Atmósfera', 'densidad del aire en h0', 'kg/m3', p.rho_aire_std, rho_isa(p.h0), 'ISA',
         'supuesto: estado de referencia ISA a nivel del mar (tabla de supuestos)', dec=4)
    fila('Atmósfera', 'g', 'm/s2', p.g, G0_ISA, 'ISA', 'gravedad estándar', dec=5)
    fila('Atmósfera', 'altura de escala de la densidad en h0', 'm', p.H_escala, Hs_isa, 'ISA',
         'exponencial con la escala local de la ISA', dec=0)
    fila('Gases', 'R universal', 'J/mol K', p.R_univ, R_CODATA, 'CODATA', 'adoptado', dec=4)
    fila('Gases', 'rho_H2 a 15 °C, 1 atm (gas ideal)', 'kg/m3', rho_h2_mod, NIST_RHO_H2_15C, 'NIST',
         'supuesto: gas ideal (tabla de supuestos)', dec=5)
    fila('Gases', 'fracción másica de O2 en el aire', '--', p.w_O2_aire, w_o2_ref, 'ISA', 'adoptado', dec=4)
    fila('Gases', 'masa molar del aire sin O2', 'g/mol', p.M_inerte_aire * 1e3, m_inerte_ref * 1e3, 'ISA',
         'adoptado', dec=3)
    fila('Gases', 'cp del N2 (27 °C)', 'J/kg K', p.cp_N2, NIST_CP_N2_27C, 'NIST', 'balance de q', dec=0)
    fila('Combustión', 'PCI del H2', 'MJ/kg', p.Q_comb_H2 / 1e6, lhv_ref / 1e6, 'CODATA', 'adoptado', dec=2)
    fila('Combustión', 'O2 por kg de H2', 'kg/kg', p.ratio_O2_H2, r_o2_ref, 'IUPAC', 'adoptado', dec=4)
    fila('Combustión', 'H2O por kg de H2', 'kg/kg', 1.0 + p.ratio_O2_H2, r_h2o_ref, 'IUPAC', 'conservación de masa',
         dec=4)
    fila('Combustión', 'calor latente del agua a 25 °C', 'kJ/kg', p.L_v_agua / 1e3, hfg25 / 1e3, 'NIST',
         'adoptado (referencia del PCI)', dec=1)
    fila('Combustión', 'PCS = PCI + (M_H2O/M_H2) L_v', 'MJ/kg', (p.Q_comb_H2 + (1.0 + p.ratio_O2_H2) * p.L_v_agua) / 1e6,
         hhv_ref / 1e6, 'CODATA', 'dato independiente: valida PCI, L_v y la estequiometría juntos', dec=2)
    fila('Combustión', 'cp del agua líquida (25 °C)', 'J/kg K', p.cp_agua_liq, NIST_SAT[25.0][3], 'NIST',
         'balance de q', dec=0)
    fila('Combustión', 'cp del vapor saturado (25 °C)', 'J/kg K', p.cp_vapor, NIST_SAT[25.0][4], 'NIST',
         'balance de q (solo el vapor que no condensa)', dec=0)
    for T in (40.0, 60.0, 100.0):
        fila('Saturación', f'p_sat a {T:g} °C (IF97)', 'kPa', p_sat(T) / 1e3, NIST_SAT[T][0], 'NIST',
             'phi, rocío y presión del condensador', tol=1e-3, dec=4)
    fila('Saturación', 'N inerte por mol de H2O (escape)', 'mol/mol', N_INERTE, n_inerte_ref, 'ISA',
         'dilución del vapor', dec=4)
    fila('Electrólisis', 'energía por kg de H2', 'kWh/kg', p.E_H2_electrolisis, E_el_iea, 'IEA',
         f'piso termoneutral {E_el_termoneutral:.1f} kWh/kg', dec=1)

    # ── 2) Supuestos que quedan y su efecto ──────────────────────────────────
    L_ideal = p.rho_aire_std / rho_h2_mod + p.ratio_O2_H2
    L_real = p.rho_aire_std / NIST_RHO_H2_15C + p.ratio_O2_H2
    q_act = q_rech(0.0, p.phi_objetivo)
    q_solo_latente = p.Q_comb_H2 + p.phi_objetivo * (1.0 + p.ratio_O2_H2) * p.L_v_agua
    q_isa = q_con(T_aire_exterior=T0_ISA - 273.15)
    dH = sens['anterior']['error_max_m'] - sens['actual']['error_max_m']
    supuestos = [
        {'supuesto': 'H2 como gas ideal', 'alternativa': 'gas real (NIST)', 'cifra': 'H2 por tonelada',
         'cambio_rel': L_ideal / L_real - 1.0},
        {'supuesto': 'aire a 1.225 kg/m3 en h0 (ISA nivel del mar)', 'alternativa': f'ISA a {p.h0:g} m',
         'cifra': 'arrastre (L no cambia)', 'cambio_rel': rho_isa(p.h0) / p.rho_aire_std - 1.0},
        {'supuesto': f'aire exterior a {p.T_aire_exterior:g} °C (día caluroso)', 'alternativa': 'ISA, 15 °C',
         'cifra': 'q', 'cambio_rel': q_isa / q_act - 1.0},
        {'supuesto': 'q con calor sensible (balance completo)', 'alternativa': 'solo latente',
         'cifra': 'q', 'cambio_rel': q_solo_latente / q_act - 1.0},
        {'supuesto': 'electrólisis al 63 % del PCI (IEA 2024)', 'alternativa': '45 kWh/kg (optimista)',
         'cifra': 'tiempo de regeneración', 'cambio_rel': ANTERIORES['E_H2_electrolisis'] / p.E_H2_electrolisis - 1.0},
        {'supuesto': f'H_s local ISA ({p.H_escala:.0f} m)', 'alternativa': '8500 m (anterior)',
         'cifra': 'excursión máx. de la simulación del artículo', 'cambio_abs_mm': 1000.0 * dH},
    ]

    # ── 3) Cifras con los valores anteriores y con los actuales ──────────────
    A = ANTERIORES
    rho_h2_ant = p.P_atmos / (A['R_univ'] / A['M_H2'] * p.T1)
    n_in_ant = (A['ratio_O2_H2'] / A['w_O2_aire'] - A['ratio_O2_H2']) / A['M_inerte_aire'] * A['M_H2']
    q_ant = A['Q_comb_H2'] + p.phi_objetivo * (1.0 + A['ratio_O2_H2']) * A['L_v_agua']   # fórmula de entonces
    antes = cifras(p.rho_aire_std, rho_h2_ant, A['ratio_O2_H2'], n_in_ant, A['E_H2_electrolisis'], q_ant)
    ahora = cifras(p.rho_aire_std, rho_h2_mod, p.ratio_O2_H2, N_INERTE, p.E_H2_electrolisis, q_act)
    # coherencia: las cifras "ahora" son las del modelo
    assert abs(ahora['phi_max_40C_1atm'] - e5.phi_max(p.T_salida_cond)) < 1e-12

    # ── 4) Masa añadida (inercia vertical): fórmula analítica frente a medición en vuelo ──
    yg = YG03_MASA_VERTICAL
    masa_anadida = {'fuente': FUENTES['YG03'], 'datos_kg': yg, 'k_modelo': p.k_masa_anadida,
                    'medida_sobre_analitica': {c: yg[c] / yg['analitica'] for c in ('restringido', 'vuelo_libre')}}
    # si la inercia real fuera r veces la analítica, los flujos dinámicos escalan con 1/sqrt(r)
    masa_anadida['efecto_flujo'] = {c: 1.0 / r ** 0.5 for c, r in masa_anadida['medida_sobre_analitica'].items()}

    salida = {'filas': filas, 'fuentes': FUENTES, 'supuestos': supuestos, 'masa_anadida': masa_anadida,
              'cifras_antes': antes, 'cifras_ahora': ahora, 'anteriores': ANTERIORES,
              'sensibilidad_escala': sens,
              'referencias_derivadas': {
                  'r_O2_H2': r_o2_ref, 'r_H2O_H2': r_h2o_ref, 'LHV_Jkg': lhv_ref, 'HHV_Jkg': hhv_ref,
                  'h_fg_25C_Jkg': hfg25, 'w_O2_aire': w_o2_ref, 'M_inerte_kgmol': m_inerte_ref,
                  'N_inerte': n_inerte_ref, 'H_escala_isa_h0_m': Hs_isa, 'rho_isa_h0': rho_isa(p.h0),
                  'E_el_iea_kWhkg': E_el_iea, 'E_el_termoneutral_kWhkg': E_el_termoneutral,
                  'efic_iea_PCI': IEA_EFIC_ELECTROLISIS, 'M_H2O_sobre_M_H2': r_h2o_ref},
              'n_filas': len(filas), 'n_dentro': sum(f['dentro'] for f in filas)}
    rutas.guardar_datos('val01_componentes', salida)
    tabla(filas)
    tabla_supuestos(supuestos)
    tabla_cifras(antes, ahora)

    print('VAL-01: validación por componente (modelo vs dato publicado)\n')
    for f in filas:
        est = 'ok    ' if f['dentro'] else 'DESVIO'
        print(f'  {est} [{f["bloque"]}] {f["magnitud"]}: modelo {f["modelo"]:.6g} vs {f["referencia"]:.6g} '
              f'{f["unidad"]} ({100 * f["desvio_rel"]:+.3f} %, criterio {100 * f["tol_rel"]:g} %) -> {f["nota"]}')
    print(f'\n{salida["n_dentro"]}/{len(filas)} dentro del criterio')
    print('\nSupuestos que quedan y su efecto:')
    for s in supuestos:
        ef = (f'{100 * s["cambio_rel"]:+.2f} %' if 'cambio_rel' in s else f'{s["cambio_abs_mm"]:+.2f} mm')
        print(f'  {s["supuesto"]} -> {s["alternativa"]}: {s["cifra"]} {ef}')
    print('\nCifras de salida (valores anteriores -> actuales):')
    for k in antes:
        print(f'  {k:<24} {antes[k]:12.5g} -> {ahora[k]:12.5g}  ({100 * (ahora[k] - antes[k]) / antes[k]:+.2f} %)')
    print('\nMasa añadida: masa vertical medida / analítica (Yamasaki & Goto 2003) y efecto en los flujos:')
    for c, r in masa_anadida['medida_sobre_analitica'].items():
        print(f'  {c:<12} {r:.3f}  -> flujo x{masa_anadida["efecto_flujo"][c]:.3f}')
    print('\nSimulación del artículo (PF 1 t, h_max = 1 m) según la altura de escala:')
    for k, s in sens.items():
        print(f'  {k:<8} H_s = {s["H_s_m"]:7.0f} m: excursión máx. {s["error_max_m"]:.4f} m, '
              f'H2 neto {s["m_H2_neto_kg"]:.3f} kg')


def _n(x, dec):
    s = f'{x:,.{dec}f}' if abs(x) >= 10000 else f'{x:.{dec}f}'
    return s.replace(',', '{,}')


def _latex(txt):
    """Nombres de la consola -> notación LaTeX (en orden: lo largo antes que lo corto)."""
    for a, b in (('rho_H2', r'$\rho_{H_2}$'), ('p_sat', r'$p_\text{sat}$'), ('L_v', r'$L_v$'),
                 ('M_H2O/M_H2', r'$M_{H_2O}/M_{H_2}$'), ('H_s', r'$H_s$'), ('kg/m3', r'kg/m$^3$'),
                 ('h0', r'$h_0$'),
                 ('cp del', r'$c_p$ del'), ('H2O', r'H$_2$O'), ('O2', r'O$_2$'), ('N2', r'N$_2$'),
                 ('H2', r'H$_2$'), ('%', r'\%')):
        txt = txt.replace(a, b)
    return txt


def _guardar(nombre, L):
    with open(rutas.tabla(nombre), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/val01_componentes.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


def tabla(filas):
    enc = r'Bloque & Magnitud & Unidad & Modelo & Publicado & Desvío & Fuente \\'
    L = [r'\begin{tabular}{l p{4.4cm} l r r r l}', r'\toprule', enc, r'\midrule']
    previo = None
    for f in filas:
        if previo is not None and f['bloque'] != previo:
            L.append(r'\addlinespace')
        b = f['bloque'] if f['bloque'] != previo else ''
        previo = f['bloque']
        d = f['dec']
        u = {'kg/m3': r'kg/m$^3$', 'm/s2': r'm/s$^2$'}.get(f['unidad'], f['unidad'])
        L.append(f"{b} & {_latex(f['magnitud'])} & {u} & {_n(f['modelo'], d)} & {_n(f['referencia'], d)} & "
                 f"{100 * f['desvio_rel']:+.2f}\\,\\% & {f['fuente']}" + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    _guardar('val01_componentes', L)


def tabla_supuestos(supuestos):
    L = [r'\begin{tabular}{p{5.0cm} p{3.2cm} p{3.6cm} r}', r'\toprule',
         r'Supuesto del modelo & Alternativa & Cifra afectada & Cambio \\', r'\midrule']
    for s in supuestos:
        ef = (f'{100 * s["cambio_rel"]:+.2f}\\,\\%' if 'cambio_rel' in s else f'{s["cambio_abs_mm"]:+.1f}~mm')
        L.append(f'{_latex(s["supuesto"])} & {_latex(s["alternativa"])} & {_latex(s["cifra"])} & {ef}' + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    _guardar('val01_supuestos', L)


def tabla_cifras(antes, ahora):
    nombres = [('L_eff', r'$L_\text{eff}$ (agua retenida)', '--', 2),
               ('H2_por_t_kg', r'H$_2$ por tonelada liberada', 'kg/t', 1),
               ('agua_por_t_kg', 'agua por tonelada liberada', 'kg/t', 0),
               ('q_MJkg', r'calor a rechazar $q$ ($\eta = 0$, $\phi = 0.99$)', 'MJ/kg', 1),
               ('phi_max_40C_1atm', r'$\phi_\text{max}$ a 40\,°C y 1 atm', '--', 3),
               ('T_rocio_1atm_C', 'rocío del escape a 1 atm', '°C', 1),
               ('P_phi_objetivo_atm', r'presión para $\phi = 0.99$ a 40\,°C', 'atm', 1),
               ('E_electrolisis_kWhkg', 'energía de electrólisis', 'kWh/kg', 1)]
    L = [r'\begin{tabular}{l l r r r}', r'\toprule',
         r'Cifra de salida & Unidad & Valores anteriores & Datos publicados & Cambio \\', r'\midrule']
    for k, txt, u, d in nombres:
        L.append(f'{txt} & {u} & {_n(antes[k], d)} & {_n(ahora[k], d)} & '
                 f'{100 * (ahora[k] - antes[k]) / antes[k]:+.1f}\\,\\%' + r' \\')
    L += [r'\bottomrule', r'\end{tabular}']
    _guardar('val01_cifras', L)


if __name__ == '__main__':
    main()
