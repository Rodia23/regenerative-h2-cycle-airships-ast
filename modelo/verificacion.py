# -*- coding: utf-8 -*-
# [EN] Automated verification of the whole pipeline; it gates the export of the numbers.
# =============================================================================
# verificacion.py — VER-01: verificación interna del modelo del paper 1
# como parte del pipeline. Corre después de los experimentos y antes de
# exportar_numeros.py; si algún chequeo falla, sale con error y correr_todo lo
# marca. Cada chequeo compara lo que la simulación o un script OBTUVO con lo que
# la física o una fórmula cerrada ESPERA, con una tolerancia declarada.
#
# Grupos:
#   C  consistencia entre experimentos (mismos datos por caminos distintos)
#   E  dinámica del condensador en las simulaciones del artículo (EXP-09)
#   F  liberación por escalones con retardo del condensador (EXP-10)
#   G  cajas negras: balance de energía y regeneración con retención real (EXP-11)
#   H  referencias de la validación (VAL-01): la implementación de p_sat reproduce
#      la tabla de verificación de IF97, los datos publicados son coherentes entre
#      sí y las constantes del modelo son coherentes entre ellas
#   I  comparación con helio (EXP-12): minimax generalizado, trabajo de compresión de gas
#      real vs ideal, densidad del He, y coherencia con EXP-06
#   J  masa añadida en la dinámica (EXP-13): inercia de equilibrio, factor de la forma
#      cerrada y bajada de la simulación del artículo = D_lim
#   K  tablas NIST-JANAF (janaf.py): nodos, ajustes de Shomate de NIST y PCI de las
#      entalpías de formación
#   L  derivaciones del paper 1 contra comprobaciones independientes (derivaciones.py)
#   M  calor por nivel de temperatura (EXP-14): calibración con FlyZero, bandas, Dalton
#   N  balance de masa por caja negra y tablas nuevas del paper 1
#   O  validación de la condensación (Ec. 5) contra experimentos H2-aire publicados (VAL-02)
#   Q  verificación numérica (VER-02): código vs solución exacta, convergencia del paso y la
#      tolerancia, incertidumbre numérica de la configuración adoptada
#   R  cotas y sensibilidades de la PIRT (VAL-03): tau con datos (FlyZero + aluminio; Ghatos),
#      arrastre, deriva térmica del gas y propagación de incertidumbre
#   P  auditoría de las fuentes de validación (fuentes.py): protocolo completo, PDF sin cambios
#      (SHA-256), datos del modelo = datos registrados, transcripciones sin editar
#   S  admisión de aire al casco (req_aire_casco.py): simulación vs ṁ/rho_H2, escalones
#   V  balance neto por tonelada (exp16_balance_neto.py): conservación de materia, identidad de la
#      cuenta de energía, coherencia con EXP-11 y balance_masa, compresión contra el DOE
#   U  operación en altura (exp15_altura.py): presión contra integración hidrostática, nivel del mar
#      = EXP-14, tendencias con la altura
#   T  irradiancia sobre el casco (irradiancia_casco.py): factor de forma contra integración
#      directa y casos analíticos, H0 contra el anexo 7 de SENAMHI, anual rearmado vs tabla 2
#   W  datos de las naves (S3): fuentes, flotación, inercia del casco y flujo que cabe en la carga
#   X  requisito de masa del rechazo y brecha (req_masa_rechazo.py): presupuesto por kg/s, curva
#      contra el flujo que cabe, candidatos con fuente, nota de la presión del 99 %
#   Z  control reactivo puro frente a la quema limpia (reactivo_vs_limpia.py): bisección (cumple el umbral y 0.1 %
#      menos no), controlador por encima de la quema reactiva ideal, razón ideal = su forma cerrada
#   Y  evaluación de credibilidad (credibilidad.py, NASA-STD-7009B): elasticidades exactas de la ecuación 7
#      contra diferencias numéricas, filas de la tabla de sensibilidad contra VAL-03 y las tablas 3, 5 y 6,
#      niveles 0-4, zona de criticidad del Apéndice D y conteo de fuentes por pedigrí
# Salidas: resultados/datos/verificacion.json, resultados/tablas/verificacion.tex
# =============================================================================

import io
import sys
import math
import parametros as p
BASE = p.NAVE_BASE['corto']     # sufijo corto del caso base (S3, 2026-09-28)
import rutas
import numpy as np
import inercia
import drop_flujo_optimo_perdrop as sim_art
import os
import janaf
import fuentes
from exp05_condensacion_tiempo_real import q_rech, p_sat

CH = []
R_AGUA = 1.0 + p.ratio_O2_H2        # kg de agua por kg de H2 (M_H2O/M_H2)
# Ajustes de Shomate de NIST Chemistry WebBook (SRD 69; Chase 1998), consulta 2026-09-24:
# especie: [((T_min, T_max) [K], (A, B, C, D, E, F, G, H))], Cp = A + Bt + Ct^2 + Dt^3 + E/t^2, t = T/1000.
# Solo para el grupo K (representación independiente de las tablas que usa janaf.py).
SHOMATE = {
    'N2': [((500.0, 2000.0), (19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900, 0.0)),
           ((2000.0, 3000.0), (35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 224.9810, 0.0))],
    'O2': [((300.0, 700.0), (31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0)),
           ((700.0, 2000.0), (30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0))],
    'Ar': [((298.0, 3000.0), (20.78600, 2.825911e-7, -1.464191e-7, 1.092131e-8, -3.661371e-8, -6.197350, 179.9990,
                               0.0))],
    'H2O': [((500.0, 1700.0), (30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, 223.3967, -241.8264)),
            ((1700.0, 3000.0), (41.96426, 8.622053, -1.499780, 0.098119, -11.15764, -272.1797, 219.7809, -241.8264))],
    'CO2': [((298.0, 1200.0), (24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224)),
            ((1200.0, 3000.0), (58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393.5224))],
    'H2': [((298.0, 1000.0), (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0)),
           ((1000.0, 2500.0), (18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133, 0.0))],
}


def chequeo(grupo, nombre, obtenido, esperado, tol_rel=None, tol_abs=None, unidad=''):
    err_abs = abs(obtenido - esperado)
    err_rel = err_abs / abs(esperado) if esperado else (0.0 if err_abs == 0 else float('inf'))
    ok = bool((tol_rel is not None and err_rel <= tol_rel) or (tol_abs is not None and err_abs <= tol_abs))
    tol = f'{100 * tol_rel:g} %' if tol_rel is not None else f'{tol_abs:g} {unidad}'.strip()
    CH.append({'grupo': grupo, 'nombre': nombre, 'obtenido': obtenido, 'esperado': esperado,
               'err_abs': err_abs, 'err_rel': err_rel, 'tolerancia': tol, 'unidad': unidad, 'ok': ok})


def condicion(grupo, nombre, cumple, detalle):
    CH.append({'grupo': grupo, 'nombre': nombre, 'obtenido': float('nan'), 'esperado': float('nan'),
               'err_abs': float('nan'), 'err_rel': float('nan'), 'tolerancia': detalle, 'unidad': '',
               'ok': bool(cumple)})


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    kWh = 3.6e6

    # ── C: consistencia entre experimentos ───────────────────────────────────
    r5 = rutas.leer_datos('exp05_condensacion_tiempo_real')
    r6 = rutas.leer_datos('exp06_requisitos_minimos')
    chequeo('C', 'q condensando todo (EXP-06 = fórmula de EXP-05)', r6['q_MJkg'],
            q_rech(0.0, r6['phi_objetivo']) / 1e6, tol_rel=1e-9, unidad='MJ/kg')
    # R-04 (revisión 2026-10-01): caso coherente a 1 atm; la quema limpia escala exactamente como 1/L
    from exp02_quema_centrada import L_ef as L_efc
    import minimax as mmc
    u1 = r6['base_1t_1atm']
    chequeo('C', 'caso a 1 atm: phi = phi_max a T_salida (EXP-05)', u1['phi'], r5['phi_max_salida'], tol_rel=1e-12)
    chequeo('C', 'caso a 1 atm: ṁ* L(phi_max) = ṁ* (cota ideal) L_eff', u1['dm_star_kgps'] * L_efc(u1['phi']),
            r6['base_1t_h1']['dm_star_kgps'] * mmc.L_eff, tol_rel=1e-9, unidad='kg/s')
    chequeo('C', 'caso a 1 atm: P_cond = ṁ* q(0, phi_max)', u1['P_cond_MW'],
            u1['dm_star_kgps'] * q_rech(0.0, u1['phi']) / 1e6, tol_rel=1e-12, unidad='MW')
    # ── E: dinámica del condensador en las simulaciones del artículo (EXP-09) ─
    r9 = rutas.leer_datos('exp09_condensacion_dinamica')
    for x in r9['casos_base'] + [y for y in r9['plataformas'] if y['modo'] != 'ideal']:
        etiqueta = x['caso'] if 'plataforma' not in x else f"{x['corto']} {x['caso']}"
        cerrado = x['agua_retenida_kg'] + x['agua_venteada_kg'] + x['vapor_pulmon_fin_kg']
        chequeo('E', f'{etiqueta}: retenida + venteada + pulmón = r x H2 quemado', cerrado,
                x['agua_producida_kg'], tol_rel=0.02, unidad='kg')
    mp = rutas.leer_datos('minimax_plataformas')['single_1t']
    for y in r9['plataformas']:
        if y['modo'] == 'ideal':
            chequeo('E', f"{y['corto']}: sin condensador = artículo (H2 neto)", y['H2_neto_kg'],
                    mp[y['plataforma']]['m_H2_neto'], tol_rel=0.02, unidad='kg')
    t0 = next(x for x in r9['casos_base'] if x['caso'] == 'tau_0s')
    chequeo('E', 'sin retardo y con capacidad sobrada retiene la fracción de equilibrio phi', t0['frac_retenida'],
            r9['phi_max'], tol_rel=0.02)

    # ── F: liberación por escalones con retardo del condensador (EXP-10) ─────
    r10 = rutas.leer_datos('exp10_escalones_condensador')
    k_H2 = 1000.0 / r10['H2_por_t_ideal_kg'] - (1.0 + p.ratio_O2_H2)      # rho_aire/rho_H2 - 1
    for x in r10['casos']:
        etiqueta = f"{x['corto']} {x['P_instalada_MW']:g} MW {x['caso']}"
        if x['modo'] != 'ideal':
            cerrado = x['agua_retenida_kg'] + x['agua_venteada_kg'] + x['vapor_pulmon_fin_kg']
            chequeo('F', f'{etiqueta}: retenida + venteada + pulmón = r x H2 quemado', cerrado,
                    x['agua_producida_kg'], tol_rel=0.02, unidad='kg')
        if x['asentado']:
            chequeo('F', f'{etiqueta}: balance final (rho_a/rho_H2 - 1) H2 neto + agua = escalón',
                    k_H2 * x['H2_neto_kg'] + x['agua_retenida_kg'], x['escalon_kg'], tol_rel=0.01, unidad='kg')
        else:
            condicion('F', f'{etiqueta}: sin asentarse figura como que no cumple h_max', not x['cumple_h_max'],
                      'lazo inestable')
        if x['caso'] == 'ventear_tau0s':
            chequeo('F', f'{etiqueta}: sin retardo retiene la fracción de equilibrio phi', x['frac_retenida'],
                    r10['phi_max'], tol_rel=0.01)
    # Física: con toda el agua a bordo (ideal) cada kg de H2 da más pesadez que venteando parte del agua,
    # así que el ideal nunca quema más. Detecta el agua fantasma del integrador (M13, 2026-09-27).
    c10 = {(x['corto'], x['P_instalada_MW'], x['caso']): x for x in r10['casos']}
    for (c, P, caso), x in c10.items():
        if caso == 'ideal':
            ven = c10[(c, P, 'ventear_tau0s')]
            condicion('F', f'{c} {P:g} MW: el ideal no quema más H2 por tonelada que ventear sin retardo',
                      x['H2_por_t_kg'] <= ven['H2_por_t_kg'],
                      f"ideal {x['H2_por_t_kg']:.2f} vs ventear {ven['H2_por_t_kg']:.2f} kg/t")

    # ── G: cajas negras, balance de energía y regeneración (EXP-11) ──────────
    import exp11_balance_energia as e11
    r11 = rutas.leer_datos('exp11_balance_energia')
    for a in r11['autoabasto']:
        etiqueta = f"{a['caso']} {a['extremo']}"
        W = a['W_comp_MJkg'] * 1e6
        chequeo('G', f'{etiqueta}: con eta = eta_min el trabajo alcanza justo (razón = 1)',
                e11.razon_autoabasto(a['eta_min'], a['phi'], a['f'], W), 1.0, tol_rel=1e-9)
        chequeo('G', f'{etiqueta}: razón por tonelada = razón por kg de H2', a['E_gen_kWh_t'] / a['E_pide_kWh_t'],
                e11.razon_autoabasto(a['eta_ref'], a['phi'], a['f'], W), tol_rel=1e-9)
    o = rutas.leer_datos('fuente_O2')
    chequeo('G', 'H2 por tonelada ideal (EXP-11) = H2 de descenso del artículo', r11['H2_por_t_ideal_kg'],
            o['m_H2_aire_kg'], tol_rel=0.005, unidad='kg')
    e10 = {(x['plataforma'], x['P_instalada_MW'], x['caso']): x
           for x in rutas.leer_datos('exp10_escalones_condensador')['casos']}
    for g in r11['regeneracion']:
        etiqueta = f"{g['corto']} {g['P_instalada_MW'] or 'forma cerrada'} {g['caso']}"
        chequeo('G', f'{etiqueta}: recuperado + perdido = quemado', g['recuperado_kg_t'] + g['perdido_kg_t'],
                g['H2_por_t_kg'], tol_rel=1e-9, unidad='kg/t')
        if g['caso'] == 'forma_cerrada':
            chequeo('G', f"{etiqueta}: días por tonelada = días del artículo", g['dias_por_t'],
                    o['regeneracion_solar'][g['plataforma']]['dias'], tol_rel=0.005, unidad='d')
        elif g['caso'] != 'ideal':
            x = e10[(g['plataforma'], g['P_instalada_MW'], g['caso'])]
            chequeo('G', f'{etiqueta}: r x recuperado = agua retenida en EXP-10',
                    R_AGUA * g['recuperado_kg_t'] * x['escalon_kg'] / 1000.0, x['agua_retenida_kg'],
                    tol_rel=0.01, unidad='kg')

    # ── I: comparación con helio (EXP-12) ────────────────────────────────────
    import minimax as mm
    r12 = rutas.leer_datos('exp12_helio')
    pf = p.NAVE_BASE
    a = mm.resolver(pf, p.m_carga, p.h_max_precision, p.D_lim_oper)
    b = mm.resolver(pf, p.m_carga, p.h_max_precision, p.D_lim_oper, rho_gas=mm.rho_H2, L_gas=mm.L_eff)
    chequeo('I', 'minimax generalizado con los datos del H2 = minimax del artículo', b['dm_star'], a['dm_star'],
            tol_rel=1e-12, unidad='kg/s')
    for x in r12['almacen']:
        chequeo('I', f"He a {x['P_MPa']:g} MPa: trabajo isotérmico de gas real (NIST) ~ gas ideal R T ln(P2/P1)",
                x['W_iso_MJkg'], x['W_iso_ideal_MJkg'], tol_rel=0.03, unidad='MJ/kg')
    chequeo('I', 'rho_He de gas ideal = NIST (15 °C, 1 atm)', r12['rho_He'], r12['rho_He_NIST'], tol_rel=1e-3,
            unidad='kg/m3')
    # revisión 2026-10-06: la maniobra limpia dura lo mismo con cualquier gas (t = raiz(6 M_I D (1+s)^3 / (g m)))
    tH_ = 1000.0 / mm.L_eff / r12['base']['H2']['dm_kgps']
    tHe_ = 1000.0 / r12['rutas'][3]['L'] / r12['base']['He']['dm_kgps']
    chequeo('I', 'la maniobra limpia dura lo mismo con He que con H2 (L se cancela)', tHe_, tH_, tol_rel=1e-6, unidad='s')
    # paso [5] de la revisión 2026-10-06: almacenamiento y devolución de la sustentación, recalculados desde los datos
    from exp12_helio import NIST_HE_15C as _He15
    _d24 = p.TANQUE_DOE24
    _kgPV = _d24['compuesto_kg'] / (_d24['P_bar'] * 1e5 * _d24['V_L'] * 1e-3)
    _cosh = r12['rutas'][3]
    chequeo('I', 'COSH: compuesto = 72.7 kg del DOE x (P V) / (700 bar x 144.2 L), 20 MPa (NIST)',
            r12['rutas'][3]['reserva_compuesto_t'], _kgPV * 20e6 * _cosh['kg_por_t'] / _He15[20.0][0] / 1e3, tol_rel=1e-9,
            unidad='t')
    for _i, _g in enumerate(r12['rutas'][:2]):
        _B, _phi = _g['kg_por_t'], 1.0 - _g['perdida_frac']
        chequeo('I', f'ruta {_i}: agua guardada = r phi B / 1000', _g['agua_m3_por_t'], R_AGUA * _phi * _B / 1000.0,
                tol_rel=1e-9, unidad='m3')
        chequeo('I', f'ruta {_i}: compuesto de los tanques de B kg a 35 MPa (NIST) con el mismo escalado P V',
                _g['reserva_compuesto_t'], _kgPV * 35e6 * _B / p.NIST_RHO_H2_15C[35.0] / 1e3, tol_rel=1e-9, unidad='t')
        chequeo('I', f'ruta {_i}: B rho_a/rho_H2 + B (r phi - 1) = 1000 kg (lo que vuelve al instante y con la '
                'electrólisis)', _B * p.rho_aire_std / _g['rho'] + _B * (R_AGUA * _phi - 1.0), 1000.0, tol_rel=1e-9,
                unidad='kg')
        chequeo('I', f'ruta {_i}: el resto que vuelve con la electrólisis = B (r phi - 1)', _g['resto_kg'],
                _B * (R_AGUA * _phi - 1.0), tol_rel=1e-9, unidad='kg')
    chequeo('I', 'COSH: descomprimir todo el He devuelve la tonelada entera', _cosh['al_instante_kg'], 1000.0,
            tol_rel=1e-9, unidad='kg')
    chequeo('I', 'devolver 1 t de golpe pide los mismos moles de H2 que de He', r12['reserva_todo']['H2_kmol'],
            r12['reserva_todo']['He_kmol'], tol_rel=1e-9, unidad='kmol')
    chequeo('I', 'esos moles = 1000 kg / M_aire (gas ideal)', r12['reserva_todo']['He_kmol'],
            r12['reserva_todo']['aire_kmol'], tol_rel=1e-3, unidad='kmol')
    # auditoría 2 (M1, M3): H2 comprimido y tiempo de devolución
    from exp12_helio import NIST_H2_15C as _H215
    for _P in (1.5, 35.0, 70.0):
        chequeo('I', f'H2 a {_P:g} MPa, 15 °C: densidad de Leachman (CoolProp) = consulta a NIST (parametros)',
                _H215[_P][0], p.NIST_RHO_H2_15C[_P], tol_rel=1e-4, unidad='kg/m3')
    _hc = r12['h2_comprimido']
    chequeo('I', 'H2 comprimido a 35 MPa: trabajo de gas real ~ gas ideal R T ln(P2/P1)', _hc['W_iso_MJkg'],
            _hc['W_iso_ideal_MJkg'], tol_rel=0.05, unidad='MJ/kg')
    chequeo('I', 'H2 comprimido: los mismos moles por tonelada que el He de COSH', _hc['kg_por_t'] / p.M_H2 / 1e3,
            r12['reserva_todo']['He_kmol'], tol_rel=1e-9, unidad='kmol')
    _Rg0 = r12['rutas'][0]['kg_por_t'] * (1.0 - r12['rutas'][0]['perdida_frac'])
    _Rg1 = r12['rutas'][1]['kg_por_t'] * (1.0 - r12['rutas'][1]['perdida_frac'])
    chequeo('I', 'devolución del resto: horas en proporción al H2 regenerado (phi = 0.99 frente a 1 atm)',
            r12['rutas'][0]['t_resto_h'], p.t_regeneracion_h * _Rg0 / _Rg1, tol_rel=1e-9, unidad='h')
    _bmI = rutas.leer_datos('balance_masa')
    _g1 = r12['rutas'][1]
    chequeo('I', 'electrólisis por tonelada de la ruta a 1 atm = balance_masa', _g1['kg_por_t'] * (1.0 - _g1['perdida_frac'])
            * p.E_H2_electrolisis, _bmI['E_el_por_t_kWh'], tol_rel=1e-6, unidad='kWh')
    h2 = r12['rutas'][0]
    chequeo('I', 'ruta H2 con agua retenida: L de EXP-12 = L(phi) de EXP-06', h2['L'], r6['L_phi'], tol_rel=1e-9)
    chequeo('I', 'ruta H2 con agua retenida: q de EXP-12 = q de EXP-06', h2['q_Jkg'] / 1e6, r6['q_MJkg'], tol_rel=1e-9,
            unidad='MJ/kg')
    for e in r12['escalones']:
        x = next(y for y in r6['pasos'] if y['corto'] == BASE and y['P_instalada_MW'] == e['P_MW'])
        chequeo('I', f"H2 a {e['P_MW']:g} MW: minutos por tonelada = tiempo de calor de EXP-06 / carga",
                e['H2_min_por_t'], x['t_calor_min'] / x['carga_t'], tol_rel=1e-6, unidad='min/t')

    # ── J: masa añadida en la dinámica (EXP-13, inercia.py) ──────────────────
    eq = mm.equilibrio(pf, p.m_carga)
    chequeo('J', 'inercia en el equilibrio = (1 + k) rho_aire V_casco (con el casco no lleno de gas)',
            inercia.inercial(eq['M_total'], p.rho_aire_std, pf['V_casco'], eq['V_gas'], pf['k']),
            (1 + pf['k']) * p.rho_aire_std * pf['V_casco'], tol_rel=1e-9, unidad='kg')
    con = mm.resolver(pf, p.m_carga, p.h_max_precision, p.D_lim_oper)['dm_star']
    sin = mm.resolver(pf, p.m_carga, p.h_max_precision, p.D_lim_oper, k_ma=0.0)['dm_star']
    chequeo('J', 'forma cerrada: m* con / sin masa añadida = 1/sqrt(1 + k)', con / sin,
            1 / math.sqrt(1 + pf['k']), tol_rel=0.01)
    ra = sim_art.correr(graficar=False, verbose=False, devolver_serie=True)
    s = ra['serie']
    bajada = p.h0 - float(np.interp(p.intervalo_multidrop - 1e-6, s['t'], s['h']))
    chequeo('J', 'simulación del artículo (caso base, 1 t): bajada al soltar = D_lim (misma inercia en fórmula y simulación)',
            bajada, p.D_lim_oper, tol_rel=0.02, unidad='m')
    condicion('J', 'simulación del artículo (caso base, 1 t) cumple h_max con la masa añadida',
              ra['error_max_m'] <= p.h_max_precision, f'excursión máx. {ra["error_max_m"]:.3f} m')

    # ── K: tablas NIST-JANAF (janaf.py; temperatura de llama de derivaciones.py) ──
    # K1 las tablas JANAF se reproducen en los nodos y coinciden con los ajustes de Shomate de
    #    NIST WebBook (representación independiente de la misma fuente, Chase 1998)
    err_nodo = 0.0
    for e, dat in janaf.DATOS.items():
        for T, cp_t, S_t, H_t in dat[3]:
            err_nodo = max(err_nodo, abs((janaf.h(e, T) - dat[2] * 1e3) / 1e3 - H_t), abs(janaf.s0(e, T) - S_t),
                           abs(janaf.cp(e, T) - cp_t))
    chequeo('K', 'JANAF: h, S y Cp interpolados = valores tabulados en todos los nodos (6 especies)', err_nodo, 0.0,
            tol_abs=1e-9, unidad='kJ/mol, J/mol K')
    err_sh = 0.0
    for e, tramos in SHOMATE.items():
        for (Tlo, Thi), (A, B, C, D, E, F, G, H) in tramos:
            for T in (Tlo + 0.25 * (Thi - Tlo), Tlo + 0.75 * (Thi - Tlo)):
                t = T / 1000.0
                cp_sh = A + B * t + C * t ** 2 + D * t ** 3 + E / t ** 2
                err_sh = max(err_sh, abs(janaf.cp(e, T) / cp_sh - 1.0))
    chequeo('K', 'JANAF interpolado vs ajuste de Shomate (NIST WebBook): Cp, 6 especies, rangos de validez',
            err_sh, 0.0, tol_abs=0.005)
    pci = (janaf.h('H2', 298.15) + 0.5 * janaf.h('O2', 298.15) - janaf.h('H2O', 298.15)) / p.M_H2
    chequeo('K', 'PCI del H2 de las entalpías de formación JANAF = Q_comb_H2 de parametros (CODATA)', pci,
            p.Q_comb_H2, tol_rel=1e-5, unidad='J/kg')
    # ── Q: verificación numérica de la simulación del artículo (VER-02) ─────
    v2 = rutas.leer_datos('ver02_numerica')
    ex = {a['config']: a for a in v2['exacta']}
    fina = ex['paso /2, rtol 1e-8']
    chequeo('Q', 'código vs solución exacta (control apagado): error relativo en v con tolerancia estricta',
            fina['err_rel_v'], 0.0, tol_abs=1e-4)
    chequeo('Q', 'código vs solución exacta (control apagado): error relativo en z con tolerancia estricta',
            fina['err_rel_z'], 0.0, tol_abs=1e-4)
    condicion('Q', 'el error contra la solución exacta baja al apretar la tolerancia (código consistente)',
              fina['err_rel_z'] < ex['anterior (rtol 1e-3)']['err_rel_z'] / 10.0,
              f"{ex['anterior (rtol 1e-3)']['err_rel_z']:.1e} -> {fina['err_rel_z']:.1e}")
    chequeo('Q', 'equilibrio antes de soltar: la nave queda quieta', fina['reposo_antes_m'], 0.0, tol_abs=1e-9,
            unidad='m')
    for caso, filas in v2['convergencia'].items():
        f = {x['config']: x for x in filas}
        chequeo('Q', f'{caso}: las dos configuraciones más finas coinciden (solución convergida)',
                max(f['paso /2, rtol 1e-8']['u_rel'].values()), 0.0, tol_abs=1e-4)
        chequeo('Q', f'{caso}: incertidumbre numérica de la configuración adoptada (excursión, H2, agua)',
                max(f['adoptada (rtol 1e-6)']['u_rel'].values()), 0.0, tol_abs=5e-3)

    # ── R: cotas y sensibilidades de la PIRT (VAL-03) ───────────────────────
    v3 = rutas.leer_datos('val03_cotas')
    ta = v3['tau']
    e14 = rutas.leer_datos('exp14_calor_por_temperatura')
    chequeo('R', 'tau: la k de EXP-14 reproduce los 10 kW/kg de FlyZero a 100 C con aire a la T de diseño',
            e14['k_fria_kWkg_por_K'][0] * (e14['T_rechazo_flyzero_C'][0] - p.T_aire_exterior), p.kWkg_radiador[0],
            tol_rel=1e-9, unidad='kW/kg')
    condicion('R', 'tau supuesto dentro de x0.5-x2 de c/k del radiador de FlyZero de aluminio',
              all(0.5 <= ta['tau_supuesto_s'] / x <= 2.0 for x in ta['tau_FlyZero_s']),
              f"supuesto {ta['tau_supuesto_s']:g} s; c/k {ta['tau_FlyZero_s'][0]:.1f}-{ta['tau_FlyZero_s'][1]:.1f} s")
    condicion('R', 'la cota de 30 s cubre el extremo rápido del intercambiador de laboratorio (Ghatos 2023)',
              ta['tau_cota_s'] >= ta['tau_Ghatos_forzada_s'][0],
              f"cota {ta['tau_cota_s']:g} s; Ghatos {ta['tau_Ghatos_forzada_s'][0]:.0f}-{ta['tau_Ghatos_forzada_s'][1]:.0f} s")
    for caso, filas in v3['cd'].items():
        dmax = max(abs(f['d_excursion_rel']) for f in filas)
        chequeo('R', f'{caso}: el arrastre (Cd 0 a 1) mueve la excursión menos de 2 % (PIRT P5, baja importancia)',
                dmax, 0.0, tol_abs=0.02)
    te = v3['termica']
    d10 = te['deriva']['10K_minimax']['dm_kg'][1]
    chequeo('R', 'gas isotérmico: deriva en la maniobra minimax con dT_eq 10 K < 5 % de la carga',
            d10 / p.m_carga, 0.0, tol_abs=0.05)
    for esc in ('ideal', 'real'):
        for k, m in v3['incertidumbre'][esc]['rango'].items():
            b = v3['incertidumbre'][esc]['base'][k]
            condicion('R', f'incertidumbre {esc}: el valor base de {k} cae dentro de su mínimo-máximo',
                      m['min'] * (1 - 1e-9) <= b <= m['max'] * (1 + 1e-9), f"{b:.4g} en [{m['min']:.4g}, {m['max']:.4g}]")
    # min-max en las esquinas (hallazgo S4): la cota ideal solo depende de la inercia, así que su rango es el de los
    # dos ensayos de Yamasaki (OAT); en el escenario real, el H2 por tonelada solo depende de la recolección
    _ri, _oi = v3['incertidumbre']['ideal']['rango']['dm_star'], v3['incertidumbre']['ideal']['oat']['masa_anadida']
    chequeo('R', 'incertidumbre ideal: mínimo-máximo del flujo = los dos ensayos de YG03 (sujeto y vuelo libre)',
            _ri['max'] - _ri['min'], abs(_oi[1]['dm_star'] - _oi[0]['dm_star']), tol_rel=1e-12, unidad='kg/s')
    _rr, _or = v3['incertidumbre']['real']['rango']['H2_t'], v3['incertidumbre']['real']['oat']['agua_retenida']
    chequeo('R', 'incertidumbre real: mínimo-máximo del H2 por tonelada = los extremos de la recolección medida',
            _rr['max'] - _rr['min'], abs(_or[1]['H2_t'] - _or[0]['H2_t']), tol_rel=1e-12, unidad='kg')
    # el factor de Yamasaki es de la masa vertical TOTAL: VAL-03 y VAL-01 (discusión, -9/+14 %) deben coincidir
    _ef1 = rutas.leer_datos('val01_componentes')['masa_anadida']['efecto_flujo']
    _oat = v3['incertidumbre']['ideal']['oat']['masa_anadida']
    _b3 = v3['incertidumbre']['ideal']['base']['dm_star']
    for _i, _c in enumerate(('restringido', 'vuelo_libre')):
        chequeo('R', f'inercia: VAL-03 escala toda la inercia con el factor de YG03 ({_c}), como VAL-01',
                _oat[_i]['dm_star'] / _b3, _ef1[_c], tol_rel=0.01)

    # ── P: auditoría de las fuentes de validación (fuentes.py, VAL-00) ───────
    import importlib
    for k, f in fuentes.FUENTES.items():
        faltan = [c for c in ('cita', 'doi', 'arbitrada', 'ambito', 'formato', 'erratas', 'consistencia', 'decision',
                              'pedigri')
                  if f.get(c) in (None, '', [])]
        if f.get('pedigri') not in (None, 'A', 'B', 'C'):
            faltan.append('pedigri A/B/C')
        if not (f.get('datos') or f.get('transcripciones')):
            faltan.append('datos')
        condicion('P', f'{k}: los 6 pasos del protocolo están llenos', not faltan,
                  'completo' if not faltan else f'falta: {faltan}')
        if f.get('archivo'):
            # los PDF de terceros no están en git: si falta, la identidad queda sin comprobar (no es falla);
            # si está y su hash no es el registrado, sí (auditoría R2-03)
            ruta = os.path.join(fuentes.RAIZ, f['archivo'])
            hay = os.path.exists(ruta)
            ok = hay and fuentes.sha256(ruta) == f['sha256']
            condicion('P', f'{k}: el PDF, si está, tiene el SHA-256 registrado', ok or not hay,
                      'mismo archivo' if ok else ('NO DISPONIBLE (fuera de git): identidad sin comprobar'
                                                  if not hay else 'EL PDF CAMBIÓ'))
            if not hay:                         # no es falla, pero tampoco aprobado: se cuenta aparte (M-7, auditoría 3)
                CH[-1]['sin_comprobar'] = True
        if f.get('modelo'):
            mod, attr = f['modelo']
            m_ = importlib.import_module(mod)
            dic = getattr(m_, attr) if attr else {n: getattr(m_, n, None) for n in f['datos']}
            attr = attr or 'constantes'
            distintos = [n for n, (v, _) in f['datos'].items() if dic.get(n) != v]
            condicion('P', f'{k}: los datos que usa el modelo ({mod}.{attr}) = los registrados', not distintos,
                      f'{len(f["datos"])} datos iguales' if not distintos else f'distintos: {distintos}')
        for nom, (mod, attr, h) in f.get('transcripciones', {}).items():
            h_act = fuentes.huella_objeto(getattr(importlib.import_module(mod), attr))
            condicion('P', f'{k}: la transcripción {nom} ({mod}.{attr}) no cambió', h_act == h,
                      'misma huella' if h_act == h else 'LA TRANSCRIPCIÓN CAMBIÓ: revisar contra la fuente')

    # ── H: referencias de la validación por componente (VAL-01) ─────────────
    # IAPWS-IF97, tabla 35: valores de verificación de la ecuación de saturación
    for T_K, ps_MPa in ((300.0, 0.353658941e-2), (500.0, 0.263889776e1), (600.0, 0.123443146e2)):
        chequeo('H', f'p_sat IF97 a {T_K:g} K = tabla 35 de IAPWS-IF97', p_sat(T_K - 273.15) / 1e6, ps_MPa,
                tol_rel=1e-8, unidad='MPa')
    v1 = rutas.leer_datos('val01_componentes')['referencias_derivadas']
    chequeo('H', 'datos publicados: PCS - PCI (CODATA) = (M_H2O/M_H2) h_fg a 25 °C (NIST)',
            v1['HHV_Jkg'] - v1['LHV_Jkg'], v1['M_H2O_sobre_M_H2'] * v1['h_fg_25C_Jkg'], tol_rel=1e-3, unidad='J/kg')
    chequeo('H', 'modelo: rho_aire = P/(R_aire T1) con R_aire de la ISA', p.rho_aire_std,
            p.P_atmos / (287.05287 * p.T1), tol_rel=1e-4, unidad='kg/m3')
    chequeo('H', 'modelo: rho_aire/rho_H2 = M_aire/M_H2 (gases ideales a igual T y P)',
            p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1)), 28.9644e-3 / p.M_H2, tol_rel=1e-4)

    # ── L: derivaciones del paper 1, con comprobaciones independientes ──────
    import numpy as _np
    import minimax as _mm
    from exp05_condensacion_tiempo_real import p_sat as _psat, N_INERTE as _N, q_rech as _q
    dv = rutas.leer_datos('derivaciones')['valores']
    chequeo('L', 'q: PCI + latente + sensible de entrada - de salida = q_rech del modelo',
            1e6 * (dv['DerQPCI'] + dv['DerQLatente'] + dv['DerQEntra'] - dv['DerQSale']), dv['q_rech_modelo'],
            tol_rel=1e-9, unidad='J/kg')
    # Quema limpia integrada numéricamente (sin la forma cerrada): quema a ṁ desde -t_a hasta que el H2
    # quemado es m/L; inercia M_I antes de soltar y M_I - m después (hallazgo C1, 2026-09-27)
    MI_L, m_L, L_L, D_L = 1e3 * dv['DerMI'], 1000.0, dv['DerLuno'], p.D_lim_oper

    def _limpia(dm, t_a, n=400001):
        T = m_L / (L_L * dm)                                    # duración total de la quema
        t = _np.linspace(-t_a, T - t_a, n)
        F = _np.where(t < 0.0, -L_L * p.g * dm * (t + t_a) / MI_L,
                      p.g * (m_L - L_L * dm * (t + t_a)) / (MI_L - m_L))
        dt = t[1] - t[0]
        v = _np.concatenate(([0.0], _np.cumsum(0.5 * (F[1:] + F[:-1]) * dt)))
        h = _np.concatenate(([0.0], _np.cumsum(0.5 * (v[1:] + v[:-1]) * dt)))
        i0 = _np.searchsorted(t, 0.0)
        return t, v, h, i0

    t_, v_, h_, i0 = _limpia(dv['DerDmStar'], dv['DerTAnticip'])
    chequeo('L', 'quema limpia integrada: baja D_lim al soltar', -float(h_[i0]), D_L, tol_rel=1e-3, unidad='m')
    chequeo('L', 'quema limpia integrada: termina en reposo con el H2 de equilibrio', float(v_[-1]), 0.0,
            tol_abs=1e-4, unidad='m/s')
    chequeo('L', 'quema limpia integrada: bajada total = (1+s) D_lim', -float(h_[-1]), dv['DerHundTotal'],
            tol_rel=1e-3, unidad='m')
    condicion('L', 'quema limpia integrada: nunca sube por encima del punto de liberación',
              float(h_[i0:].max()) <= float(h_[i0]) + 1e-9, f'máx. después de soltar {h_[i0:].max():+.4f} m')
    chequeo('L', 'quema limpia integrada: t de la quema después de soltar = s t_anticip', float(t_[-1]),
            dv['DerTReactiva'], tol_rel=1e-3, unidad='s')
    # revisión 2026-09-30 (R-07): el paso del apéndice A, |v0| t_r / 3, contra la bajada integrada después de soltar
    chequeo('L', 'quema limpia integrada: bajada después de soltar = |v0| t_r / 3 (apéndice A)',
            float(h_[i0] - h_[-1]), dv['DerHundReactiva'], tol_rel=1e-3, unidad='m')

    # Dos lados (I3 de la auditoría): con ṁ ±1 % la quema que termina en reposo con el H2 de equilibrio
    # necesita hundirse más de D_lim (ṁ menor) o menos (ṁ mayor); se busca t_a con v_final = 0.
    from scipy.optimize import brentq as _brentq
    def _bajada_previa(dm):
        t_a = _brentq(lambda ta: _limpia(dm, ta, 20001)[1][-1], 0.2 * dv['DerTAnticip'], 5.0 * dv['DerTAnticip'])
        tt, vv, hh, ii = _limpia(dm, t_a, 200001)
        return -float(hh[ii])
    condicion('L', 'dos lados: con 0.99 ṁ* la quema limpia necesita hundirse más de D_lim (no cabe)',
              _bajada_previa(0.99 * dv['DerDmStar']) > D_L, f'{_bajada_previa(0.99 * dv["DerDmStar"]):.4f} m')
    condicion('L', 'dos lados: con 1.01 ṁ* la quema limpia cabe en D_lim',
              _bajada_previa(1.01 * dv['DerDmStar']) < D_L, f'{_bajada_previa(1.01 * dv["DerDmStar"]):.4f} m')
    chequeo('L', 'quema limpia: H2 quemado = m/L', dv['DerHdosQuemado'], dv['DerBuno'], tol_rel=1e-9, unidad='kg')

    # Cota con reinyección: quema máxima hasta el pico; el pico integrado = h_max
    dmc = dv['DerCotaDm']
    t_ac = (6.0 * MI_L * D_L / (L_L * p.g * dmc)) ** (1.0 / 3.0)
    tc = _np.linspace(-t_ac, 12.0 * t_ac, 800001)
    Fc = _np.where(tc < 0.0, -L_L * p.g * dmc * (tc + t_ac) / MI_L,
                   p.g * (m_L - L_L * dmc * (tc + t_ac)) / (MI_L - m_L))
    dtc = tc[1] - tc[0]
    vc = _np.concatenate(([0.0], _np.cumsum(0.5 * (Fc[1:] + Fc[:-1]) * dtc)))
    hc = _np.concatenate(([0.0], _np.cumsum(0.5 * (vc[1:] + vc[:-1]) * dtc)))
    ic = _np.searchsorted(tc, 0.0)
    ipk = ic + int(_np.argmax(hc[ic:]))
    chequeo('L', 'cota con reinyección integrada: pico = h_max', float(hc[ipk]), p.h_max_precision, tol_rel=1e-3,
            unidad='m')
    chequeo('L', 'cota con reinyección integrada: H2 quemado de más al pico', 100.0 * (L_L * dmc * (tc[ipk] + t_ac) / m_L - 1.0),
            dv['DerCotaSobrequema'], tol_abs=0.5, unidad='%')
    chequeo('L', 'forma cerrada del documento = minimax.resolver (m*)', dv['DerDmStar'], dv['resolver']['dm_star'],
            tol_rel=1e-9, unidad='kg/s')
    chequeo('L', 'cota del documento = minimax.resolver (cota)', dv['DerCotaDm'], dv['resolver']['dm_cota'],
            tol_rel=1e-9, unidad='kg/s')
    chequeo('L', 'L(1) del documento = L_eff del modelo', dv['DerLuno'], _mm.L_eff, tol_rel=1e-9)
    chequeo('L', 'H2 por tonelada = H2 de descenso del artículo', dv['DerBuno'],
            rutas.leer_datos('fuente_O2')['m_H2_aire_kg'], tol_rel=1e-6, unidad='kg')
    for phi, P, nom in ((dv['DerPhiCuarenta'], p.P_atmos, 'phi_max a 1 atm'),
                        (p.phi_objetivo, 1e3 * dv['DerPReq'], 'presión para phi objetivo')):
        chequeo('L', f'Dalton: tras condensar, presión parcial del vapor = p_sat(40 C) ({nom})',
                (1.0 - phi) / ((1.0 - phi) + _N) * P, _psat(p.T_salida_cond), tol_rel=1e-9, unidad='Pa')
    for f, suf in zip(p.perdida_rad, ('Hoy', 'Futuro')):
        e = dv['DerEtaMin' + suf]
        chequeo('L', f'con eta = eta_min el trabajo paga justo el aire (f = {f})', e * p.Q_comb_H2,
                f * _q(e, dv['DerPhiCuarenta']), tol_rel=1e-9, unidad='J/kg')
    MI = 1e3 * dv['DerMI']
    # (autor, 2026-10-09) con la rigidez real del lazo: el zeta declarado es el amortiguamiento real (antes rho g)
    rg = (p.rho_aire_std - p.P_atmos / (p.R_H2 * p.T1)) * p.g
    chequeo('L', 'K_d da el amortiguamiento real: (rho - rho_H2) g K_d / (2 sqrt(M_I (rho - rho_H2) g K_p)) = zeta',
            rg * dv['DerKd'] / (2.0 * _np.sqrt(MI * rg * p.Kp)), p.zeta, tol_rel=1e-9)
    # revisión 2026-09-30 (R-11): periodo natural con la rigidez real, directo de M_I, K_p y las densidades
    _rgr = (p.rho_aire_std - p.P_atmos / (p.R_H2 * p.T1)) * p.g
    chequeo('L', 'periodo natural del lazo = 2 pi sqrt(M_I / ((rho_aire - rho_H2) g K_p))',
            dv['DerPeriodo'], 2.0 * _np.pi * _np.sqrt(MI / (_rgr * p.Kp)) / 60.0, tol_rel=1e-9, unidad='min')
    Gm, T0 = p.gradiente_ISA, p.T_ref_atm
    expo = p.g / (p.R_aire * Gm) - 1.0
    rho = lambda z: (1.0 - Gm * z / T0) ** expo
    dz = 1e-3
    chequeo('L', 'H_s = -rho/(d rho/dh) de la densidad ISA derivada numéricamente',
            -rho(0.0) / ((rho(dz) - rho(-dz)) / (2 * dz)), dv['DerHs'], tol_rel=1e-6, unidad='m')
    # 2026-10-13: el texto da H_s en km (macro ValidHsKm, una decimal); sale de la misma p.H_escala que ValidHs
    chequeo('L', 'H_s del modelo (macros ValidHs y ValidHsKm) = H_s de la ISA derivada numéricamente', p.H_escala,
            -rho(0.0) / ((rho(dz) - rho(-dz)) / (2 * dz)), tol_rel=1e-6, unidad='m')
    sl = rutas.leer_datos('exp11_balance_energia')['sumidero_LN2']
    chequeo('L', 'Q por tonelada del documento = el de EXP-11', dv['DerQPorT'], sl['Q_por_t_GJ'], tol_rel=1e-6,
            unidad='GJ')
    paso = next(x['paso_max_kg'] for x in rutas.leer_datos('exp06_requisitos_minimos')['pasos']
                if x['plataforma'] == p.nombre_dirigible and x['h_m'] == p.h_max_seguridad and x['P_instalada_MW'] == 10.0)
    chequeo('L', 'escalón por la ley m^{3/2} ~ escalón exacto de EXP-06 (caso base, 10 MW)', dv['DerEscalonAprox'],
            paso, tol_rel=0.05, unidad='kg')
    # pasos de la versión 8 del documento contra lo que calculan los scripts del modelo
    _bm = rutas.leer_datos('balance_masa')
    _t14 = rutas.leer_datos('exp14_calor_por_temperatura')['tina']
    _b16 = rutas.leer_datos('exp16_balance_neto')
    chequeo('L', 'O2 que quedaría a bordo por kg de H2 = r - 1', dv['DerOdos'], dv['DerR'] - 1.0, tol_rel=1e-12)
    chequeo('L', 'tina: agua para guardar la banda fría (10 K) = EXP-14', dv['DerTina'], _t14['agua_t_por_t']['10.0'],
            tol_rel=1e-9, unidad='t')
    chequeo('L', 'quema de referencia = P_ref / q(eta ref, phi 1 atm) de balance_masa', dv['DerQuemaRef'],
            _bm['mdot_kgps'], tol_rel=1e-9, unidad='kg/s')
    chequeo('L', 'electrolizador: B phi E_el / 24 h = balance_masa', dv['DerPpem'] * 1e3, _bm['P_pem_W'], tol_rel=1e-9,
            unidad='W')
    chequeo('L', 'kg de tanque por kg de H2 = balance_masa', dv['DerTanquePorKg'], _bm['tanque_por_kg_H2'],
            tol_rel=1e-12)
    # tanques a 350 bar desde el diseño 2023 del DOE (autor 2026-09-29): datos cruzados entre fuentes
    _d24, _rh = p.TANQUE_DOE24, p.NIST_RHO_H2_15C
    _Vm3 = _d24['V_L'] / 1e3
    chequeo('L', 'tanque 2023: V (rho_NIST 70 MPa - rho 1.5 MPa) a 15 C = H2 útil de la tabla 2 del Record 24006',
            _Vm3 * (_rh[70.0] - _rh[1.5]), _d24['util_kg'], tol_rel=0.01, unidad='kg')
    chequeo('L', 'tanque 2023: V rho_NIST 70 MPa a 15 C = H2 guardado de la tabla 2 del Record 24006',
            _Vm3 * _rh[70.0], _d24['guardado_kg'], tol_rel=0.01, unidad='kg')
    _t350 = p.tanque_a_presion(p.P_tanque_bar)
    condicion('L', 'tanque 2023: el liner por diferencia es positivo y menor que el de 2019 (HDPE de 0.5 cm, 8.0 kg; '
                   'en 2023, PA de 0.38 cm)', 0.0 < _t350['liner_ref_kg'] < 8.0, f"{_t350['liner_ref_kg']:.2f} kg")
    _anl = fuentes.FUENTES['HUA10']['datos']['wt_ANL_350_700'][0]
    chequeo('L', 'tanque: razón 350/700 bar con el diseño 2023 = la de ANL-10/24 (otro análisis, 2010)',
            p.wt_COPV / p.tanque_a_presion(700.0)['w'], _anl[0] / _anl[1], tol_rel=0.02)
    chequeo('L', 'compresor al ritmo del electrolizador = EXP-16', dv['DerCompFlujo'], _b16['compresor']['flujo_kgph'],
            tol_rel=1e-9, unidad='kg/h')
    chequeo('L', 'electricidad de compresión típica = extremo bajo de EXP-16', dv['DerEcomp'], _b16['e_comp_kWhkg'][0],
            tol_rel=1e-12, unidad='kWh/kg')
    from exp12_helio import w_iso_he as _wHe
    chequeo('L', 'He a 20 MPa: trabajo de gas ideal ~ gas real (NIST) de EXP-12', dv['DerWHeIdeal'] * 1e6, _wHe(20.0),
            tol_rel=0.05, unidad='J/kg')
    # casco pasivo, paneles, altura, presión y gas isotérmico (derivaciones de la versión 8)
    chequeo('L', 'casco pasivo del caso base: área x h x dT (calma) = macro del documento', dv['DerAtotalBase']
            * dv['DerHConvCalma'] * dv['DerDTCasco'] / 1e6, dv['DerCasco' + BASE][0], tol_rel=1e-12, unidad='MW')
    _pf = next(x for x in _bm['plataformas'] if x['corto'] == BASE)
    chequeo('L', 'paneles a-Si del caso base en invierno: días por tonelada = balance_masa', dv['DerDiasBase'],
            _pf['pv']['a-Si']['dias_por_t']['invierno'], tol_rel=1e-9)
    chequeo('L', 'paneles a-Si del caso base: masa = balance_masa', dv['DerMasaPvBase'] * 1e3, _pf['pv']['a-Si']['masa_kg'],
            tol_rel=1e-9, unidad='kg')
    chequeo('L', 'electrólisis por tonelada del documento = balance_masa', dv['DerElecT'], _bm['E_el_por_t_kWh'],
            tol_rel=1e-9, unidad='kWh')
    _f15 = next(x for x in rutas.leer_datos('exp15_altura')['filas'] if x['h_m'] == 1000.0)
    chequeo('L', 'altura 1 km: presión del perfil = EXP-15', dv['DerPAlt'] * 1e3, _f15['P_Pa'], tol_rel=1e-9, unidad='Pa')
    chequeo('L', 'altura 1 km: densidad relativa = EXP-15', dv['DerRhoRelAlt'], _f15['rho_rel'], tol_rel=1e-9)
    chequeo('L', 'altura 1 km: phi_max = EXP-15', dv['DerPhiAlt'], _f15['phi'], tol_rel=1e-9)
    _f5 = next(x for x in rutas.leer_datos('exp14_calor_por_temperatura')['filas'] if x['P_atm'] == 5.0)
    chequeo('L', '5 atm: fin de la condensación = EXP-14', dv['DerTFinCinco'], _f5['T_fin_C'], tol_rel=1e-9, unidad='C')
    chequeo('L', '5 atm: trabajo de compresión del escape = EXP-14', dv['DerWCinco'], _f5['W_comp_MJkg'], tol_rel=1e-9,
            unidad='MJ/kg')
    # trabajo mínimo de compresión (M-1, auditoría 3): la forma cerrada contra la integral numérica de los moles de GAS
    # (1 + N sin saturar; N P / (P - p_s) saturado) en ln P, sin usar la fórmula
    from exp05_condensacion_tiempo_real import N_INERTE as _NL
    _e6L = rutas.leer_datos('exp06_requisitos_minimos')
    for _PL, _TL, _WL, _nomL in ((_e6L['P_req_atm'] * p.P_atmos, p.T_salida_cond, _e6L['W_comp_MJkg'], 'phi 0.99 a 40 C'),
                                 (5.0 * p.P_atmos, _f5['T_fin_C'], _f5['W_comp_MJkg'], '5 atm, sin saturar al inicio')):
        _psL = p_sat(_TL)
        _lnP = np.linspace(math.log(p.P_atmos), math.log(_PL), 200001)
        _Pg = np.exp(_lnP)
        _ng = np.where(_Pg <= _psL * (1.0 + _NL), 1.0 + _NL, _NL * _Pg / (_Pg - _psL))
        _Wn = np.trapezoid(_ng, _lnP) * p.R_univ * (_TL + 273.15) / p.M_H2 / 1e6
        chequeo('L', f'compresión mínima ({_nomL}) = integral numérica de la fase gaseosa', _WL, _Wn, tol_rel=1e-6,
                unidad='MJ/kg')
    from exp05_condensacion_tiempo_real import phi_max as _phmx
    chequeo('L', '5 atm: la ecuación (5) en el fin de la condensación devuelve la misma phi (Dalton al revés)',
            _phmx(dv['DerTFinCinco'], 5.0 * p.P_atmos), dv['phi14'], tol_rel=1e-6)
    # F1 (2026-10-13): la ec. (5) es la escala de un movimiento de tirón limitado, D_lim = a^3 / (6 j^2 (1 + s)^3), con
    # a = m g / M_I el salto de aceleración y j = L g m* / M_I el tirón de la quema (sección 4)
    import minimax as _mmF
    _rF = _mmF.resolver(p.nombre_dirigible, p.m_carga, p.h_max_precision, p.D_lim_oper)
    _aF = p.m_carga * p.g / _rF['M_inercial']
    _jF = _mmF.L_eff * p.g * _rF['dm_star'] / _rF['M_inercial']
    chequeo('L', 'ec. (5) como tirón limitado: a^3 / (6 j^2 (1 + s)^3) = D_lim (caso base)',
            _aF ** 3 / (6.0 * _jF ** 2 * (1.0 + _rF['s']) ** 3), p.D_lim_oper, tol_rel=1e-12, unidad='m')
    _te3 = rutas.leer_datos('val03_cotas')['termica']
    for _k in ('5K_minimax', '5K_pasos_10MW_1t'):
        _x = _te3['deriva'][_k]
        chequeo('L', f'gas isotérmico ({_k}): (M/T1) dT (1 - e^(-t/tau)) = deriva de VAL-03', dv['DerKgPorK'] * _x['dT_K']
                * (1.0 - math.exp(-_x['duracion_s'] / _te3['tau_th_s'][1])), _x['dm_kg'][1], tol_rel=1e-9, unidad='kg')
    # la entrega de referencia (A2, 2026-10-13): la nave cargada al salir pesa vacía + a bordo + el gas de EXP-22
    _E22 = rutas.leer_datos('exp22_energia_potencia')['E']
    _te3e = _te3['entrega']
    _Me = p.PLATAFORMAS[_E22['nave']]['m_vacio'] + _E22['carga0_t'] * 1e3 + _E22['gas_kg'][0]
    chequeo('L', 'deriva de la entrega: el gas de equilibrio al salir (VAL-03) = el de EXP-22', _te3e['m_H2_kg'],
            _E22['gas_kg'][0], tol_rel=1e-6, unidad='kg')
    chequeo('L', 'deriva de la entrega: la duración = n sueltas x t por suelta de EXP-22 (6 h)', _te3e['duracion_s'],
            _E22['n'] * _E22['t_paso_s'], tol_rel=1e-12, unidad='s')
    for _k in ('5K_entrega', '10K_entrega'):
        _x = _te3['deriva'][_k]
        chequeo('L', f'gas isotérmico ({_k}): (M/T1) dT (1 - e^(-t/tau)) con la nave cargada al salir', _Me / p.T1
                * _x['dT_K'] * (1.0 - math.exp(-_x['duracion_s'] / _te3e['tau_th_s'][1])), _x['dm_kg'][1], tol_rel=1e-6,
                unidad='kg')
    chequeo('L', 'tiempo térmico de la nave cargada: a Ra alto tau escala como V^(1/3) desde el caso base',
            _te3e['tau_th_s'][1], _te3['tau_th_s'][1] * (_te3e['m_H2_kg'] / _te3['m_H2_kg']) ** (1.0 / 3.0), tol_rel=0.01,
            unidad='s')
    condicion('L', 'deriva de la entrega: en horas el gas llega a la temperatura nueva (t > 5 tau, h x0.5)',
              _te3e['duracion_s'] > 5.0 * _te3e['tau_th_s'][2],
              f"{_te3e['duracion_s'] / 3600:.1f} h frente a 5 x {_te3e['tau_th_s'][2] / 60:.0f} min")
    condicion('L', 'deriva de la entrega: mayor que una suelta (1 t) con 5 K, del orden de varias',
              _te3['deriva']['5K_entrega']['dm_kg'][1] > _E22['s_t'] * 1e3,
              f"{_te3['deriva']['5K_entrega']['dm_kg'][1]:.0f} kg frente a {_E22['s_t'] * 1e3:.0f} kg")

    # ── M: calor por nivel de temperatura (EXP-14) ──────────────────────────
    import exp14_calor_por_temperatura as _e14
    r14 = rutas.leer_datos('exp14_calor_por_temperatura')
    for T_rej, kw, k in zip(r14['T_rechazo_flyzero_C'], (10.0, 20.0), r14['k_fria_kWkg_por_K']):
        chequeo('M', f'la ley kW/kg = k (T - T_aire) reproduce FlyZero a {T_rej:g} °C', k * (T_rej - p.T_aire_exterior),
                kw, tol_rel=1e-9, unidad='kW/kg')
    for x in r14['filas']:
        et = f"{x['P_atm']:g} atm, eta {x['eta']:.2f}, {x['tec']}"
        chequeo('M', f'{et}: banda caliente + banda fría = calor total', x['q_caliente_MJkg'] + x['q_fria_MJkg'],
                x['q_total_MJkg'], tol_rel=1e-9, unidad='MJ/kg')
        if x['eta'] == 0.0 and x['tec'] == '2026':
            ps = _e14.p_sat(x['T_fin_C'])
            chequeo('M', f'{et}: al fin de la banda fría condensa la misma phi (Dalton)',
                    1.0 - _e14.N_INERTE * ps / (x['P_atm'] * p.P_atmos - ps), r14['phi_obj'], tol_rel=1e-6)
    x1 = next(x for x in r14['filas'] if x['P_atm'] == 1.0 and x['eta'] == 0.0 and x['tec'] == '2026')
    chequeo('M', 'a 1 atm la banda fría termina a la temperatura de salida del modelo', x1['T_fin_C'],
            p.T_salida_cond, tol_abs=1e-6, unidad='°C')
    chequeo('M', 'a 1 atm y eta = 0 el calor total = q_rech de EXP-05', x1['q_total_MJkg'] * 1e6,
            q_rech(0.0, r14['phi_obj']), tol_rel=1e-9, unidad='J/kg')
    # (2026-10-13) la banda fría de los dos casos de retención (auditoría, objeción 5.1)
    _bdV = r14['bandas']
    chequeo('M', 'banda fría a la presión que retiene 0.99: q(0, 0.99) = q_rech de EXP-05',
            _bdV['objetivo']['q_MJkg'] * 1e6, q_rech(0.0, p.phi_objetivo), tol_rel=1e-9, unidad='J/kg')
    _psV = _e14.p_sat(_bdV['objetivo']['T_fin_C'])
    chequeo('M', 'a esa presión la banda termina donde condensa 0.99 (Dalton)',
            1.0 - _e14.N_INERTE * _psV / (_bdV['objetivo']['P_atm'] * p.P_atmos - _psV), p.phi_objetivo, tol_rel=1e-6)
    for _cl, _bV in _bdV.items():
        _dqV = _e14.q_total(0.0, _bV['phi'], _bV['T_fin_C']) - _e14.q_total(0.0, 0.0, _bV['T_rocio_C'])
        chequeo('M', f'banda fría ({_cl}): integrada a lo largo de la banda = diferencia de entalpías del balance de la '
                'ec. (4) entre el rocío y el fin (independiente; difieren por L_v(T) contra L_v(25 °C) + cp)',
                _bV['Q_banda_MJkg'] * 1e6, _dqV, tol_rel=0.005, unidad='J/kg')
    condicion('M', 'a la presión que retiene 0.99 la banda fría es una parte mayor del calor que a 1 atm',
              _bdV['objetivo']['frac'] > _bdV['1atm']['frac'],
              f"{100 * _bdV['objetivo']['frac']:.1f} % > {100 * _bdV['1atm']['frac']:.1f} %")

    # ── N: balance de masa y tablas nuevas del paper 1 ──────────────────────
    bm = rutas.leer_datos('balance_masa')
    chequeo('N', 'balance de masa: total = suma de las cajas (mínimo)', bm['total_kg'][0],
            sum(c['min_kg'] for c in bm['cajas']), tol_rel=1e-12, unidad='kg')
    chequeo('N', 'balance de masa: total = suma de las cajas (máximo)', bm['total_kg'][1],
            sum(c['max_kg'] for c in bm['cajas']), tol_rel=1e-12, unidad='kg')
    chequeo('N', 'quema de referencia x q = potencia de rechazo de referencia', bm['mdot_kgps'] * bm['q_Jkg'],
            bm['P_ref_W'], tol_rel=1e-12, unidad='W')
    chequeo('N', 'electrolizador: recupera en 24 h el H2 de una tonelada', bm['P_pem_W'] * 24.0 / 1e3,
            bm['recuperado_por_t_kg'] * p.E_H2_electrolisis, tol_rel=1e-12, unidad='kWh')
    dv2 = rutas.leer_datos('derivaciones')['valores']
    import minimax as _mm2
    d_pf = p.NAVE_BASE
    _r5 = _mm2.resolver(d_pf, 1000.0, p.h_max_seguridad, p.D_lim_oper)
    chequeo('N', 'tabla potencia/escalón: caso base 1 t = P_rej de las derivaciones', _r5['dm_star'] * q_rech(0.0, p.phi_objetivo) / 1e6,
            dv2['DerPrej'], tol_rel=1e-9, unidad='MW')
    # tabla 6 con el flujo IDEAL del caso base (autor, 2026-09-27)
    bi = bm['ideal']
    for j, lado in ((0, 'mínimo'), (1, 'máximo')):
        chequeo('N', f'tabla 6 (flujo ideal): total = suma de las cajas ({lado})', bi['total_kg'][j],
                sum(c['min_kg' if j == 0 else 'max_kg'] for c in bi['cajas']), tol_rel=1e-12, unidad='kg')
    chequeo('N', 'tabla 6: quema ideal = ṁ* de la simulación del caso base (caso base, 1 t)', bi['mdot_kgps'],
            rutas.leer_datos('minimax_plataformas')['single_1t'][p.nombre_dirigible]['dm_opt'], tol_rel=1e-9, unidad='kg/s')
    chequeo('N', 'tabla 6: calor sin exportar trabajo = tabla de potencia por escalón (caso base, 1 t)', bi['P_rej_eta0_W'] / 1e6,
            dv2['DerPrej'], tol_rel=1e-9, unidad='MW')
    chequeo('N', 'tabla 6: conversión (extremo liviano) = eje / kW/kg del Rotax al despegue', bi['cajas'][0]['min_kg'],
            bi['P_eje_W'] / max(w / m for w, m in p.MOTORES_ALTERNATIVOS), tol_rel=1e-12, unidad='kg')
    chequeo('N', 'tabla 6: rechazo (extremo liviano) = calor / 10 kW/kg (FlyZero 2035, sistema)', bi['cajas'][1]['min_kg'],
            bi['P_rej_W'] / (max(p.kWkg_TMS_sistema) * 1e3), tol_rel=1e-12, unidad='kg')
    # arreglo con quemador (eta = 0): mismas funciones, sin motor (auditoría R2-02)
    bq = bi['quemador']
    for j, lado in ((0, 'mínimo'), (1, 'máximo')):
        chequeo('N', f'tabla 6 con quemador: total = suma de las cajas ({lado})', bq['total_kg'][j],
                sum(c['min_kg' if j == 0 else 'max_kg'] for c in bq['cajas']), tol_rel=1e-12, unidad='kg')
    chequeo('N', 'tabla 6 con quemador: calor = el de la tabla de potencia por escalón (eta 0)', bq['P_rej_W'],
            bi['P_rej_eta0_W'], tol_rel=1e-12, unidad='W')
    condicion('N', 'tabla 6: el eje del motor cubre los ventiladores con margen (el motor es una elección)',
              min(bi['eje_sobra_W']) > 0.0, f"sobra {bi['eje_sobra_W'][0] / 1e6:.0f}-{bi['eje_sobra_W'][1] / 1e6:.0f} MW")
    from exp05_condensacion_tiempo_real import phi_max as _phm
    chequeo('N', 'tabla de ambiente: con el aire de diseño condensa phi_max(T_salida)',
            _phm(p.T_aire_exterior + (p.T_salida_cond - p.T_aire_exterior)), _phm(p.T_salida_cond), tol_rel=1e-12)

    # ── O: validación de la condensación contra experimentos publicados (VAL-02) ──
    v2 = rutas.leer_datos('val02_condensacion')
    chequeo('O', 'Ec. 5 reproduce la figura 4 de Zakhvatkin et al. (solo cambia Tetens por IF97)', v2['dif_max_curvas'],
            0.0, tol_abs=0.005)
    coherentes = v2['experimentos'][:5]          # el 6.o tiene una AH de saturación que no cuadra con su T
    desv = [abs(x['esperada_g'] / x['esperada_zak_recalc_g'] - 1.0) for x in coherentes]
    chequeo('O', 'agua esperada = su apéndice A recalculado (experimentos 1-5, peor caso)', max(desv), 0.0, tol_abs=0.03)

    # ── S: admisión de aire al casco (req_aire_casco.py) ─────────────────────
    ra = rutas.leer_datos('req_aire_casco')
    vs = ra['verif_sim']
    chequeo('S', 'caso base simulado: máximo |dV_gas/dt| = ṁ*/rho_H2', vs['max_dVdt_sim_m3s'], vs['aire_formula_m3s'],
            tol_rel=1e-3, unidad='m3/s')
    import minimax as _mm3
    chequeo('S', 'aire por tonelada = (1000/L)/rho_H2', ra['volumen_por_t_m3'], 1000.0 / _mm3.L_eff / ra['rho_H2'],
            tol_rel=1e-12, unidad='m3')
    q99 = q_rech(0.0, p.phi_objetivo)
    for x in ra['filas']:
        chequeo('S', f"{x['corto']}: escalón de {ra['P_escalon_MW']:g} MW -> aire = P/(q rho_H2)", x['aire_escalon_m3s'],
                ra['P_escalon_MW'] * 1e6 / (q99 * ra['rho_H2']), tol_rel=0.01, unidad='m3/s')

    # ── T: irradiancia sobre el casco (irradiancia_casco.py) ─────────────────
    import numpy as _npT
    import irradiancia_casco as ic
    th = _npT.linspace(-_npT.pi / 2, _npT.pi / 2, 200001)
    chequeo('T', 'difusa isotrópica: promedio de (1 + cos b)/2 en la mitad superior = 1/2 + 1/pi',
            float(_npT.trapezoid((1 + _npT.cos(th)) / 2, th) / _npT.pi), ic.F_DIFUSA, tol_rel=1e-8)
    chequeo('T', 'suelo: promedio de (1 - cos b)/2 en la mitad superior = 1/2 - 1/pi',
            float(_npT.trapezoid((1 - _npT.cos(th)) / 2, th) / _npT.pi), ic.F_SUELO, tol_rel=1e-7)
    rng = _npT.random.default_rng(7)
    peor = 0.0
    for _ in range(25):       # directa: forma cerrada vs integración de max(0, n·s) en la mitad superior
        s = rng.normal(size=3); s[2] = abs(s[2]); s /= _npT.linalg.norm(s)
        b = rng.uniform(0, _npT.pi); eje = _npT.array([_npT.sin(b), _npT.cos(b), 0.0])
        h = _npT.cross([0.0, 0.0, 1.0], eje)
        n = _npT.outer(_npT.sin(th), h) + _npT.outer(_npT.cos(th), [0.0, 0.0, 1.0])
        num = float(_npT.trapezoid(_npT.clip(n @ s, 0, None), th) / _npT.pi)
        peor = max(peor, abs(ic.directa_media_semicilindro(s, eje) - num))
    chequeo('T', 'directa: (|s_p| + s_z)/pi = integración directa (25 soles y ejes al azar, peor caso)', peor, 0.0,
            tol_abs=1e-6)
    chequeo('T', 'ecuador en equinoccio, eje N-S (sol en el plano de la sección): R_b = (pi + 2)/(2 pi)',
            ic.Rb_dia(0.0, 80, 0.0), (_npT.pi + 2) / (2 * _npT.pi), tol_rel=2e-3)
    chequeo('T', 'ecuador en equinoccio, eje E-W (sol en el plano del eje): R_b = 2/pi', ic.Rb_dia(0.0, 80, _npT.pi / 2),
            2 / _npT.pi, tol_rel=2e-3)
    for n_, mes, H0_fig in ((162, 'junio', 28.1), (344, 'diciembre', 40.0)):
        chequeo('T', f'H0 a 12 S en {mes} = anexo 7 de SENAMHI (leído de la figura)', ic.H0_dia(-12.0, n_) / 1e6, H0_fig,
                tol_rel=0.02, unidad='MJ/m2')
    ir = rutas.leer_datos('irradiancia_casco')
    chequeo('T', 'anual rearmado (anexo 6 x H0) = tabla 2 de SENAMHI (lectura de figura +-0.01)',
            ir['anual_reconstruido'], ir['anual_tabla'], tol_rel=0.03, unidad='kWh/m2/día')
    for est, e in ir['estaciones'].items():
        condicion('T', f'{est}: el factor de forma cae entre toda directa y toda difusa',
                  min(e['F_fd0'], e['F_fd1']) - 1e-9 <= e['F'] <= max(e['F_fd0'], e['F_fd1']) + 1e-9,
                  f"{e['F']:.3f} en [{e['F_fd0']:.3f}, {e['F_fd1']:.3f}]")
    for x in rutas.leer_datos('balance_masa')['plataformas']:
        dpt = x['pv']['a-Si']['dias_por_t']
        chequeo('T', f"{x['corto']}: días invierno / verano = irradiancia verano / invierno", dpt['invierno'] / dpt['verano'],
                ir['estaciones']['verano']['efectiva'] / ir['estaciones']['invierno']['efectiva'], tol_rel=1e-9)
    import exp14_calor_por_temperatura as _e14
    condicion('M', 'banda caliente de EXP-14 = 20 kW/kg, el único número de FlyZero (decisión del autor)',
              tuple(_e14.KWKG_CALIENTE) == (20.0, 20.0), f'{_e14.KWKG_CALIENTE}')

    # ── U: operación en altura (exp15_altura.py) ─────────────────────────────
    import exp15_altura as _e15
    a15 = rutas.leer_datos('exp15_altura')
    r14 = rutas.leer_datos('exp14_calor_por_temperatura')
    hU = 3000.0
    n = 30000
    PU, TU0 = p.P_atmos, p.T_aire_exterior + 273.15
    for i in range(n):             # dP/dh = -P g / (R T(h)), integración independiente (punto medio)
        hm = (i + 0.5) * hU / n
        PU *= math.exp(-p.g / (p.R_aire * (TU0 - p.gradiente_ISA * hm)) * hU / n)
    chequeo('U', 'presión a 3 km: forma cerrada = integración hidrostática del perfil', _e15.perfil(hU)[1], PU,
            tol_rel=1e-9, unidad='Pa')
    f0 = a15['filas'][0]
    chequeo('U', 'a nivel del mar, H2 por tonelada = EXP-14', f0['H2_por_t_kg'], r14['H2_por_t_kg'], tol_rel=1e-12,
            unidad='kg')
    for i, tec in enumerate(('2026', '2035')):
        x = next(r for r in r14['filas'] if r['P_atm'] == 1.0 and r['eta'] == 0.0 and r['tec'] == tec)
        chequeo('U', f'a nivel del mar, masa de rechazo por kW = EXP-14 ({tec})', f0['kg_por_kW'][i], x['kg_por_kW'],
                tol_rel=1e-9, unidad='kg/kW')
    fs = a15['filas']
    condicion('U', 'con la altura sube el agua retenible y bajan el H2 por tonelada y la sustentación',
              all(b['phi'] > a['phi'] and b['H2_por_t_kg'] < a['H2_por_t_kg'] and b['rho_rel'] < a['rho_rel']
                  for a, b in zip(fs, fs[1:])), f"phi {fs[0]['phi']:.3f} -> {fs[-1]['phi']:.3f}")

    # ── V: balance neto por tonelada (exp16_balance_neto.py) ─────────────────
    b16 = rutas.leer_datos('exp16_balance_neto')
    e11 = rutas.leer_datos('exp11_balance_energia')
    bmV = rutas.leer_datos('balance_masa')
    # con quemador (hallazgo M7): los ventiladores rechazan q(0, phi) en vez de q(eta, phi), y la neta los suma entera
    fB = next(f for f in b16['filas'] if f['caso'] == 'B')
    phiB = fB['materia']['H2_regenerado'] / fB['materia']['H2_quemado']
    chequeo('V', 'B con quemador: ventiladores / ventiladores con motor = q(0, phi) / q(eta, phi)',
            fB['E_ventiladores_quemador'][0] / fB['E_ventiladores'][0], q_rech(0.0, phiB) / q_rech(p.rend_gen, phiB),
            tol_rel=1e-9)
    chequeo('V', 'B con quemador: neta = electrólisis + compresión + ventiladores (extremo alto)',
            fB['E_externa_neta_quemador'][1],
            fB['E_electrolisis'] + max(fB['E_compresion']) + max(fB['E_ventiladores_quemador']), tol_rel=1e-12,
            unidad='kWh/t')
    for f in b16['filas']:
        m = f['materia']
        et = f"{f['caso']}{' ' + f['corto'] if f['corto'] else ''}"
        chequeo('V', f'{et}: materia, H2 + aire = agua retenida + vapor venteado + inertes',
                m['agua_retenida'] + m['vapor_venteado'] + m['inertes'], m['H2_quemado'] + m['aire'], tol_rel=1e-12,
                unidad='kg/t')
        chequeo('V', f'{et}: regeneración, agua retenida = H2 regenerado + O2 liberado',
                m['H2_regenerado'] + m['O2_liberado'], m['agua_retenida'], tol_rel=1e-12, unidad='kg/t')
        chequeo('V', f'{et}: externa neta = electrólisis + compresión - excedente (extremos)',
                f['E_externa_neta'][1], f['E_electrolisis'] + max(f['E_compresion']) - min(f['excedente']),
                tol_rel=1e-12, unidad='kWh/t')
    fB = next(f for f in b16['filas'] if f['caso'] == 'B')
    for a in (x for x in e11['autoabasto'] if x['caso'] == '1atm'):
        i = list(p.perdida_rad).index(a['f'])
        chequeo('V', f"caso B = EXP-11 (1 atm, f = {a['f']}): trabajo generado", fB['W_generado'], a['E_gen_kWh_t'],
                tol_rel=1e-9, unidad='kWh/t')
        chequeo('V', f"caso B = EXP-11 (1 atm, f = {a['f']}): ventiladores (f de avión, cota alta)",
                fB['avion']['E_ventiladores'][i], a['E_pide_kWh_t'], tol_rel=1e-9, unidad='kWh/t')
    # R-03 (revisión 2026-10-01): ventiladores con la f por caída de presión = f / f_avión x los de EXP-11
    from balance_masa import f_caida_presion as _fcpV
    for fV, EV in zip(_fcpV(), fB['E_ventiladores']):
        chequeo('V', f'caso B: ventiladores con la f por caída de presión ({fV:.4f}) = f x calor (EXP-11)', EV,
                fV * fB['avion']['E_ventiladores'][0] / p.perdida_rad[0], tol_rel=1e-9, unidad='kWh/t')
    condicion('V', 'caso B: la neta con la f por caída de presión no pasa de la de avión',
              fB['E_externa_neta'][1] <= fB['avion']['E_externa_neta'][1],
              f"{fB['E_externa_neta'][1]:.0f} <= {fB['avion']['E_externa_neta'][1]:.0f} kWh/t")
    chequeo('V', 'caso B: electrólisis = balance_masa', fB['E_electrolisis'], bmV['E_el_por_t_kWh'], tol_rel=1e-9,
            unidad='kWh/t')
    chequeo('V', 'caso B: excedente / electrólisis = lo que la quema cubre en balance_masa (extremo alto)',
            max(fB['excedente']) / fB['E_electrolisis'], max(bmV['cubre']), tol_rel=1e-9)
    chequeo('V', 'caso B, f de avión: excedente / electrólisis = cubre_avion de balance_masa (extremo alto)',
            max(fB['avion']['excedente']) / fB['E_electrolisis'], max(bmV['cubre_avion']), tol_rel=1e-9)
    regC = {(g['corto'], g['P_instalada_MW']): g for g in e11['regeneracion'] if g['caso'] == 'ventear_tau5s'}
    for f in (x for x in b16['filas'] if x['caso'] == 'C'):
        g = regC[(f['corto'], 10.0)]
        chequeo('V', f"caso C {f['corto']}: H2 quemado = EXP-11 (10 MW, tau 5 s)", f['materia']['H2_quemado'],
                g['H2_por_t_kg'], tol_rel=1e-9, unidad='kg/t')
        chequeo('V', f"caso C {f['corto']}: H2 regenerado = EXP-11 (10 MW, tau 5 s)", f['materia']['H2_regenerado'],
                g['recuperado_kg_t'], tol_rel=1e-9, unidad='kg/t')
    W_ideal = p.R_H2 * 300.0 * math.log(350.0 / 20.0) / 3.6e6
    chequeo('V', 'compresión 20 -> 350 bar: gas ideal isotérmico (300 K) vs DOE con gas real (NIST)', W_ideal,
            p.W_iso_20_350, tol_rel=0.10, unidad='kWh/kg')
    condicion('V', 'electricidad de compresión típica dentro del rango medido por el DOE (1.7-6.4 kWh/kg)',
              1.7 <= b16['e_comp_kWhkg'][0] <= 6.4 and 1.7 <= b16['e_comp_kWhkg'][1] <= 6.4,
              f"{b16['e_comp_kWhkg'][0]:.2f} y {b16['e_comp_kWhkg'][1]:.2f} kWh/kg")
    gc = fuentes.FUENTES['GC25']['datos']['e_comp_kWhkg'][0]
    condicion('V', 'compresión: el rango del modelo (DOE) se superpone con el de Giuffrida & Colbertaldo 2025',
              b16['e_comp_kWhkg'][0] <= gc[1] and gc[0] <= b16['e_comp_kWhkg'][1],
              f"modelo {b16['e_comp_kWhkg'][0]:.1f}-{b16['e_comp_kWhkg'][1]:.1f}; revisión {gc[0]:.1f}-{gc[1]:.1f} kWh/kg")
    condicion('V', 'el compresor pedido no supera el flujo del equipo de referencia con ~20 bar de entrada (su masa es '
              'cota superior)', b16['compresor']['cota_valida'],
              f"{b16['compresor']['flujo_kgph']:.2f} kg/h <= {min(p.COMPRESOR_REF['flujo_kgph'])} kg/h")
    condicion('V', 'compresión: el compresor va al ritmo del electrolizador (comprime todo lo regenerado en el mismo día)',
              abs(b16['compresor']['flujo_kgph'] * p.t_regeneracion_h
                  - next(f for f in b16['filas'] if f['caso'] == 'B')['materia']['H2_regenerado']) < 1e-9,
              f"{b16['compresor']['flujo_kgph']:.3f} kg/h x {p.t_regeneracion_h:g} h = H2 regenerado del caso B")

    # ── W: datos de las naves (hallazgo S3, autor 2026-09-28): fuentes, flotación e inercia del casco ──
    import minimax as mmW
    FU = fuentes.FUENTES
    gy, doc, dur, ea, fw = (FU[k]['datos'] for k in ('GY75', 'DOC37', 'DUR24', 'EASANT', 'FWH26'))
    esperado = {
        'LZ 129 Hindenburg': {'V_casco': gy['LZ129_V_total_ft3'][0] * p.FT3, 'm_vacio': gy['LZ129_vacio_lb'][0] * p.LB,
                              'carga_max': doc['LZ129_carga_lb'][0] * p.LB,
                              'L': gy['LZ129_L_ft'][0] * p.FT, 'D': gy['LZ129_D_ft'][0] * p.FT},
        'Flying Whales LCA60T': {'V_casco': fw['FW_V_helio_m3'][0], 'm_vacio': fw['FW_vacio_t'][0] * 1e3,
                                 'carga_max': fw['FW_carga_t'][0] * 1e3, 'L': fw['FW_L_m'][0], 'D': fw['FW_D_m'][0]},
        'ZR-3 Los Angeles': {'V_casco': gy['ZR3_V_total_ft3'][0] * p.FT3, 'm_vacio': gy['ZR3_vacio_lb'][0] * p.LB,
                             'carga_max': dur['ZR3_carga_kg'][0],
                             'L': gy['ZR3_L_ft'][0] * p.FT, 'D': gy['ZR3_D_ft'][0] * p.FT},
        'Zeppelin NT': {'V_casco': ea['NT_V_envolvente_m3'][0], 'L': ea['NT_L_m'][0], 'D': ea['NT_D_m'][0]},
    }
    condicion('W', 'las naves del modelo son las de las fuentes registradas (y ninguna más)',
              set(esperado) == set(p.PLATAFORMAS), ', '.join(sorted(p.PLATAFORMAS)))
    for n, e in esperado.items():
        d = p.PLATAFORMAS[n]
        distintos = [k for k, v in e.items() if abs(d[k] - v) > 1e-9 * max(1.0, abs(v))]
        condicion('W', f"{d['corto']}: datos de parametros.PLATAFORMAS = los de fuentes.py", not distintos,
                  f"{len(e)} datos iguales" if not distintos else f"distintos: {distintos}")
        eq = mmW.equilibrio(n, p.m_carga)
        condicion('W', f"{d['corto']}: con la liberación de referencia el gas cabe en el casco",
                  eq['V_gas'] <= d['V_casco'] * (1 + 1e-9), f"llenado {100 * eq['llenado']:.1f} %")
        condicion('W', f"{d['corto']}: peso vacío + carga máxima flota con H2 en el casco",
                  d['m_vacio'] + d['carga_max'] <= d['sustentacion'] * (1 + 1e-9),
                  f"{(d['m_vacio'] + d['carga_max']) / 1e3:.1f} t <= {d['sustentacion'] / 1e3:.1f} t")
        chequeo('W', f"{d['corto']}: inercia en equilibrio = (1 + k) rho_aire V_casco", eq['M_I'],
                (1 + d['k']) * p.rho_aire_std * d['V_casco'], tol_rel=1e-9, unidad='kg')
        # la inercia no depende del peso de la nave: con 30 % menos de peso vacío es la misma (entra aire)
        chequeo('W', f"{d['corto']}: la inercia no depende del peso vacío (-30 %)",
                mmW.equilibrio(dict(d, m_vacio=0.7 * d['m_vacio']), p.m_carga)['M_I'], eq['M_I'], tol_rel=1e-9,
                unidad='kg')
        chequeo('W', f"{d['corto']}: área de arrastre = (pi/4) L D", d['Area'], math.pi / 4 * d['L'] * d['D'],
                tol_rel=1e-12, unidad='m2')
    for f_, k_ in ((3.0, 0.804), (4.0, 0.860), (6.0, 0.917)):
        chequeo('W', f'k de Lamb en esbeltez {f_:g} = el de Tuckerman 1926 evaluado', p.k_lamb(f_), k_, tol_abs=5e-4)
    chequeo('W', 'k de referencia (esbeltez 4) = el registrado en fuentes.TK26 (0.86, a dos cifras)',
            round(p.k_masa_anadida, 2), FU['TK26']['datos']['k_masa_anadida'][0], tol_abs=1e-12)
    # soltar la carga y quemar m/L deja fuerza neta cero y la misma inercia (el gas que sale lo reemplaza aire)
    pb = p.NAVE_BASE
    eq0 = mmW.equilibrio(pb, p.m_carga)
    B = p.m_carga / mmW.L_eff
    M1 = eq0['M_total'] - p.m_carga - B + (1 + p.ratio_O2_H2) * B
    Vg1 = eq0['V_gas'] - B / mmW.rho_H2
    chequeo('W', 'caso base: soltar 1 t y quemar m/L deja la fuerza neta en cero', (p.rho_aire_std * Vg1 - M1) / M1,
            0.0, tol_abs=1e-12)
    chequeo('W', 'caso base: soltar 1 t y quemar m/L deja la inercia en (1 + k) rho V_casco',
            inercia.inercial(M1, p.rho_aire_std, pb['V_casco'], Vg1, pb['k']), eq0['M_I'], tol_rel=1e-9, unidad='kg')
    # flujo que cabe en la carga útil (brecha de masa como flujo, autor 2026-09-28)
    biW = rutas.leer_datos('balance_masa')['ideal']
    for nom, c in biW['cabe_en_carga'].items():
        # (2026-10-14, objeción 1.8) la tonelada que se suelta también va en la carga útil
        chequeo('W', f'{nom}: cajas fijas + fracción x lo que escala con el flujo + la tonelada = carga útil del caso base',
                c['fijo_kg'] + c['fraccion'] * c['escala_kg'] + biW['m_kg'], pb['carga_max'], tol_rel=1e-12, unidad='kg')
        chequeo('W', f'{nom}: flujo que cabe = fracción x flujo ideal', c['mdot_kgps'],
                c['fraccion'] * biW['mdot_kgps'], tol_rel=1e-12, unidad='kg/s')
    condicion('W', 'el flujo que cabe en la carga es menor que el requerido en los dos extremos',
              all(c['fraccion'] < 1.0 for c in biW['cabe_en_carga'].values()),
              ', '.join(f"{k}: {c['mdot_kgps']:.2f} kg/s" for k, c in biW['cabe_en_carga'].items()))

    # ── entrega completa por escalones (I-1, auditoría 3): con la carga útil COMPLETA de cada nave, y en el caso base
    # con 10 MW el tiempo lo fija el calor: carga x calor por tonelada / P (antes salía con 30 t en vez de 60) ──
    from exp02_quema_centrada import L_ef as _LefI
    _r6I = rutas.leer_datos('exp06_requisitos_minimos')
    for _x in _r6I['pasos']:
        chequeo('W', f"{_x['corto']}, {_x['P_instalada_MW']:g} MW, h {_x['h_m']:g} m: la entrega usa la carga útil completa",
                _x['carga_t'] * 1e3, p.PLATAFORMAS[_x['plataforma']]['carga_max'], tol_rel=1e-12, unidad='kg')
    _xb = next(x for x in _r6I['pasos'] if x['plataforma'] == p.nombre_dirigible and x['P_instalada_MW'] == 10.0
               and x['h_m'] == p.h_max_seguridad)
    _tq = p.NAVE_BASE['carga_max'] / _LefI(p.phi_objetivo) * q_rech(0.0, p.phi_objetivo) / 10e6
    chequeo('W', 'caso base, 10 MW: duración de la entrega = max(pasos x intervalo, carga x calor por t / P)',
            _xb['t_total_min'] * 60.0, max(_xb['n_pasos'] * p.intervalo_multidrop, _tq), tol_rel=1e-9, unidad='s')

    # ── X: requisito de masa del rechazo y brecha con lo que existe (req_masa_rechazo.py; ruta nueva 2026-09-28) ──
    from balance_masa import masa_buffer
    from exp05_condensacion_tiempo_real import phi_max as phi_maxX, T_rocio as T_rocioX, N_INERTE as N_INX
    rq = rutas.leer_datos('req_masa_rechazo')
    rX, cX = rq['requisitos'], rq['candidatos']
    biX = rutas.leer_datos('balance_masa')['ideal']
    BX = rX['presupuesto_t_por_kgps'] * 1e3
    qX = rX['q_MW_por_kgps'] * 1e6
    chequeo('X', 'presupuesto x flujo ideal + cajas fijas = carga útil del caso base',
            BX * rX['mdot_kgps'] + rX['fijas_t'] * 1e3 + rX['suelta_t'] * 1e3, p.NAVE_BASE['carga_max'], tol_rel=1e-12,
            unidad='kg')
    chequeo('X', 'el presupuesto descuenta la tonelada que se suelta (objeción 1.8)', rX['suelta_t'], 1.0, tol_abs=1e-12,
            unidad='t')
    chequeo('X', 'sin ventiladores, el rechazo pide calor / presupuesto', rX['kWkg_pide_f0'], qX / 1e3 / BX,
            tol_rel=1e-12, unidad='kW/kg')
    chequeo('X', 'en f_max la batería sola llena el presupuesto',
            masa_buffer(rX['f_max'] * qX, biX['t_buf_s'], rX['kWkg_bateria']), BX, tol_rel=1e-9, unidad='kg por kg/s')
    # las dos lecturas coinciden: con el FlyZero 2035 (10 kW/kg, f = 15 %) el flujo que cabe es el de balance_masa
    fz = next(c for c in cX if c['clave'] == 'FZ_sistema_2035')
    xX = BX / (qX / 1e3 / fz['kWkg'][0] + masa_buffer(fz['f'] * qX, biX['t_buf_s'], rX['kWkg_bateria']))
    chequeo('X', 'presupuesto por kg/s y flujo que cabe (quemador, clases livianas) dicen lo mismo', xX,
            biX['cabe_en_carga']['quemador_liviano']['fraccion'], tol_rel=1e-9)
    condicion('X', 'con la batería que mueve los ventiladores, ningún candidato con fuente cumple el requisito',
              not any(c['cumple'] for c in cX), f"{len(cX)} candidatos, ninguno cumple")
    condicion('X', 'con motor no hay reparto: su clase más liviana ya pesa más que el presupuesto',
              rX['conversion_min_t_por_kgps'] > rX['presupuesto_t_por_kgps'],
              f"{rX['conversion_min_t_por_kgps']:.1f} > {rX['presupuesto_t_por_kgps']:.1f} t por kg/s")
    chequeo('X', 'ACC: masa total = suma de la tabla S7 (2879 ton)', rq['acc_masa_t'], 2879.0, tol_abs=1e-9, unidad='t')
    chequeo('X', 'ACC: haces + estructura = 66 % del peso (texto del paper, p. 9)',
            (rq['acc_haces_t'] + p.ACC_masas_t['estructura']) / rq['acc_masa_t'], 0.66, tol_abs=0.005)
    chequeo('X', 'ACC: calor / ventiladores = COP impreso en la tabla 3 (97)', p.ACC_Q_MWth / p.ACC_vent_MWe, 97.0,
            tol_abs=0.5)
    chequeo('X', 'requisitos: agua a retener = m* (1 + r) phi', rX['agua_kgps'],
            rX['mdot_kgps'] * R_AGUA * rX['phi'], tol_rel=1e-12, unidad='kg/s')
    # (2026-10-13) el requisito a 1 atm, sin el compresor del escape (auditoría, objeción 5.2)
    _u1X = rX['una_atm']
    _r6X = rutas.leer_datos('exp06_requisitos_minimos')['base_1t_1atm']
    chequeo('X', 'a 1 atm: la tasa = la de EXP-06 con phi_max en L (ReqDmUnaAtm)', _u1X['mdot_kgps'],
            _r6X['dm_star_kgps'], tol_rel=1e-3, unidad='kg/s')
    chequeo('X', 'a 1 atm: presupuesto x tasa + cajas fijas = carga útil del caso base',
            _u1X['presupuesto_t_por_kgps'] * 1e3 * _u1X['mdot_kgps'] + rX['fijas_t'] * 1e3 + rX['suelta_t'] * 1e3,
            p.NAVE_BASE['carga_max'],
            tol_rel=1e-12, unidad='kg')
    _q1X = q_rech(0.0, _u1X['phi'])
    for _fX, _kX in zip(rX['f_fisico'], _u1X['kWkg_pide_ffis']):
        chequeo('X', f'a 1 atm: requisito con f = {_fX:.3f} = calor / (presupuesto - batería que mueve los ventiladores)',
                _kX, _q1X / 1e3 / (_u1X['presupuesto_t_por_kgps'] * 1e3
                                    - masa_buffer(_fX * _q1X, biX['t_buf_s'], rX['kWkg_bateria'])), tol_rel=1e-12,
                unidad='kW/kg')
    # (2026-10-14, suplementario A2) las tablas S7 y S8 a 1 atm: cada fila desde las fórmulas, no desde evaluar()
    chequeo('X', 'a 1 atm: sin ventiladores, el rechazo pide calor / presupuesto (fila de los rígidos de S8)',
            _u1X['kWkg_pide_f0'], _q1X / 1e3 / (_u1X['presupuesto_t_por_kgps'] * 1e3), tol_rel=1e-12, unidad='kW/kg')
    chequeo('X', 'a 1 atm: en f_max la batería sola llena el presupuesto',
            masa_buffer(_u1X['f_max'] * _q1X, biX['t_buf_s'], rX['kWkg_bateria']), _u1X['presupuesto_t_por_kgps'] * 1e3,
            tol_rel=1e-9, unidad='kg por kg/s')
    chequeo('X', 'a 1 atm: el calor por kg de agua = calor / ((1 + r) phi_max)', _u1X['q_agua_MJkg'] * 1e6,
            _q1X / (R_AGUA * _u1X['phi']), tol_rel=1e-12, unidad='J/kg')
    _rigX = next(c for c in cX if c['clave'] == 'rigidos')
    _rig1X = next(fl for fl in _u1X['filas'] if fl['clave'] == 'rigidos')
    chequeo('X', 'tabla S8 a 1 atm: lo que rinden los rígidos escala con el calor por kg de agua (1 atm / 0.99)',
            _rig1X['kWkg'][0] / _rigX['kWkg'][0], _u1X['q_agua_MJkg'] / rX['q_agua_MJkg'], tol_rel=1e-12, unidad='')
    _accX = next(fl for fl in _u1X['filas'] if fl['clave'] == 'ACC_sistema')
    chequeo('X', 'tabla S8 a 1 atm: el extremo alto de la brecha = pide (batería, f del ACC) / ACC entero '
            '= ReqTempVecesUnaAtm', _accX['brecha_rango'][1], _u1X['brecha_temp_ciclo'][1], tol_rel=1e-12, unidad='')
    chequeo('X', 'tabla S7 a 1 atm: presupuesto por kg/s = el de 0.99 x tasa de 0.99 / tasa de 1 atm',
            _u1X['tabla_s7']['presupuesto_t_por_kgps'], rX['presupuesto_t_por_kgps'] * rX['mdot_kgps'] / _u1X['mdot_kgps'],
            tol_rel=1e-12, unidad='t por kg/s')
    condicion('X', 'a 1 atm, como con 0.99: ninguna clase que condensa a la temperatura del ciclo cumple',
              all(fl['brecha_rango'][0] > 1.0 for fl, c in zip(_u1X['filas'], cX) if c['temp_ciclo']),
              f"brecha mínima {_u1X['brecha_temp_ciclo'][0]:.2f}")
    condicion('X', 'a 1 atm el requisito es mayor que con el calor de 0.99 sin compresor (más calor y menos presupuesto)',
              min(_u1X['kWkg_pide_ffis']) > min(rX['kWkg_pide_ffis']),
              f"{min(_u1X['kWkg_pide_ffis']):.2f} > {min(rX['kWkg_pide_ffis']):.2f} kW/kg")
    _fwSX = next(nv for nv in rutas.leer_datos('fig_presupuesto_masa')['apilado']['naves'] if nv['corto'] == 'FW')
    condicion('X', 'a 14 atm el compresor del escape solo ya pesa más que el presupuesto de lo que escala, por kg/s',
              _fwSX['sin']['t'][1] / _fwSX['sin']['mdot_kgps'] > rX['presupuesto_t_por_kgps'],
              f"{_fwSX['sin']['t'][1] / _fwSX['sin']['mdot_kgps']:.1f} > {rX['presupuesto_t_por_kgps']:.1f} t por kg/s")
    chequeo('X', 'requisitos: fracción retenida a 1 atm y 40 °C = la de EXP-11 (1 atm)', rX['phi_1atm'],
            next(a['phi'] for a in rutas.leer_datos('exp11_balance_energia')['autoabasto'] if a['caso'] == '1atm'),
            tol_rel=1e-12)
    condicion('X', 'nota de la tabla: a la presión de EXP-06 se retiene el 99 % a 40 °C',
              phi_maxX(p.T_salida_cond, rX['P_99_atm'] * p.P_atmos) >= rX['phi'] - 1e-9,
              f"phi = {phi_maxX(p.T_salida_cond, rX['P_99_atm'] * p.P_atmos):.4f} a {rX['P_99_atm']:.0f} atm")
    chequeo('X', 'rocío a esa presión: p_sat(T) = x_H2O P (Dalton)', p_sat(rX['T_rocio_99_C']),
            rX['P_99_atm'] * p.P_atmos / (1.0 + N_INX), tol_rel=1e-6, unidad='Pa')
    chequeo('X', 'tabla 6: lo que deja el presupuesto al rechazo = presupuesto - batería del quemador (f de 2035)',
            rX['queda_rechazo_t_por_kgps'] * 1e3,
            BX - masa_buffer(p.perdida_rad[1] * qX, biX['t_buf_s'], rX['kWkg_bateria']), tol_rel=1e-9,
            unidad='kg por kg/s')
    chequeo('X', 'tabla 6: con eso el rechazo pide el kW/kg de la curva en f de 2035',
            qX / 1e3 / (rX['queda_rechazo_t_por_kgps'] * 1e3), rX['kWkg_pide_f2035'], tol_rel=1e-9, unidad='kW/kg')
    sbX = rX['sens_bateria']
    chequeo('X', 'sensibilidad de la batería: celda Amprius SA88 = pulso de 30 s x voltaje nominal / masa (AMP26)',
            p.kWkg_celda_30s_mejor, 153.0 * 3.45 / 96.5, tol_rel=1e-12, unidad='kW/kg')
    chequeo('X', 'sensibilidad de la batería: paquete = celda x fracción celda/paquete del X-57',
            sbX['kWkg_paquete'], p.kWkg_celda_30s_mejor * p.f_celda_paquete, tol_rel=1e-12, unidad='kW/kg')
    condicion('X', 'sensibilidad de la batería: con la mejor celda baja el requisito, pero sigue sobre el sistema 2035',
              sbX['kWkg_pide_f2035'] < rX['kWkg_pide_f2035'] and sbX['kWkg_pide_f2035'] > p.kWkg_TMS_sistema[1],
              f"{sbX['kWkg_pide_f2035']:.1f} kW/kg (antes {rX['kWkg_pide_f2035']:.1f}); cumplen: {sbX['cumplen']}")
    # aun con el mejor empaquetado de 25 autos eléctricos (LOB20, 0.742) la mejor celda no alcanza el sistema 2035
    _kwb = p.kWkg_celda_30s_mejor * max(p.f_celda_paquete_bev)
    _mb = masa_buffer(p.perdida_rad[1] * qX, biX['t_buf_s'], _kwb)
    condicion('X', 'sensibilidad de la batería: con el mejor empaquetado (LOB20) el requisito sigue sobre el sistema 2035',
              qX / 1e3 / (BX - _mb) > p.kWkg_TMS_sistema[1],
              f"{qX / 1e3 / (BX - _mb):.1f} kW/kg con {_kwb:.2f} kW/kg de paquete")
    # la turbina que solo mueve los ventiladores (I-2, auditoría 3) y la brecha por temperatura (I-3)
    tvX = rX['turbina_ventiladores']
    _kwt = max(w / 1e3 / m for w, m in p.TURBOEJE_TCDS)
    chequeo('X', 'turbina para los ventiladores: la más liviana = potencia / masa del TCDS', tvX['kWkg'][1], _kwt,
            tol_rel=1e-12, unidad='kW/kg')
    chequeo('X', 'turbina para los ventiladores: el rechazo pide q / (presupuesto - f q / kW/kg de la turbina)',
            tvX['kWkg_pide_f2035'][0], qX / 1e3 / (BX - p.perdida_rad[1] * qX / 1e3 / _kwt), tol_rel=1e-12,
            unidad='kW/kg')
    chequeo('X', 'turbina para los ventiladores: H2 extra = f q / (eta LHV)', tvX['H2_extra_kgps_por_kgps'],
            p.perdida_rad[1] * qX / (p.rend_gen * p.Q_comb_H2), tol_rel=1e-12, unidad='kg/s por kg/s')
    _Tr1 = T_rocioX(p.P_atmos)[0]
    condicion('X', 'tabla 7: las clases de avión rechazan sobre el rocío a 1 atm y las del ciclo condensan bajo él',
              all((c['temp_ciclo'] == (c['T'] == 'exhaust' or float(c['T']) <= _Tr1)) for c in cX),
              f"rocío {_Tr1:.1f} °C; " + ', '.join(f"{c['clave']} {c['T']}" for c in cX))
    _cic = [c for c in cX if c['temp_ciclo']]
    condicion('X', 'ninguna clase que condensa a la temperatura del ciclo cumple, ni con la turbina',
              all(c['brecha_rango'][0] > 1.0 for c in _cic),
              ', '.join(f"{c['clave']} {c['brecha_rango'][0]:.1f}x" for c in _cic))
    chequeo('X', 'brecha a la temperatura del ciclo: extremo favorable = min(turbina / lo mejor que rinde)',
            rX['brecha_temp_ciclo'][0],
            min(c['pide_turbina'] / max(c.get('kWkg_favorable', 0.0), c['kWkg'][1]) for c in _cic), tol_rel=1e-12)
    gyX = fuentes.FUENTES['GY75']['datos']
    rg = next(c for c in cX if c['clave'] == 'rigidos')
    chequeo('X', 'rígidos: equipo por kg/s de agua del LZ 129 = masa (p. 20) / consumo a 50 kn (p. 19)',
            rg['t_por_kgps']['LZ129'], gyX['agua_lb'][0]['LZ129'] / (gyX['consumo_lbh'][0]['LZ129'] / 3600.0),
            tol_rel=1e-12, unidad='kg por kg/s')

    # ── Z: control reactivo puro frente a la quema limpia (reactivo_vs_limpia.py) ──
    import reactivo_vs_limpia as rvz
    rz = rutas.leer_datos('reactivo_vs_limpia')
    for nomz, xz in rz['naves'].items():
        cz = xz['corto']
        rvz.sim.usar_nave(nomz, rz['m_carga_kg'])
        chequeo('Z', f'{cz}: tasa de la quema limpia = la de la tabla 3 (minimax_plataformas)', xz['dm_limpia'],
                rutas.leer_datos('minimax_plataformas')['single_1t'][nomz]['dm_opt'], tol_rel=1e-9, unidad='kg/s')
        condicion('Z', f'{cz}: la quema limpia no sube sobre h0', xz['subida_limpia_m'] < 1e-3,
                  f"sube {xz['subida_limpia_m']:.2e} m")
        for kz, hz in xz['objetivos_m'].items():
            dz = xz['reactivo_sim'][kz]
            condicion('Z', f'{cz}: reactivo {kz}: la tasa hallada cumple el umbral y 0.1 % menos no lo cumple',
                      rvz.subida(dz) <= hz < rvz.subida(dz / 1.001), f"{dz:.3f} kg/s para {hz:.3f} m")
            condicion('Z', f'{cz}: reactivo {kz}: el controlador no baja de la quema reactiva ideal',
                      dz >= xz['reactivo_ideal'][kz] * (1 - 1e-9),
                      f"sim {dz:.2f} >= ideal {xz['reactivo_ideal'][kz]:.2f} kg/s")
            chequeo('Z', f'{cz}: reactivo ideal / quema limpia = raiz(4 (1+s)^3 D_lim / (s^2 h)) ({kz})',
                    xz['razon_ideal'][kz], xz['razon_ideal_formula'][kz], tol_rel=1e-6)
        # barrido de la figura del suplementario (2026-10-12)
        chequeo('Z', f'{cz}: excursión de la quema limpia, forma cerrada de la ec. (5) = la simulada',
                xz['excursion_limpia_formula_m'], xz['excursion_limpia_m'], tol_rel=0.015, unidad='m')
        _sz = xz['barrido']['subida_reactivo_m']
        condicion('Z', f'{cz}: barrido reactivo: la subida baja con la tasa instalada',
                  all(b <= a * (1 + 1e-6) for a, b in zip(_sz, _sz[1:])), f'{_sz[0]:.1f} -> {_sz[-1]:.3f} m')
        for kz in ('precision', 'seguridad'):
            # el punto de umbral cae sobre la curva del barrido (interpolación log-log)
            _hz = float(np.exp(np.interp(np.log(xz['reactivo_sim'][kz]), np.log(rz['tasas_barrido_kgps']),
                                         np.log(_sz))))
            chequeo('Z', f'{cz}: barrido reactivo: en la tasa del umbral ({kz}) la curva pasa por el umbral',
                    _hz, xz['objetivos_m'][kz], tol_rel=0.05, unidad='m')

    # (revisión 2026-10-01) hélices existentes de la discusión: energía con la fuerza lineal a cero, por integración
    fz = rutas.leer_datos('helices_vectoriales')
    import helices_vectoriales as hvz
    for casoz in (fz['referencia'], fz['bajo']):
        for hz in casoz['helices']:
            ts = np.linspace(0.0, casoz['t_quema_s'], 20001)
            Pz = np.array([hvz.potencia(fz['F_N'] * (1.0 - t / casoz['t_quema_s']), hz['A_m2']) for t in ts])
            chequeo('Z', f"hélices, {casoz['P_rej_W']/1e6:.0f} MW, discos {hz['A_m2']:.0f} m2: energía = P_max t / 2.5 "
                    f"(integración numérica)", hz['E_J'], float(np.trapezoid(Pz, ts)), tol_rel=1e-6, unidad='J')
            chequeo('Z', f"hélices, {casoz['P_rej_W']/1e6:.0f} MW, discos {hz['A_m2']:.0f} m2: fracción de la energía "
                    f"química = E / (m/L PCI)", hz['E_frac_quim'],
                    hz['E_J'] / (fz['m_kg'] / fz['L'] * p.Q_comb_H2), tol_rel=1e-12)
    chequeo('Z', 'hélices: la quema de 1 t dura lo mismo que los escalones de VAL-03 (10 MW)', fz['referencia']['t_quema_s'],
            rutas.leer_datos('val03_cotas')['termica']['dur_pasos_s'], tol_rel=1e-9, unidad='s')
    chequeo('Z', 'hélices: la duración escala como 1/P_rej', fz['bajo']['t_quema_s'] * p.freno_P_bajo,
            fz['referencia']['t_quema_s'] * p.P_ref_masa, tol_rel=1e-12, unidad='J')
    for Az, Pgz in zip(p.freno_A_discos, fz['escalon_grande']['P_W']):
        chequeo('Z', f'hélices, discos {Az:.0f} m2: la potencia escala como F^1.5 (escalón grande / 1 t)', Pgz,
                hvz.potencia(fz['F_N'], Az) * (fz['escalon_grande']['m_kg'] / fz['m_kg']) ** 1.5, tol_rel=1e-12, unidad='W')
    for Az, Egz, hz in zip(p.freno_A_discos, fz['escalon_grande']['E_J'], fz['referencia']['helices']):
        chequeo('Z', f'hélices, discos {Az:.0f} m2: con el mismo rechazo la energía escala como m^2.5 (escalón grande / 1 t)',
                Egz, hz['E_J'] * (fz['escalon_grande']['m_kg'] / fz['m_kg']) ** 2.5, tol_rel=1e-12, unidad='J')
    # cualquier freno: techo de respuesta, con la inercia por la forma cerrada del casco y la subida libre integrada
    import inercia as inz
    for nz, rz in fz['respuesta'].items():
        dz = p.PLATAFORMAS[nz]
        MIz = inz.inercial_casco(dz['V_casco'], dz['k'])
        chequeo('Z', f"freno, {dz['corto']}: M_I del equilibrio = (1 + k) rho V_casco", rz['M_I_kg'], MIz, tol_rel=1e-9,
                unidad='kg')
        tz = np.linspace(0.0, rz['tau_techo_s'], 20001)
        vz = rz['m_kg'] * p.g / (MIz - rz['m_kg']) * tz
        chequeo('Z', f"freno, {dz['corto']}: la nave libre sube h_max en el techo de respuesta (integración numérica)",
                float(np.trapezoid(vz, tz)), p.h_max_precision, tol_rel=1e-9, unidad='m')
    condicion('Z', 'freno: la nave más chica exige la respuesta más rápida (tau ~ raiz(V/m))',
              fz['respuesta']['Zeppelin NT']['tau_techo_s'] < fz['respuesta']['Flying Whales LCA60T']['tau_techo_s'],
              f"NT {fz['respuesta']['Zeppelin NT']['tau_techo_s']:.2f} s < LCA60T "
              f"{fz['respuesta']['Flying Whales LCA60T']['tau_techo_s']:.2f} s")
    for cz in fz['cota_todas_verticales']:
        chequeo('Z', f"cota todas verticales ({cz['P_inst_W']/1e6:g} MW, {cz['A_m2']:.0f} m2): P(m_max g) = P_inst",
                hvz.potencia(cz['m_max_kg'] * p.g, cz['A_m2']), cz['P_inst_W'], tol_rel=1e-12, unidad='W')

    # (revisión 2026-10-06) el freno como 4.ª caja negra (freno.py; sección nueva del paper 1)
    fr = rutas.leer_datos('freno')
    import minimax as mmF
    for fF in fr['filas']:
        eF = f"freno, {fF['corto']} {fF['m_kg']/1e3:.0f} t, {fF['P_W']/1e6:.0f} MW"
        MIF = mmF.equilibrio(fF['nave'], fF['m_kg'])['M_I']
        sF = math.sqrt(1.0 - fF['m_kg'] / MIF)
        # la ecuación 7 escrita de nuevo con m_e y la s de la liberación completa: su rechazo es P
        dmF = math.sqrt(p.g * fF['m_e_kg'] ** 3 / (6.0 * MIF * fr['L_tasa'] ** 2 * p.D_lim_oper * (1.0 + sF) ** 3))
        chequeo('Z', f"{eF}: la quema limpia de m_e (ec. 7, s de la liberación completa) pide el rechazo instalado",
                dmF * fr['Q_Jkg'], fF['P_W'], tol_rel=1e-9, unidad='W')
        chequeo('Z', f"{eF}: rechazo en tiempo real del escalón = ṁ*(m) q (minimax.limpio)", fF['P_rt_W'],
                mmF.limpio(MIF, fF['m_kg'], p.D_lim_oper)['dm_star'] * fr['Q_Jkg'], tol_rel=1e-12, unidad='W')
        chequeo('Z', f"{eF}: tiempo de recuperación x P = calor del escalón, (m/L) q", fF['t_rec_s'] * fF['P_W'],
                fF['m_kg'] / fr['L_calor'] * fr['Q_Jkg'], tol_rel=1e-12, unidad='J')
        tsF = np.linspace(0.0, fF['t_rec_s'], 20001)
        for AF, EF in zip(p.freno_A_discos, fF['E_hel_J']):
            PzF = np.array([hvz.potencia(fF['F_N'] * (1.0 - t / fF['t_rec_s']), AF) for t in tsF])
            chequeo('Z', f"{eF}, discos {AF:.0f} m2: energía del freno = P_max t / 2.5 (integración numérica)", EF,
                    float(np.trapezoid(PzF, tsF)), tol_rel=1e-6, unidad='J')
        condicion('Z', f"{eF}: el freno toma lo que el rechazo no compensa en tiempo real (0 < m_e < m)",
                  0.0 < fF['m_e_kg'] < fF['m_kg'], f"m_e {fF['m_e_kg']:.1f} kg de {fF['m_kg']:.0f}")
    f10F = next(x for x in fr['filas'] if x['corto'] == 'FW' and x['m_kg'] == 1000.0 and x['P_W'] == p.P_ref_masa)
    f2F = next(x for x in fr['filas'] if x['corto'] == 'FW' and x['m_kg'] == 1000.0 and x['P_W'] == p.freno_P_bajo)
    chequeo('Z', 'freno: 1 t con 10 MW dura lo mismo que en helices_vectoriales (discusión)', f10F['t_rec_s'],
            fz['referencia']['t_quema_s'], tol_rel=1e-9, unidad='s')
    chequeo('Z', 'freno: haces de ACC con 10 MW = los de helices_vectoriales', f10F['masa_acc_haces_t'],
            fz['referencia']['rechazo']['acc_haces_t'], tol_rel=1e-12, unidad='t')
    chequeo('Z', 'freno: haces de ACC con 2 MW = los de helices_vectoriales', f2F['masa_acc_haces_t'],
            fz['bajo']['rechazo']['acc_haces_t'], tol_rel=1e-12, unidad='t')
    # G(mu): la subida con la fuerza (1 + mu) m g mandada a tope y retardo de primer orden, integrada (a = tau = 1)
    from scipy.integrate import solve_ivp as ivpF
    muF = fr['mu']
    paradaF = lambda t, y: y[1]
    paradaF.terminal, paradaF.direction = True, -1
    solF = ivpF(lambda t, y: [y[1], 1.0 - (1.0 + muF) * (1.0 - math.exp(-t))], (0.0, 100.0), [0.0, 0.0],
                events=paradaF, rtol=1e-12, atol=1e-14, max_step=0.01)
    for nF, rF in fr['respuesta'].items():
        chequeo('Z', f"freno, {p.PLATAFORMAS[nF]['corto']}: G(mu) = subida máxima / (a tau^2), integrando la subida",
                solF.y_events[0][0][0], rF['G'], tol_rel=1e-6)
        aF = 1000.0 * p.g / (rF['M_I_kg'] - 1000.0)
        chequeo('Z', f"freno, {p.PLATAFORMAS[nF]['corto']}: plazo de primer orden = raiz(h_max / (a G))",
                rF['tau_primer_orden_s'], math.sqrt(p.h_max_precision / (aF * rF['G'])), tol_rel=1e-9, unidad='s')
    for sF_ in fr['simulacion']:
        cF = p.PLATAFORMAS[sF_['nave']]['corto']
        condicion('Z', f"freno, simulación {cF}, 1 t con 10 MW: la subida no pasa de h_max",
                  sF_['sube_m'] <= p.h_max_precision, f"{sF_['sube_m']:.3f} m <= {p.h_max_precision} m")
        chequeo('Z', f"freno, simulación {cF}, 1 t con 10 MW: H2 quemado = m/L", sF_['H2_kg'], sF_['H2_esperado_kg'],
                tol_rel=0.01, unidad='kg')
        condicion('Z', f"freno, simulación {cF}, 1 t con 10 MW: energía del freno a menos de 10 % de la forma cerrada",
                  abs(sF_['razon_E'] - 1.0) < 0.10, f"simulación / cerrada = {sF_['razon_E']:.3f}")
        condicion('Z', f"freno, simulación {cF}: el freno no satura", not sF_['satura'], str(sF_['satura']))
    # (revisión 2026-10-06, ronda 2) anticipar con el freno: la ecuación 7 con la tasa de fuerza F' = L g ṁ*
    for nA, aA in fr['anticipacion'].items():
        cA = p.PLATAFORMAS[nA]['corto']
        fA = next(x for x in fr['filas'] if x['nave'] == nA and x['m_kg'] == aA['m_kg'] and x['P_W'] == aA['P_W'])
        chequeo('Z', f"freno anticipado, {cA}: rampa = L g ṁ* (ecuación 7)", aA['rampa_Nps'],
                fr['L_tasa'] * p.g * fA['dm_star_kgps'], tol_rel=1e-12, unidad='N/s')
        chequeo('Z', f"freno anticipado, {cA}: el escalón entero llega en t_nat", aA['t_F_s'], fA['t_nat_s'],
                tol_rel=1e-9, unidad='s')
        mA, MA, FpA, t1A = aA['m_kg'], aA['M_I_kg'], aA['rampa_Nps'], aA['t_previa_s']
        s1A = ivpF(lambda t, y: [y[1], -FpA * t / MA], (0.0, t1A), [0.0, 0.0], rtol=1e-12, atol=1e-14, max_step=0.01)
        chequeo('Z', f"freno anticipado, {cA}: la rampa previa hunde D_lim (integración numérica)", -s1A.y[0][-1],
                p.D_lim_oper, tol_rel=1e-6, unidad='m')
        s2A = ivpF(lambda t, y: [y[1], (mA * p.g - FpA * t) / (MA - mA)], (t1A, aA['t_F_s']),
                   [s1A.y[0][-1], s1A.y[1][-1]], rtol=1e-12, atol=1e-14, max_step=0.01)
        chequeo('Z', f"freno anticipado, {cA}: bajada total (1+s) D_lim y en reposo al llegar al escalón",
                -s2A.y[0][-1], aA['bajada_m'], tol_rel=1e-5, unidad='m')
        condicion('Z', f"freno anticipado, {cA}: velocidad nula al llegar al escalón y nunca sube",
                  abs(s2A.y[1][-1]) < 1e-6 and max(s1A.y[0].max(), s2A.y[0].max()) <= 1e-12,
                  f"v final {s2A.y[1][-1]:.2e} m/s")
    from helices_vectoriales import potencia as potenciaP
    for fP in fr['filas']:
        PhP = [potenciaP(fP['m_kg'] * p.g, A) / potenciaP(fP['F_N'], A) for A in p.freno_A_discos]
        chequeo('Z', f"freno {fP['corto']} {fP['m_kg']/1e3:.0f} t, {fP['P_W']/1e6:.0f} MW: sin quema previa, potencia "
                f"x (m/(m - m_e))^1.5", fP['factor_P_sin_previa'], PhP[0], tol_rel=1e-9)
    tz = fr['tau_nt']
    condicion('Z', 'freno, NT: anticipar admite más retardo que esperar (|dh| <= h_max)',
              tz['anticipa_s'] > tz['espera_s'], f"{tz['anticipa_s']:.2f} s > {tz['espera_s']:.2f} s")
    condicion('Z', 'freno, NT: con ese retardo, esperar pasa h_max y anticipar no',
              max(tz['espera_con_tau_ant']) > p.h_max_precision >= max(tz['anticipa_con_tau_ant']),
              f"esperando {max(tz['espera_con_tau_ant']):.2f} m, anticipando {max(tz['anticipa_con_tau_ant']):.2f} m")
    condicion('Z', 'freno, NT: sin quema previa y con el retardo base, |dh| <= h_max',
              max(tz['sin_previa_tau_base']) <= p.h_max_precision, f"{max(tz['sin_previa_tau_base']):.2f} m")
    # f por caída de presión (R-03) y aire de enfriamiento
    rqz = rutas.leer_datos('req_masa_rechazo')['requisitos']
    for dpz, fz_ in zip(p.dp_rechazo_Pa, rqz['f_fisico']):
        chequeo('X', f'f por caída de presión ({dpz:.0f} Pa) = dp / (eta rho cp dT)', fz_,
                dpz / (p.eta_ventilador * p.P_atmos / (p.R_aire * (p.T_aire_exterior + 273.15)) * p.cp_aire *
                       p.dT_aire_enfr), tol_rel=1e-12)
    condicion('X', 'requisito con el f por caída de presión entre f = 0 y el f de FlyZero 2035',
              rqz['kWkg_pide_f0'] <= min(rqz['kWkg_pide_ffis']) <= max(rqz['kWkg_pide_ffis']) <= rqz['kWkg_pide_f2035'],
              f"{rqz['kWkg_pide_f0']:.2f} <= {rqz['kWkg_pide_ffis'][0]:.2f}-{rqz['kWkg_pide_ffis'][1]:.2f} <= "
              f"{rqz['kWkg_pide_f2035']:.2f} kW/kg")
    chequeo('X', 'aire de enfriamiento = P_rej (quemador) / (cp dT)', rqz['aire_enfr_kgps'],
            rutas.leer_datos('balance_masa')['ideal']['P_rej_eta0_W'] / (p.cp_aire * p.dT_aire_enfr), tol_rel=1e-9,
            unidad='kg/s')
    # humedad (R-07): con HR = 0 reproduce el modelo seco; balance de agua; el aire húmedo baja el H2 por tonelada
    import val03_cotas as v3z
    sz = v3z.humedo(0.0)
    chequeo('R', 'humedad: con HR = 0 el calor es q_rech(0, phi_max)', sz['q_MJkg'] * 1e6,
            q_rech(0.0, sz['phi_max']), tol_rel=1e-12, unidad='J/kg')
    from exp02_quema_centrada import L_ef as L_efz
    chequeo('R', 'humedad: con HR = 0 la pesadez es L(phi_max)', sz['L'], L_efz(sz['phi_max']), tol_rel=1e-12)
    for hz in rutas.leer_datos('val03_cotas')['humedad']['humedo']:
        condicion('R', f"humedad {100*hz['HR']:.0f} %: el agua que entra condensa y baja el H2 por tonelada",
                  hz['agua_extra_kg_por_kgH2'] > 0 and hz['dH2_t_rel'] < 0,
                  f"+{hz['agua_extra_kg_por_kgH2']:.2f} kg/kg, H2/t {100*hz['dH2_t_rel']:+.1f} %")

    # ── Y: evaluación de credibilidad (credibilidad.py, CRED-01, NASA-STD-7009B) ──
    import credibilidad as cr
    import minimax as mmY
    from exp02_quema_centrada import L_ef as L_efY
    cY = rutas.leer_datos('credibilidad')
    Eb = cY['elasticidades_base']
    for x in ('m', 'M_I', 'L', 'D'):
        chequeo('Y', f'elasticidad de ṁ* respecto de {x}: analítica = diferencias centrales (caso base)',
                Eb['numerica'][x], Eb['analitica'][x], tol_rel=1e-6)
    for c_, e_ in cY['elasticidades_naves'].items():
        _rY = cr.requisitos(nave=next(v for v in p.PLATAFORMAS.values() if v['corto'] == c_))
        _nY = cr.elasticidades(_rY['M_I'], p.m_carga, p.D_lim_oper, L_efY(1.0))['numerica']
        chequeo('Y', f'{c_}: elasticidad de ṁ* respecto de la carga, analítica = numérica', _nY['m'], e_['m'], tol_rel=1e-6)
    rbY = cY['robustez']
    dY = rutas.leer_datos('dlim_flujo')['plataformas'][p.nombre_dirigible]
    chequeo('Y', 'ṁ* del caso base de CRED-01 = el de la tabla 3 (dlim_flujo)', rbY['base']['dm_star'],
            dY['dm_star_kgps'], tol_rel=1e-9, unidad='kg/s')
    chequeo('Y', 'P_rej del caso base de CRED-01 = la de la tabla 6 con eta 0 (balance_masa)', rbY['base']['P_MW'],
            rutas.leer_datos('balance_masa')['ideal']['P_rej_eta0_W'] / 1e6, tol_rel=1e-6, unidad='MW')
    filY = {f['clave']: f['requisitos'] for f in rbY['filas']}
    v3Y = rutas.leer_datos('val03_cotas')['incertidumbre']
    for k_ in ('dm_star', 'P_MW'):
        chequeo('Y', f'fila de la inercia = rango de VAL-03, cota ideal ({k_}, mínimo)',
                min(r[k_] for r in filY['M_I']), v3Y['ideal']['rango'][k_]['min'], tol_rel=1e-9)
        chequeo('Y', f'fila de la inercia = rango de VAL-03, cota ideal ({k_}, máximo)',
                max(r[k_] for r in filY['M_I']), v3Y['ideal']['rango'][k_]['max'], tol_rel=1e-9)
    for i_, r_ in enumerate(v3Y['real']['oat']['agua_retenida']):
        for k_ in ('dm_star', 'P_MW', 'H2_t'):
            chequeo('Y', f'fila del agua retenida = VAL-03, recolección medida extremo {i_ + 1} ({k_})',
                    filY['phi'][2 - i_][k_], r_[k_], tol_rel=1e-9)
    chequeo('Y', 'fila de D_lim: el flujo escala exactamente como D_lim^(-1/2) (ecuación 7)',
            filY['D'][1]['dm_star'] / filY['D'][0]['dm_star'], (rbY['D_max'] / p.D_lim_oper) ** -0.5, tol_rel=1e-12)
    _P20 = mmY.resolver(p.NAVE_BASE, 20000.0, p.h_max_seguridad, p.D_lim_oper)['dm_star'] * q_rech(0.0, p.phi_objetivo) / 1e6
    chequeo('Y', 'fila de la carga: P_rej con 20 t = la de la tabla 5', filY['m'][1]['P_MW'], _P20, tol_rel=1e-12,
            unidad='MW')
    chequeo('Y', 'fila del volumen del casco: ṁ* con el volumen alternativo = elasticidad de la inercia integrada',
            filY['V'][1]['dm_star'] / filY['V'][0]['dm_star'],
            cr.requisitos(fk=p.V_casco_FW_alternativo / p.NAVE_BASE['V_casco'])['dm_star'] / filY['V'][0]['dm_star'],
            tol_rel=1e-9)
    niv = [f['nivel'] for f in cY['capacidad'] + cY['resultados']]
    condicion('Y', 'niveles enteros de 0 a 4 en los 11 factores (5 de capacidad, 6 de resultados)',
              len(cY['capacidad']) == 5 and len(cY['resultados']) == 6 and all(n_ in range(5) for n_ in niv),
              ', '.join(str(n_) for n_ in niv))
    condicion('Y', 'criticidad: la matriz del Apéndice D da la zona amarilla (consecuencia II-III, influencia 4)',
              cY['criticidad']['zonas'] == ['Y', 'Y'], '/'.join(cY['criticidad']['zonas']))
    condicion('Y', 'uso: nivel 1 por definición en el primer uso del M&S (tabla 15 de la norma)',
              next(f['nivel'] for f in cY['resultados'] if f['clave'] == 'Uso') == 1, 'primer uso')
    _fu = [f for f in fuentes.FUENTES.values() if f.get('ambito') == 'paper 1']
    chequeo('Y', 'fuentes del paper 1 contadas por pedigrí = registro', sum(cY['fuentes_paper1'][g] for g in 'ABC'),
            len(_fu), tol_abs=0)

    # ── Tablas 3 y 9 y figura de los esferoides de la copia (puntos E, H e I, 2026-10-08) ──
    # se leen las celdas de las tablas generadas y se comparan con el modelo por otro camino
    import re as _re

    def _celdas(nombre):
        filas = {}
        for lin in open(rutas.tabla(nombre), encoding='utf-8'):
            lin = lin.strip()
            if ' & ' in lin and lin.endswith(r'\\'):
                c = [x.strip() for x in lin[:-2].split(' & ')]
                filas[c[0]] = c[1:]
        return filas

    def _num(celda):
        t = celda.replace('{,}', '').replace('$-$', '-')
        return float(_re.search(r'[-+]?\d+(\.\d+)?', t).group(0))

    import minimax as _mmE
    _m1 = rutas.leer_datos('minimax_plataformas')['single_1t']
    _t3 = _celdas('art_cota_ideal_rev2')
    _kr3 = p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1))
    _nom3 = {'Hind': r'LZ~129 Hindenburg', 'FW': 'Flying Whales LCA60T', 'ZR': 'ZR-3 Los Angeles', 'ZNT': r'Zeppelin~NT'}
    for nom, d in p.PLATAFORMAS.items():
        x, fila = _m1[nom], _t3[_nom3[d['corto']]]
        # 2026-10-15: entra la columna P_rej = dm* q (segunda), y las de después se corren una posición
        chequeo('Z', f"{d['corto']}: tabla 3, P_rej = dm* q (q de EXP-06 a phi objetivo)", _num(fila[1]),
                x['dm_opt'] * r6['q_MJkg'], tol_abs=0.51, unidad='MW')
        chequeo('Z', f"{d['corto']}: tabla 3, escalón procesado = t_anticip/alpha* = (1+s) t_anticip", _num(fila[3]),
                x['t_anticipo'] / x['alpha_star'], tol_abs=0.051, unidad='s')
        chequeo('Z', f"{d['corto']}: tabla 3, velocidad al soltar = 3 D_lim/t_anticip", _num(fila[4]),
                3.0 * p.D_lim_oper / x['t_anticipo'], tol_abs=0.00051, unidad='m/s')
        # ANTES (hasta 2026-10-08): la tabla 3 traía la columna "Tanks return (s)" y se leía fila[6]; el autor la
        # quitó y el dato quedó solo para el LCA60T en la leyenda (valAbortoS), que se chequea abajo.
        # chequeo('Z', f"{d['corto']}: tabla 3, devolución si falla = (L/(rho_a/rho_H2)) t_anticip", _num(fila[6]),
        #         _mmE.L_eff / _kr3 * x['t_anticipo'], tol_abs=0.051, unidad='s')
    _fwf = _t3[_nom3['FW']]
    _tnat = next(x for x in rutas.leer_datos('freno')['filas'] if x['corto'] == 'FW' and x['m_kg'] == 1000.0)['t_nat_s']
    chequeo('Z', 'tabla 3, LCA60T: escalón procesado = el del freno (FrenoTnatFW)', _num(_fwf[3]), _tnat, tol_abs=0.051,
            unidad='s')
    _ropE = _mmE.resolver(next(d for d in p.PLATAFORMAS.values() if d['corto'] == 'FW'), 1000.0, p.h_max_precision,
                          p.D_lim_oper)
    # leyenda de la tabla 3 (valAbortoS): la devolución si falla la liberación, por dos caminos, el json del minimax
    # de las plataformas y el resolver (que es el que usa exportar_numeros.py)
    _nomFW = next(n for n, d in p.PLATAFORMAS.items() if d['corto'] == 'FW')
    chequeo('Z', 'tabla 3 (leyenda), LCA60T: devolución si falla = AbortoS (pesadez / (rho_a/rho_H2) / dm*)',
            _mmE.L_eff / _kr3 * _m1[_nomFW]['t_anticipo'],
            _mmE.L_eff * _ropE['dm_star'] * _ropE['t_anticip'] / _kr3 / _ropE['dm_star'], tol_abs=0.051, unidad='s')

    _t9 = _celdas('art_helio_rev2')
    _r12E = rutas.leer_datos('exp12_helio')
    _calor = [_num(c) for c in _t9['Heat to reject']]
    _pot = [_num(c) for c in _t9[r'Power, at once (14~s)']] if r'Power, at once (14~s)' in _t9 else None
    if _pot is None:
        _pot = [_num(c) for k_, v_ in _t9.items() if k_.startswith('Power, at once') for c in v_]
    for i_, (q_, P_) in enumerate(zip(_calor, _pot)):
        chequeo('I', f'tabla 9, columna {i_ + 1}: potencia de golpe = calor / tiempo de la maniobra limpia', P_,
                q_ * 1e3 / _tnat, tol_rel=0.03, tol_abs=0.6, unidad='MW')
    _ven = _t9['Energy in this stage'][3]
    chequeo('I', 'tabla 9: energía del H2 venteado = kg venteados x energía de electrólisis', _num(_ven),
            _num(_t9['Gas processed'][3]) * p.E_H2_electrolisis / 1e3, tol_abs=0.051, unidad='MWh')

    # críticas externas (autor, 2026-10-11): el freno a 1 atm sin compresor y el lazo con la inercia medida
    _c99v, _c1v = _r12E['rutas'][0], _r12E['rutas'][1]
    _tr1v = 1000.0 / _c1v['L'] * _c1v['q_Jkg'] / p.P_ref_masa
    chequeo('I', 'recuperación de 1 t con 10 MW a 1 atm: (1 t / L) q / P = calor por tonelada / P', _tr1v,
            _c1v['calor_GJ_por_t'] * 1e9 / p.P_ref_masa, tol_rel=1e-6, unidad='s')
    condicion('I', 'a 1 atm, sin compresor, la recuperación tarda más y se pierde más H2 que con phi = 0.99',
              _tr1v > 1000.0 / _c99v['L'] * _c99v['q_Jkg'] / p.P_ref_masa
              and _c1v['perdido_kg_por_t'] > _c99v['perdido_kg_por_t'], f'{_tr1v / 60:.1f} min a 1 atm')
    _fma = rutas.leer_datos('val01_componentes')['masa_anadida']['medida_sobre_analitica']
    condicion('I', 'con la masa vertical medida (Yamasaki & Goto) el lazo sigue sobreamortiguado: zeta / raiz(f) > 1',
              p.zeta / math.sqrt(max(_fma.values())) > 1.0, f'zeta = {p.zeta}, f = {_fma}')
    # compresor del escape escalado con el paquete del SGT-800 (fuentes.SGT800): trabajo isentrópico de una etapa
    from exp05_condensacion_tiempo_real import N_INERTE as _NIV, AIRE_ESTEQ as _AEV, p_sat as _psV
    _r6V = rutas.leer_datos('exp06_requisitos_minimos')
    _nV = _NIV / p.M_H2 * p.P_atmos / (p.P_atmos - _psV(p.T_salida_cond))      # inerte + vapor saturado [mol/kg H2]
    _CV = (_AEV - p.ratio_O2_H2) * p.cp_N2 + (_nV - _NIV / p.M_H2) * p.M_H2O * p.cp_vapor
    _kV, _T1V, _PRV = _nV * p.R_univ / _CV, p.T_salida_cond + 273.15, _r6V['P_req_atm']
    _wsV = _CV * _T1V * (_PRV ** _kV - 1.0)
    _NPV = 20000                                                    # integral de n R T(P) dP / P sobre la isentrópica
    _intV = sum(_nV * p.R_univ * _T1V * (1.0 + (_PRV - 1.0) * (i + 0.5) / _NPV) ** (_kV - 1.0) * (_PRV - 1.0) / _NPV
                for i in range(_NPV))
    chequeo('I', 'compresor del escape: trabajo isentrópico cerrado = integral de V dP sobre la isentrópica', _wsV,
            _intV, tol_rel=1e-6, unidad='J/kg')
    condicion('I', 'el trabajo isentrópico de una etapa supera el mínimo isotérmico del modelo (EXP-06)',
              _wsV > _r6V['W_comp_MJkg'] * 1e6, f'{_wsV / 1e6:.1f} > {_r6V["W_comp_MJkg"]:.1f} MJ/kg')
    _s8V, _kaV = p.SGT800, (p.gamma_aire - 1.0) / p.gamma_aire
    _wS8V = [(me - hr * 1e3 * P * 1e3 / 3600.0 / p.PCI_CH4) * p.cp_aire * (p.T_ISO_C + 273.15) * (pr ** _kaV - 1.0)
             for P, pr, me, hr in zip(_s8V['P_MWe'], _s8V['PR'], _s8V['m_escape_kgps'], _s8V['heat_rate_kJkWh'])]
    _dm10V = _c99v['kg_por_t'] / (1000.0 / _c99v['L'] * _c99v['q_Jkg'] / p.P_ref_masa)
    _mxV = _s8V['masa_kg'] * _dm10V * _wsV / min(_wS8V)                  # el más pesado de los cuatro [kg]
    _haV = rutas.leer_datos('helices_vectoriales')['referencia']['rechazo']['acc_haces_t'] * 1e3
    condicion('I', 'compresor del escape a 10 MW con el doble de kg/kW del SGT-800 < haces de ACC de 10 MW',
              2.0 * _mxV < _haV, f'{2 * _mxV / 1e3:.1f} t < {_haV / 1e3:.0f} t')
    _cuV = p.PLATAFORMAS[p.nombre_dirigible]['carga_max']
    condicion('I', 'el compresor de la cota ideal escalado igual no cabe en la carga útil del caso base',
              _s8V['masa_kg'] * _r6V['base_1t_h1']['dm_star_kgps'] * _wsV / max(_wS8V) > _cuV,
              f'carga útil {_cuV / 1e3:.0f} t')
    # el motor a 10 MW, ya descontados los ventiladores, entrega más que el pico de las hélices (sin batería)
    from balance_masa import f_caida_presion as _fcpV
    _WcV = rutas.leer_datos('exp06_requisitos_minimos')['W_comp_MJkg'] * 1e6     # compresor del escape (14 atm)
    _netV = min(eta * p.Q_comb_H2 * _dm10V - fv * (_c99v['q_Jkg'] - eta * p.Q_comb_H2 + _WcV) * _dm10V - _WcV * _dm10V
                for eta in p.eta_motor_H2 for fv in (min(_fcpV()), max(p.perdida_rad)))
    _helV = max(max(x['P_hel_W']) for x in rutas.leer_datos('freno')['filas']
                if x['corto'] == 'FW' and x['m_kg'] == 1000.0)
    condicion('I', 'freno a 10 MW: la potencia neta del motor (menor) supera el pico de las hélices (mayor), 1 t',
              _netV > _helV, f'{_netV / 1e6:.1f} MW > {_helV / 1e6:.2f} MW')
    # EXP-18: el LCA60T de la figura de alturas es la misma simulación que exp09 (su figura salió del paper)
    _e18V = rutas.leer_datos('exp18_por_nave')['agua']['FW']
    _e09V = {x['caso']: x for x in rutas.leer_datos('exp09_condensacion_dinamica')['casos_base']}
    for _k18, _k09 in (('10', 'ventear_10MW'), ('100', 'ventear_100MW'), ('inf', 'ventear_infMW')):
        chequeo('I', f'EXP-18, LCA60T, {_k18} MW: la subida = la de exp09 ({_k09}), el caso base',
                _e18V[_k18]['dh_arriba_m'], _e09V[_k09]['dh_arriba_m'], tol_abs=1e-6, unidad='m')
    # carta de diseño (EXP-18, parte 4): el punto de 1 t coincide con la barra HDM del desglose (fig_presupuesto_masa).
    # (2026-10-13) el compresor de H2 en unidades de catálogo: el techo es la suelta única más grande que cabe, la misma
    # que en el reparto (EXP-20) y en la carga contra el tiempo (EXP-22); con 1 t, una unidad (b = b' + una unidad)
    _ddV = rutas.leer_datos('exp18_por_nave')['diseno']
    _apMaxV = {nv['corto']: nv['max'] for nv in rutas.leer_datos('fig_presupuesto_masa')['apilado']['naves']}
    _r20V = rutas.leer_datos('exp20_reparto_sueltas')
    _cuV, _bpV = _ddV['compresor'], _ddV['b_sin_comp_t_por_t']
    _muV, _capV = _cuV['kg'] / 1e3, _cuV['cap_kgph'] / _cuV['kgph_por_t']
    _uniV = lambda m: max(1, int(np.ceil(m / _capV - 1e-9)))
    condicion('I', 'compresor de H2: una unidad de catálogo alcanza para regenerar 1 t en t_regeneracion_h',
              _capV >= 1.0, f'{_capV:.2f} t por unidad')
    chequeo('I', "compresor de H2: b = b' + una unidad (las fijas de 1 t, tabla 4)", _ddV['b_t_por_t'],
            _bpV + _muV, tol_rel=1e-12, unidad='t/t')
    for _c, _n in _ddV['naves'].items():
        _m = _n['limite_t']
        condicion('I', f"carta de diseño, {_c}: el techo cabe en la carga útil con sus unidades de compresor",
                  _m * (1.0 + _bpV) + _uniV(_m) * _muV <= _n['carga_t'] * (1.0 + 1e-12),
                  f"{_m * (1.0 + _bpV) + _uniV(_m) * _muV:.3f} <= {_n['carga_t']:.3f} t")
        condicion('I', f"carta de diseño, {_c}: el techo es el más grande (llena la carga útil o la capacidad de sus "
                  'unidades)', abs(_m * (1.0 + _bpV) + _uniV(_m) * _muV - _n['carga_t']) < 1e-9
                  or abs(_m - _uniV(_m) * _capV) < 1e-9, f'{_m:.3f} t, {_uniV(_m)} unidades')
        chequeo('I', f'carta de diseño, {_c}: techo = la suelta única más grande del reparto (EXP-20)', _m,
                _r20V['naves'][_c]['unica_max_t'], tol_rel=1e-12, unidad='t')
        if _apMaxV[_c] is None:
            condicion('I', f'carta de diseño, {_c}: la curva no llega a 1 t, como en el desglose', _n['t_1t_s'] is None,
                      f"techo {_n['limite_t']:.2f} t")
        else:
            chequeo('I', f'carta de diseño, {_c}: el punto de 1 t = la recuperación de la barra HDM', _n['t_1t_s'] / 60.0,
                    _apMaxV[_c]['t_rec_min'], tol_rel=1e-9, unidad='min')
    # EXP-20, reparto en sueltas (2026-10-12): con 1 t el compresor es una unidad, así que la recuperación más corta de
    # una suelta de 1 t es el punto de 1 t de la carta de diseño
    for _c, _n in _r20V['naves'].items():
        for _f in _n['sueltas']:
            _nu = round(_f['comp_t'] / _muV)
            condicion('I', f"reparto, {_c}, sueltas de {_f['s_t']:g} t: las unidades del compresor alcanzan y son las "
                      'mínimas', _nu * _capV >= _f['s_t'] * (1.0 - 1e-12) and (_nu == 1 or (_nu - 1) * _capV < _f['s_t']),
                      f'{_nu} unidades')
        _f1 = next((f for f in _n['sueltas'] if abs(f['s_t'] - 1.0) < 1e-9), None)
        if _f1 is not None and _ddV['naves'][_c]['t_1t_s'] is not None:
            chequeo('I', f'reparto, {_c}: la suelta de 1 t más rápida = el punto de 1 t de la carta de diseño',
                    _f1['t_min_s'], _ddV['naves'][_c]['t_1t_s'], tol_rel=1e-9, unidad='s')
    _fwR = {f['s_t']: f for f in _r20V['naves']['FW']['sueltas']}
    chequeo('I', 'reparto, LCA60T 1 t: la potencia de las hélices al soltar = k_min Q s / t_max (borde pesimista)',
            _fwR[1.0]['P_hel_pes_W'], _r20V['k'][0] * _r20V['Q_J_por_t'] * 1.0 / _fwR[1.0]['t_hel_pes_s'],
            tol_rel=1e-12, unidad='W')
    condicion('I', 'reparto, LCA60T: con sueltas de 1 t de 10 min las hélices alcanzan aun en el caso pesimista',
              _fwR[1.0]['t_hel_pes_s'] > 600.0, f"{_fwR[1.0]['t_hel_pes_s'] / 60:.0f} min")
    # el neto del motor de EXP-20 (k Q por tonelada) = el de EXP-21 sin las hélices (compresor del escape incluido)
    _k21 = rutas.leer_datos('exp21_balance_energia')['filas'][1]['kWh']
    chequeo('I', 'reparto: neto del motor por tonelada (k_min Q) = EXP-21, trabajo - ventiladores - compresor (pesimista)',
            _r20V['k'][0] * _r20V['Q_J_por_t'] / 3.6e6, _k21['trabajo_motor'][0] - max(p.perdida_rad)
            * _k21['calor_rechazado'][0] - _k21['compresor_escape'], tol_rel=1e-6, unidad='kWh')
    # (2026-10-14, autor) la cota pesimista de la entrega: el sistema ACC entero en vez de sus haces. a sube en la
    # diferencia de las dos masas específicas de la tabla 5; la ec. 12 en T reproduce la carga por suelta; con él, el
    # LZ 129 no alcanza con sueltas de 1 t de 30 min pero sí con 1 h (el texto de 5.3)
    _accV = _r20V['acc_entero']
    _candV = {c['clave']: c for c in rutas.leer_datos('req_masa_rechazo')['candidatos']}
    _apV = rutas.leer_datos('fig_presupuesto_masa')['apilado']
    chequeo('I', 'reparto con ACC entero: a - a(haces) = 1/(kW/kg del sistema) - 1/(kW/kg de los haces)',
            _accV['a_t_por_MW'] - _r20V['a_t_por_MW'],
            1.0 / _candV['ACC_sistema']['kWkg'][0] - 1.0 / _candV['ACC_haces']['kWkg'][0], tol_rel=1e-9, unidad='t/MW')
    chequeo('I', 'reparto con ACC entero: a Q = a x calor por tonelada', _accV['aQ_s'],
            _accV['a_t_por_MW'] * _r20V['Q_J_por_t'] / 1e6, tol_rel=1e-12, unidad='s')
    _wV, _b1V = _r20V['w_t_por_t'], _r20V['b_sin_comp_t_por_t']
    _Cfw = _r20V['naves']['FW']['carga_t']
    _Cp1V = _Cfw - _fwR[1.0]['comp_t'] - _b1V
    _e22V =rutas.leer_datos('exp22_energia_potencia')['E']
    _T6V = _e22V['n'] * _e22V['t_paso_s']
    _N6V = _Cp1V * _T6V / ((1.0 + _wV) * _T6V + _accV['aQ_s'])
    chequeo('I', 'reparto con ACC entero, LCA60T: la ec. 12 en la duración de la entrega = la carga por suelta con '
            't_rec = T / n', _N6V, (_Cp1V - _accV['aQ_s'] * _N6V / _T6V) / (1.0 + _wV), tol_rel=1e-9, unidad='t')
    condicion('I', 'reparto con ACC entero, LCA60T: sueltas de 1 t de 10 min no caben (la más rápida tarda más)',
              _accV['naves']['FW']['t_min_s'] > 600.0, f"{_accV['naves']['FW']['t_min_s'] / 60:.1f} min")
    _thV = _accV['naves']['Hind']['t_min_s']
    condicion('I', 'reparto con ACC entero, LZ 129: sueltas de 1 t de 30 min no caben y de 1 h sí',
              1800.0 < _thV < 3600.0 and _accV['naves']['Hind']['N_t']['3600'] >= 1.0, f'{_thV / 60:.1f} min')
    chequeo('I', 'reparto con ACC entero, LCA60T: la suelta de 1 t más rápida = la barra del ACC entero (desglose)',
            _accV['naves']['FW']['t_min_s'] / 60.0,
            next(n for n in _apV['naves'] if n['corto'] == 'FW')['max_acc']['t_rec_min'], tol_rel=1e-9, unidad='min')
    # EXP-22, la entrega de referencia de la sección 5.3 (2026-10-12): el reparto de EXP-20 con el balance por tonelada
    # de EXP-21; el H2 de los gases, recalculado aquí con el equilibrio de la nave antes y después de la entrega
    import minimax as _mmV
    _eV = rutas.leer_datos('exp22_energia_potencia')['E']
    _m21 = rutas.leer_datos('exp21_balance_energia')['filas'][1]['masa_kg']
    chequeo('I', 'entrega de referencia: equipo + carga de negocio = carga útil del LCA60T', _eV['equipo_t'] + _eV['N_t'],
            p.PLATAFORMAS[_eV['nave']]['carga_max'] / 1e3, tol_rel=1e-12, unidad='t')
    chequeo('I', 'entrega de referencia: H2 quemado = sueltas x H2 por tonelada (EXP-21)', _eV['total']['H2_kg'],
            _eV['n'] * _eV['s_t'] * _m21['H2_quemado'], tol_rel=1e-12, unidad='kg')
    _cfV = _eV['carga0_t'] - _eV['n'] * (_eV['s_t'] - _eV['s_t'] * _m21['agua_retenida'] / 1e3)
    chequeo('I', 'entrega de referencia: el H2 de los gases baja lo quemado (equilibrio, el agua queda como lastre)',
            _mmV.equilibrio(_eV['nave'], _eV['carga0_t'] * 1e3)['m_gas'] - _mmV.equilibrio(_eV['nave'], _cfV * 1e3)['m_gas'],
            _eV['total']['H2_kg'], tol_rel=1e-9, unidad='kg')
    chequeo('I', 'entrega de referencia: regeneración = sueltas x (electrólisis + compresión) por tonelada (EXP-21)',
            _eV['total']['regen_kWh'][1], _eV['n'] * _eV['s_t'] * (_k21['electrolisis'] + _k21['compresion_H2'][1]),
            tol_rel=1e-12, unidad='kWh')
    chequeo('I', 'entrega de referencia: días con el PEM del equipo = sueltas x t_regeneracion_h',
            _eV['t_regen_pem_h'], _eV['n'] * _eV['s_t'] * p.t_regeneracion_h, tol_rel=1e-12, unidad='h')
    condicion('I', 'entrega de referencia: las hélices al soltar caben en la potencia neta del motor (pesimista)',
              _eV['paso']['P_hel_W'][1] < _eV['paso']['motor_kWh'][0] * 3.6e6 / _eV['t_paso_s'],
              f"{_eV['paso']['P_hel_W'][1] / 1e6:.2f} < {_eV['paso']['motor_kWh'][0] * 3.6e6 / _eV['t_paso_s'] / 1e6:.2f} MW")
    # EXP-22 en la sección 6 (2026-10-13): el piso del costo por tonelada recuperada son los tanques del H2 de una
    # tonelada (balance_masa), y con sol tarda más de una semana en costar menos de lo que recupera; con la red, menos
    _e22V = rutas.leer_datos('exp22_energia_potencia')
    # objeciones del revisor (2026-10-14): 1.9, la tasa con phi_objetivo; 3.1, comprimir H2 en una entrega; 3.21, el
    # compuesto de los tanques casi no depende de la presión
    import minimax as _mmP
    from exp02_quema_centrada import L_ef as _LP
    _r12P = rutas.leer_datos('exp12_helio')
    chequeo('I', 'objeción 1.9: L con phi = 1 (minimax) / L con phi_objetivo (exp12) = la razón con L_ef de EXP-02',
            _mmP.L_eff / _r12P['rutas'][0]['L'], _LP(1.0) / _LP(p.phi_objetivo), tol_rel=1e-9, unidad='')
    _vmaP = rutas.leer_datos('val01_componentes')['masa_anadida']['efecto_flujo']
    condicion('I', 'objeción 1.9: la tasa con phi_objetivo sube menos que la incertidumbre de la inercia (vuelo libre)',
              _LP(1.0) / _LP(p.phi_objetivo) - 1.0 < 1.0 - _vmaP['vuelo_libre'],
              f"{100 * (_LP(1.0) / _LP(p.phi_objetivo) - 1.0):.2f} % < {100 * (1.0 - _vmaP['vuelo_libre']):.0f} %")
    _compE = _r12P['h2_comprimido']['reserva_compuesto_t'] * _e22V['E']['n'] * _e22V['E']['s_t']
    chequeo('I', 'objeción 3.1: compuesto por kg de H2 a 35 MPa, el mismo en el H2 comprimido y en los tanques del ciclo '
            '(densidad de Leachman y de NIST)', _r12P['h2_comprimido']['reserva_compuesto_t']
            / _r12P['h2_comprimido']['kg_por_t'], _r12P['rutas'][0]['reserva_compuesto_t'] / _r12P['rutas'][0]['kg_por_t'],
            tol_rel=1e-4, unidad='t/kg')
    condicion('I', 'objeción 3.1: comprimir H2 pide más compuesto que todo el equipo del ciclo para la misma entrega',
              _compE > _e22V['E']['equipo_t'], f"{_compE:.1f} t > {_e22V['E']['equipo_t']:.1f} t")
    # objeción 3.18: el O2 venteado al regenerar = agua retenida - H2 regenerado (masa del agua), y los tanques no
    # alcanzan a compensar con gas comprimido toda la sustentación que devuelve
    _bmP = rutas.leer_datos('balance_masa')
    _tqP = next(c for c in _bmP['cajas'] if c['caja'] == 'tanques_por_t')
    for _fP in rutas.leer_datos('exp21_balance_energia')['filas'][1:]:
        _mP = _fP['masa_kg']
        chequeo('I', f"objeción 3.18, {_fP['nombre']}: O2 venteado = agua retenida - H2 regenerado = H2 regenerado x r",
                _mP['agua_retenida'] - _mP['H2_regenerado'], _mP['H2_regenerado'] * p.ratio_O2_H2, tol_rel=1e-9,
                unidad='kg')
        _capP = _tqP['max_kg'] / _bmP['tanque_por_kg_H2']                       # [kg de H2 por t]
        _kP = p.ratio_O2_H2 / _r12P['h2_comprimido']['L']                       # gas a comprimir por kg regenerado
        condicion('I', f"objeción 3.18, {_fP['nombre']}: el regenerado y el gas que compensa su O2 no caben juntos en "
                  'los tanques (en vuelo no se regenera todo con solo los tanques como compensación)',
                  _mP['H2_regenerado'] < _capP < _mP['H2_regenerado'] * (1.0 + _kP),
                  f"{_mP['H2_regenerado']:.1f} < {_capP:.1f} < {_mP['H2_regenerado'] * (1.0 + _kP):.1f} kg por t")
    # el piso de la ec. 13 con tanques también para el gas que compensa el O2 (14 atm): el de EXP-22 x (1 + r / L)
    _pisoCP = _e22V['I4']['piso_t_por_t'] * (1.0 + p.ratio_O2_H2 / _r12P['h2_comprimido']['L'])
    condicion('I', 'objeción 3.18: con tanques también para el gas que compensa el O2, el piso supera la tonelada que '
              'recupera', _pisoCP > 1.0, f'{_pisoCP:.2f} t por t')
    _almP = [a['tanque_compuesto_t_por_t'] for a in _r12P['almacen']]
    condicion('I', 'objeción 3.21: el compuesto del He casi no depende de la presión (5-30 MPa, menos de 15 %)',
              max(_almP) / min(_almP) < 1.15, f'{min(_almP):.3f}-{max(_almP):.3f} t/t')
    chequeo('I', 'regeneración: piso del costo por tonelada = tanque por kg (balance_masa) x H2 regenerado por t',
            _e22V['I4']['piso_t_por_t'], rutas.leer_datos('balance_masa')['tanque_por_kg_H2']
            * _m21['H2_regenerado'] / 1e3, tol_rel=1e-12, unidad='t/t')
    condicion('I', 'regeneración: cuesta menos de 1 t por t recuperada antes de una semana con la red y después con a-Si',
              max(_e22V['I4']['fuentes']['red']['T_equilibrio_d']) < 7.0 < min(_e22V['I4']['fuentes']['a-Si']['T_equilibrio_d']),
              f"red {max(_e22V['I4']['fuentes']['red']['T_equilibrio_d']):.1f} d, a-Si "
              f"{min(_e22V['I4']['fuentes']['a-Si']['T_equilibrio_d']):.0f} d")
    condicion('I', 'entrega de referencia: aun optimista, el motor en vuelo paga menos de un cuarto de la regeneración',
              _eV['total']['neto_kWh'][1] < 0.25 * _eV['total']['regen_kWh'][0],
              f"{_eV['total']['neto_kWh'][1] / _eV['total']['regen_kWh'][0]:.2f}")
    # (2026-10-13) los límites de las hélices de la sección 5.3 salen de EXP-22 a 14 atm (lo que el texto afirma)
    _lsV = _e22V['limites_sueltas_s']
    condicion('I', 'hélices (EXP-22, 14 atm): con sueltas de 1 t de 10 min alcanzan aun en el caso pesimista',
              _lsV['1']['potencia_pes'] > 600.0, f"{_lsV['1']['potencia_pes'] / 60:.0f} min")
    condicion('I', 'hélices (EXP-22, 14 atm): con sueltas de 10 t no alcanzan en el pesimista y sí en el optimista',
              _lsV['10']['potencia_pes'] < _fwR[10.0]['t_min_s'] <= _lsV['10']['potencia_opt'],
              f"{_lsV['10']['potencia_pes'] / 60:.0f} < {_fwR[10.0]['t_min_s'] / 60:.0f} <= "
              f"{_lsV['10']['potencia_opt'] / 60:.0f} min")
    condicion('I', 'hélices (EXP-22, 14 atm): con sueltas de 20 t no alcanzan ni en el optimista (amarre)',
              _lsV['20']['potencia_opt'] < _fwR[20.0]['t_min_s'],
              f"{_lsV['20']['potencia_opt'] / 3600:.1f} < {_fwR[20.0]['t_min_s'] / 3600:.1f} h")
    _N1h = [_fwR[s]['N_t']['3600'] for s in sorted(_fwR)]
    condicion('I', 'reparto, LCA60T a 1 h: la carga de negocio baja con el tamaño de la suelta',
              all(y < x for x, y in zip(_N1h, _N1h[1:])), ', '.join(f'{x:.1f}' for x in _N1h))
    # (D, 2026-10-13) con el tanque de agua, que va con toda la carga entregada: N/C = (1 - f (b + aQ/t)) / (1 + b_w)
    _wV = _r20V['w_t_por_t']
    for _fr in (0.01, 0.02, 0.05, 0.1, 0.2):
        _ts = _fr * _r20V['aQ_s'] / (1.0 - _fr * (1.0 + _r20V['b_lineal_t_por_t'] + _wV))
        chequeo('I', f'reparto universal, sueltas del {100 * _fr:g} %: la curva nace sobre la suelta única (ec. 10)',
                (1.0 - _fr * (_r20V['b_lineal_t_por_t'] + _r20V['aQ_s'] / _ts)) / (1.0 + _wV), _fr, tol_rel=1e-12,
                unidad='')
    # el tanque de agua (D): la entrega cumple N (1 + b_w) = C - m_c - s (b' + a Q / t) en toda nave y toda suelta
    for _c, _nv in _r20V['naves'].items():
        for _f in _nv['sueltas']:
            for _x, _N in _f['N_t'].items():
                if _N >= _f['s_t']:
                    chequeo('I', f"reparto {_c}, {_f['s_t']:g} t cada {float(_x) / 60:.0f} min: N (1 + b_w) + equipo = C",
                            _N * (1.0 + _wV) + _f['comp_t'] + _f['s_t'] * (_r20V['b_sin_comp_t_por_t']
                                                                         + _r20V['aQ_s'] / float(_x)),
                            _nv['carga_t'], tol_rel=1e-12, unidad='t')
    # B2: el margen de diseño, recalculado con la ecuación del reparto y los datos de las cajas (independiente de EXP-20)
    _mgV = _r20V['margen']
    _CfV = p.PLATAFORMAS['Flying Whales LCA60T']['carga_max'] / 1e3
    for _nom, _fm, _fq in (('N_t', 1.0, 1.0), ('N_masa_t', 1.1, 1.0), ('N_calor_t', 1.0, 1.1)):
        chequeo('I', f'B2, margen ({_nom}): N (1 + f_m b_w) = C - f_m m_c - s (f_m b\' + f_m f_q a Q / t), 1 t de 10 min',
                _mgV[_nom], (_CfV - _fm * _r20V['comp_t'] - 1.0 * (_fm * _r20V['b_sin_comp_t_por_t']
                                                                  + _fm * _fq * _r20V['aQ_s'] / 600.0))
                / (1.0 + _fm * _wV), tol_rel=1e-12, unidad='t')
    condicion('I', 'B2: con 10 % de margen de masa o de calor, la entrega de 1 t de 10 min sigue llevando carga',
              min(_mgV['N_masa_t'], _mgV['N_calor_t']) > 1.0,
              f"{_mgV['N_masa_t']:.1f} t (masa), {_mgV['N_calor_t']:.1f} t (calor)")
    _taV = rutas.leer_datos('balance_masa')['tanque_agua']
    chequeo('I', 'tanque de agua: el agua que guarda por tonelada = la que retiene cada suelta de la entrega (EXP-22)',
            _taV['agua_por_t_kg'], _eV['paso']['agua_kg'] / _eV['s_t'], tol_rel=1e-6, unidad='kg')
    chequeo('I', 'tanque de agua: b_w de EXP-20 = la clase liviana de balance_masa (tela de 30 oz/yd2)', _wV * 1e3,
            _taV['kg_por_t'][0], tol_rel=1e-12, unidad='kg/t')
    chequeo('I', 'tanque de agua: tela de las dos caras de 14 x 14 ft x gramaje / agua a la capacidad (independiente)',
            _taV['kg_por_kg_agua'][0], 2 * (14 * 0.3048) ** 2 * 30 * 0.028349523125 / 0.83612736
            / _taV['agua_cap_kg'], tol_rel=1e-12, unidad='kg/kg')
    condicion('I', 'tanque de agua: la almohada llena (12.5 x 12.5 x 4 ft, MIL-T-53029C fig. 1) contiene su capacidad',
              (12.5 * 0.3048) ** 2 * 4 * 0.3048 > 3000 * 3.785411784e-3,
              f"{(12.5 * 0.3048) ** 2 * 4 * 0.3048:.1f} m3 de caja > {3000 * 3.785411784e-3:.1f} m3")
    from exp05_condensacion_tiempo_real import rho_liq_sat as _rls
    for _T, _rho in ((6.96957, 999.857), (47.6831, 989.026), (61.1165, 982.583), (99.8843, 958.432)):
        chequeo('L', f'densidad del agua líquida saturada a {_T:g} °C (Wagner & Pruss) = NIST WebBook (IAPWS-95)',
                _rls(_T), _rho, tol_rel=2e-5, unidad='kg/m3')
    chequeo('I', 'la mitad de la entrega: a Q / (1 + b_w) redondea igual que a Q (las dos en 0.1 h)',
            round(_r20V['aQ_s'] / (1.0 + _wV) / 3600, 1), round(_r20V['aQ_s'] / 3600, 1), tol_abs=1e-9, unidad='h')
    # la figura del desglose con el NT soltando 0.2 t (2026-10-12): la barra llena la carga útil. (2026-10-13) con el
    # compresor en unidades lleva una unidad entera, así que ya no recupera como la ZR-3 con 1 t
    _apfV = {nv['corto']: nv for nv in rutas.leer_datos('fig_presupuesto_masa')['apilado_fig']['naves']}
    _ap1V = {nv['corto']: nv for nv in rutas.leer_datos('fig_presupuesto_masa')['apilado']['naves']}
    chequeo('I', 'desglose, NT con 0.2 t: la barra HDM suma la carga útil', sum(_apfV['ZNT']['max']['t']),
            _apfV['ZNT']['carga_t'], tol_rel=1e-9, unidad='t')
    chequeo('I', "desglose, NT con 0.2 t: las fijas = b' x 0.2 t + una unidad del compresor", _apfV['ZNT']['max']['t'][3],
            _bpV * _apfV['ZNT']['masa_t'] + _muV, tol_rel=1e-12, unidad='t')
    for _c in ('Hind', 'FW', 'ZR'):
        chequeo('I', f'desglose, {_c}: la figura usa la misma barra de 1 t que las macros', sum(_apfV[_c]['max']['t']),
                sum(_ap1V[_c]['max']['t']), tol_rel=1e-12, unidad='t')
    # presupuesto apilado de la figura de la sección 5.2 (fig_presupuesto_masa)
    _apV = rutas.leer_datos('fig_presupuesto_masa')['apilado']
    _rqV = rutas.leer_datos('req_masa_rechazo')
    _kwV = {c['clave']: c['kWkg'][0] for c in _rqV['candidatos']}
    _fwV = next(nv for nv in _apV['naves'] if nv['corto'] == 'FW')
    chequeo('I', 'presupuesto apilado: el compresor del escape sin freno, LCA60T = el de la cota ideal escalado (macro '
            'CompMasaBase)', _fwV['sin']['t'][1], _s8V['masa_kg'] * _r6V['base_1t_h1']['dm_star_kgps'] * _wsV
            / max(_wS8V) / 1e3, tol_rel=1e-9, unidad='t')
    chequeo('I', 'presupuesto apilado: el rechazo con 10 MW = 10 MW / kW/kg de los haces de ACC',
            _fwV['con']['t'][0], p.P_ref_masa / 1e3 / _kwV['ACC_haces'] / 1e3, tol_rel=1e-9, unidad='t')
    chequeo('I', 'presupuesto apilado: el motor con 10 MW = t por kg/s del más liviano x la tasa del freno',
            _fwV['con']['t'][2], _rqV['requisitos']['conversion_t_por_kgps'][0] * _dm10V, tol_rel=1e-6, unidad='t')
    for nv in _apV['naves']:
        if nv['max'] is not None:
            chequeo('I', f"presupuesto apilado, {nv['nombre']}: la barra del rechazo más grande que cabe suma la carga "
                    'útil', sum(nv['max']['t']), nv['carga_t'], tol_rel=1e-9, unidad='t')
        else:
            condicion('I', f"presupuesto apilado, {nv['nombre']}: no cabe nada porque las fijas y la tonelada ya superan "
                      'la carga útil', _apV['fijas_t'] + 1.0 > nv['carga_t'],
                      f"{_apV['fijas_t'] + 1.0:.2f} t > {nv['carga_t']:.2f} t")
    # tercera barra (autor, 2026-10-12): HDM a 1 atm, sin el compresor del escape, dentro de la misma carga útil
    for nv in _apV['naves']:
        _x1 = nv['max_una_atm']
        if _x1 is None:
            condicion('I', f"presupuesto apilado a 1 atm, {nv['nombre']}: tampoco cabe nada (el techo no depende de la "
                      'presión)', nv['max'] is None, 'sin barra HDM a 14 atm ni a 1 atm')
            continue
        chequeo('I', f"presupuesto apilado a 1 atm, {nv['nombre']}: la barra suma la carga útil", sum(_x1['t']),
                nv['carga_t'], tol_rel=1e-9, unidad='t')
        condicion('I', f"presupuesto apilado a 1 atm, {nv['nombre']}: sin compresor del escape", _x1['t'][1] == 0.0,
                  f"{_x1['t'][1]} t")
        chequeo('I', f"presupuesto apilado a 1 atm, {nv['nombre']}: la recuperación = el punto de 1 t de la curva de 1 atm "
                'de la carta de diseño', _x1['t_rec_min'], _ddV['naves'][nv['corto']]['t_1t_s_una_atm'] / 60.0,
                tol_rel=1e-9, unidad='min')
    # carta de diseño a 1 atm, sin el compresor del escape (panel derecho, comparación; la base sigue a 14 atm)
    _uV = _ddV['una_atm']
    _rt12V = rutas.leer_datos('exp12_helio')['rutas']
    if '1 atm' not in _rt12V[1]['ruta']:
        sys.exit('rutas de exp12: la segunda ya no es la de 1 atm')
    chequeo('I', 'carta de diseño a 1 atm: Q / 10 MW = la recuperación de 1 t a 1 atm de la sección 5.2 (ruta de 1 atm de '
            'la tabla 7)', _uV['Q_MJ_por_t'] * 1e6 / p.P_ref_masa, 1000.0 / _rt12V[1]['L'] * _rt12V[1]['q_Jkg']
            / p.P_ref_masa, tol_rel=1e-9, unidad='s')
    chequeo('I', 'carta de diseño a 1 atm: a = los haces por MW y el motor de la barra de 10 MW (sin el compresor), el '
            'motor a la tasa de q(0, phi_max)', _uV['a_t_por_MW'], (_fwV['con']['t'][0] + _fwV['con']['t'][2]
            * _rt12V[0]['q_Jkg'] / _rt12V[1]['q_Jkg']) / (p.P_ref_masa / 1e6), tol_rel=1e-9, unidad='t/MW')
    for _c, _n in _ddV['naves'].items():
        if _n['t_1t_s_una_atm'] is None:
            condicion('I', f'carta de diseño a 1 atm, {_c}: tampoco llega a 1 t (el techo no depende de la presión)',
                      _n['t_1t_s'] is None, f"techo {_n['limite_t']:.2f} t")
            continue
        chequeo('I', f'carta de diseño a 1 atm, {_c}: en su punto de 1 t cabe justo 1 t', _n['carga_t'] / (
            1.0 + _ddV['b_t_por_t'] + _uV['aQ_s'] / _n['t_1t_s_una_atm']), 1.0, tol_rel=1e-12, unidad='t')
        chequeo('I', f'carta de diseño a 1 atm, {_c}: t(1 atm) / t(14 atm) = a Q (1 atm) / a Q (14 atm), igual en toda nave',
                _n['t_1t_s_una_atm'] / _n['t_1t_s'], _uV['razon_t'], tol_rel=1e-12, unidad='')
    # (2026-10-13) la entrega a 1 atm (auditoría, 5.6): a Q de EXP-20 = la de la carta de diseño a 1 atm
    chequeo('I', 'entrega a 1 atm: a Q de EXP-20 = la de la carta de diseño a 1 atm (EXP-18)', _r20V['una_atm']['aQ_s'],
            _uV['aQ_s'], tol_rel=1e-12, unidad='s')
    condicion('I', 'entrega a 1 atm: a Q menor que a 14 atm, así que en el mismo tiempo entrega más',
              _r20V['una_atm']['aQ_s'] < _r20V['aQ_s'],
              f"{_r20V['una_atm']['aQ_s'] / 3600:.2f} h < {_r20V['aQ_s'] / 3600:.2f} h")

    _fe = rutas.leer_datos('fig_esferoides')
    for nom, d in p.PLATAFORMAS.items():
        n_ = _fe['naves'][d['corto']]
        chequeo('W', f"{d['corto']}: figura, volumen del esferoide = pi l D^2 / 6", n_['V_esferoide_m3'],
                math.pi * d['L'] * d['D'] ** 2 / 6.0, tol_rel=1e-12, unidad='m3')
        chequeo('W', f"{d['corto']}: figura, k del punto = k de la tabla 2 (Lamb en su esbeltez)", n_['k'], d['k'],
                tol_abs=5e-4)
    chequeo('W', 'figura: el esferoide discontinuo del LCA60T encierra su V_casco',
            math.pi / 6.0 * _fe['naves']['FW']['L_equiv_m'] * _fe['naves']['FW']['D_equiv_m'] ** 2,
            _fe['naves']['FW']['V_casco_m3'], tol_rel=1e-9, unidad='m3')

    # ── Salidas ──────────────────────────────────────────────────────────────
    n_ok = sum(c['ok'] for c in CH)
    # conservación del agua en las simulaciones del artículo (grupo E)
    cons = [c['err_rel'] for c in CH if c['grupo'] == 'E' and 'retenida + venteada' in c['nombre']]
    # qué tipo de chequeo es cada uno (M-6, auditoría 3): descriptivo, para no presentar el total como si fueran todos
    # comparaciones independientes; los de diferencia nula incluyen identidades y regresiones de la misma fórmula
    n_cond = sum(1 for c in CH if c['obtenido'] != c['obtenido'])
    n_nulos = sum(1 for c in CH if c['obtenido'] == c['obtenido'] and c['err_rel'] < 1e-9)
    rutas.guardar_datos('verificacion', {'n_total': len(CH), 'n_ok': n_ok, 'n_fallas': len(CH) - n_ok,
                                         'n_condiciones': n_cond, 'n_diferencia_nula': n_nulos,
                                         'n_con_tolerancia': len(CH) - n_cond - n_nulos,
                                         'n_pdf_sin_comprobar': sum(1 for c in CH if c.get('sin_comprobar')),
                                         'sello_datos': rutas.sello_datos(),
                                         'err_rel_max_conservacion': max(cons), 'chequeos': CH})
    tabla()
    for c in CH:
        est = 'ok   ' if c['ok'] else 'FALLA'
        det = (f'obtenido {c["obtenido"]:.6g} vs esperado {c["esperado"]:.6g} (tol {c["tolerancia"]})'
               if c['obtenido'] == c['obtenido'] else c['tolerancia'])
        print(f'  {est} [{c["grupo"]}] {c["nombre"]}: {det}')
    print(f'\n{n_ok}/{len(CH)} chequeos ok ({n_cond} condiciones, {n_nulos} de diferencia nula, '
          f'{len(CH) - n_cond - n_nulos} con tolerancia; PDF sin comprobar: '
          f"{sum(1 for c in CH if c.get('sin_comprobar'))})")
    if n_ok < len(CH):
        sys.exit(1)


def tabla():
    def num(x):
        if x != x:
            return '--'
        return f'{x:.4g}'
    enc = r'Grupo & Chequeo & Obtenido & Esperado & Tolerancia & Estado \\'
    L = [r'\begin{longtable}{c p{7.2cm} r r l c}', r'\toprule', enc, r'\midrule', r'\endfirsthead',
         r'\toprule', enc, r'\midrule', r'\endhead', r'\bottomrule', r'\endfoot']
    previo = None
    for c in CH:
        if previo is not None and c['grupo'] != previo:
            L.append(r'\addlinespace')
        previo = c['grupo']
        nombre = c['nombre'].replace('_', r'\_').replace('%', r'\%')
        tol = c['tolerancia'].replace('%', r'\%')
        L.append(f"{c['grupo']} & {nombre} & {num(c['obtenido'])} & {num(c['esperado'])} & {tol} & "
                 f"{'ok' if c['ok'] else r'\textbf{FALLA}'}" + r' \\')
    L += [r'\end{longtable}']
    with open(rutas.tabla('verificacion'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/verificacion.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
