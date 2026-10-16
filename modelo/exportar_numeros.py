# -*- coding: utf-8 -*-
# [EN] Writes resultados/numeros.tex, the LaTeX macros with every number of the article.
# =============================================================================
# exportar_numeros.py — arma resultados/numeros.tex a partir de resultados/datos/
#
# Cada cifra del modelo que aparece en el artículo se define UNA vez
# aquí, con su redondeo. El .tex escribe \valAlphaStar en vez de 0.497: "cita"
# el resultado igual que cita la bibliografía. Si una constante cambia, basta
# correr correr_todo.py y compilar.
#
# Reglas:
#   - nombres de macro: prefijo "val" + solo letras (LaTeX no admite dígitos)
#   - sufijos por plataforma: Hind, FW, ZR, ZNT (parametros.PLATAFORMAS[...]['corto']); el caso base
#     (parametros.nombre_dirigible) lleva además el sufijo Base (alias_base, al final)
#   - separador de miles "{,}" solo desde 10 000 (p. ej. 107{,}481)
# =============================================================================

import re
import sys
import io
import math
import parametros as p
import rutas

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
BASE = p.NAVE_BASE['corto']     # sufijo corto del caso base (FW desde 2026-09-28)

MACROS = []            # (nombre, texto, comentario)


def fmt(x, dec):
    """Redondeo único: dec decimales; "{,}" como separador de miles desde 10 000."""
    if x != x:                      # NaN: sin solución verificada
        return '--'
    s = f'{x:,.{dec}f}' if abs(x) >= 10000 else f'{x:.{dec}f}'
    return s.replace(',', '{,}')


def val(nombre, x, dec, comentario):
    nombre = 'val' + nombre
    if not re.fullmatch(r'[A-Za-z]+', nombre):
        sys.exit(f'Nombre de macro inválido (solo letras): {nombre}')
    if any(m[0] == nombre for m in MACROS):
        sys.exit(f'Macro repetida: {nombre}')
    MACROS.append((nombre, fmt(x, dec), comentario))


def rango(nombre, a, b, dec, comentario):
    """Macro con un rango 'a--b' (p. ej. datos 2026 -- metas 2035)."""
    nombre = 'val' + nombre
    if not re.fullmatch(r'[A-Za-z]+', nombre):
        sys.exit(f'Nombre de macro inválido (solo letras): {nombre}')
    if any(m[0] == nombre for m in MACROS):
        sys.exit(f'Macro repetida: {nombre}')
    ta, tb = fmt(a, dec), fmt(b, dec)
    MACROS.append((nombre, ta if ta == tb else f'{ta}--{tb}', comentario))   # rango degenerado: un solo número


def seccion(titulo):
    MACROS.append((None, None, titulo))


def alias_base(pref, suf=''):
    """La macro por nave del caso base, también con sufijo Base: el texto cita \\val<pref>Base<suf> y no
    depende de qué nave es el caso base (hasta 2026-09-28 el sufijo era PF)."""
    fuente, nombre = 'val' + pref + BASE + suf, 'val' + pref + 'Base' + suf
    if any(m[0] == nombre for m in MACROS):
        sys.exit(f'Macro repetida: {nombre}')
    MACROS.append((nombre, next(m[1] for m in MACROS if m[0] == fuente), f'= \\{fuente} (caso base)'))


def mismo(a, b, que, tol=1e-9):
    """Chequeo cruzado: la misma cifra calculada por dos scripts debe coincidir."""
    if abs(a - b) > tol * max(1.0, abs(a)):
        sys.exit(f'INCONSISTENCIA en {que}: {a} vs {b}')


# ── Parámetros de entrada (Tabla 1) ──────────────────────────────────────────
seccion('Parámetros del modelo (parametros.py) — Tabla 1')
val('Zeta',              p.zeta, 2,                 'coef. de amortiguamiento del PD')
# la ganancia usa rho_aire; la rigidez real por volumen de gas es (rho_aire - rho_H2) g: el H2 quemado sale con su
# masa (hallazgo M14). zeta real = zeta raiz(1 - rho_H2/rho_aire)
# ANTES (hasta 2026-10-09): ZetaReal = zeta raiz(1 - rho_H2/rho_aire), porque K_d usaba rho_aire g; ahora K_d usa la
# rigidez real y el amortiguamiento real es zeta mismo.
val('ZetaReal', p.zeta, 2, 'amortiguamiento real del lazo (= zeta: K_d con la rigidez real)')
# rigidez por altura del modelo (declarada, autor 2026-09-29): la sustentación sigue a rho_aire(h) con V_gas fijo,
# rho_aire V_gas g / H_escala, frente a la del lazo (rho_aire - rho_H2) g K_p, en el caso base con 1 t
import minimax as _mmr
_eqr = _mmr.equilibrio(p.NAVE_BASE, p.m_carga)
_k_alt = p.rho_aire_std * _eqr['V_gas'] * p.g / p.H_escala
_k_lazo = (p.rho_aire_std - p.P_atmos / (p.R_H2 * p.T1)) * p.g * p.Kp
val('RigAlturaPct', 100 * _k_alt / _k_lazo, 0, 'rigidez por altura del modelo / rigidez del lazo, caso base [%]')
val('Kp',                p.Kp, 0,                   'ganancia proporcional [m3/m]')
val('KmEmergencia',      p.K_mass_control, 0,       'K_m fase emergencia [kg/s/m3]')
val('KmTrim',            p.K_mass_trim, 1,          'K_m fase trim [kg/s/m3]')
val('DmMaxEstandar',     p.max_flujo_H2, 0,         'capacidad del esquema estándar [kg/s]')
val('DmProactivoEstandar', abs(p.dm_H2_proactivo), 0, 'pulso proactivo estándar [kg/s]')
val('TAnticipEstandar',  p.t_anticip_estandar, 1,   'ventana del pulso estándar [s]')
val('Dlim',              p.D_lim_oper, 2,           'D_lim operacional [m]')
val('HmaxPrec',          p.h_max_precision, 0,      'umbral de precisión [m]')
val('HmaxSeg',           p.h_max_seguridad, 0,      'umbral de seguridad [m]')
val('Hcero',             p.h0, 0,                   'altitud inicial [m]')
val('Cd',                p.Cd, 1,                   'coeficiente de arrastre')
val('EtaGen',            p.rend_gen, 2,             'eficiencia del generador')
val('EElectrolisis',     p.E_H2_electrolisis, 0,    'energía de electrólisis [kWh/kg]')
val('EtaPanel',          p.eta_panel, 2,            'eficiencia fotovoltaica')
val('FPaneles',          p.A_paneles_fraccion, 2,   'fracción de área con paneles')
val('TCondObjetivo',     p.t_cond_objetivo, 0,      'ventana de condensación [min]')
val('IntervaloMultidrop', p.intervalo_multidrop, 0, 'tiempo entre entregas [s]')
val('TAireExterior',     p.T_aire_exterior, 0,      'aire exterior del caso de diseño del condensador [°C]')
val('TRefComb',          p.T_ref_comb, 0,           'temperatura de referencia del PCI y de L_v [°C]')
val('FRadHoy',           p.perdida_rad[0], 2,       'potencia de ventilación por kW rechazado, radiador 2026 (FlyZero)')
val('FRadFuturo',        p.perdida_rad[1], 2,       'potencia de ventilación por kW rechazado, radiador 2035 (FlyZero)')
rango('KWkgRad',         p.kWkg_radiador[0], p.kWkg_radiador[1], 0, 'radiador de aire, kW rechazados por kg, 2026--2035 (FlyZero)')
rango('KWkgTMS',         p.kWkg_TMS_sistema[0], p.kWkg_TMS_sistema[1], 0,
      'sistema de gestión térmica completo, kW rechazados por kg de sistema, 2026--2035 (FlyZero p. 4; caja negra)')
# ANTES (revisión 2026-09-30, I-04/N-4: con un decimal, 0.7-1.3, un lector obtiene 63 kW/kg y 40 t en vez de 56 y 43):
# rango('KWkgBat',         p.kWkg_bateria[0], p.kWkg_bateria[1], 1, 'potencia específica del buffer a nivel paquete [kW/kg]')
rango('KWkgBat',         p.kWkg_bateria[0], p.kWkg_bateria[1], 2, 'potencia específica del buffer a nivel paquete [kW/kg]')
rango('KWkgBatCelda',    p.kWkg_bateria_celda[0], p.kWkg_bateria_celda[1], 0, 'potencia específica de celda Li-ion [kW/kg]')
val('FCeldaPaquete',     100 * p.f_celda_paquete, 0, 'masa de celda / masa de paquete, NASA X-57 [%]')
val('WhkgBatCelda',      p.Whkg_bateria_celda, 0, 'energía específica de celda Li-ion [Wh/kg]')
val('DTBufferAgua',      p.dT_buffer_agua, 0,       'calentamiento del agua si se guardara el calor en vez de rechazarlo [K]')
val('WhkgBat',           p.Whkg_bateria, 0,         'energía específica de la batería, a nivel de paquete [Wh/kg]')

# ── Quema limpia y D_lim (apendice_dlim_flujo.py) ────────────────────────────
d = rutas.leer_datos('dlim_flujo')
seccion('Quema limpia a D_lim operacional (dlim_flujo.json)')
val('AlphaStar', d['alpha_star'], 1, 'alpha* de la ley general (s = 1): 1/2, igual para todas')
for nom, x in d['plataformas'].items():
    c = x['corto']
    val('MStar' + c,      x['dm_star_kgps'], 2,            f'ṁ* quema limpia [kg/s] — {nom}')
    val('MTotal' + c,     x['M_total_kg'], 0,              f'M_total con 1 t de carga [kg] — {nom}')
    val('HundTotal' + c,  x['hundimiento_total_m'], 2,     f'bajada total (1+s) D_lim [m] — {nom}')
    val('MCota' + c,      x['dm_cota_reinyeccion_kgps'], 2, f'cota con reinyección (quema máxima hasta el pico, 1 m) [kg/s] — {nom}')
for _pref in ('MStar', 'MTotal', 'HundTotal', 'MCota'):
    alias_base(_pref)
_pfd = d['plataformas'][p.nombre_dirigible]
val('CotaBajaPct', 100 * (1 - _pfd['dm_cota_reinyeccion_kgps'] / _pfd['dm_star_kgps']), 0,
    'cuánto baja el flujo instalado la cota con reinyección, caso base [%]')
val('CotaSobrequema', 100 * _pfd['sobrequema_cota'], 0, 'H2 quemado de más en la cota con reinyección, caso base [%]')

# ── Fuente de O2 y regeneración solar (comparacion_O2_corrected.py) ──────────
o = rutas.leer_datos('fuente_O2')
seccion('Fuente de O2 y regeneración solar (fuente_O2.json)')
val('RhoAireHdos',   o['rho_aire_sobre_rho_H2'], 1, 'rho_aire / rho_H2')
val('Leff',          o['L_eff'], 1,                 'L_eff = rho_aire/rho_H2 + r_O2/H2 (O2 del aire que queda a bordo en el agua)')
val('HdosTanques',   o['m_H2_tanques_kg'], 1,       'H2 por entrega con O2 de tanques [kg]')
val('HdosDescenso',  o['m_H2_aire_kg'], 1,          'H2 por entrega con O2 del aire [kg]')
val('ReduccionOdos', o['reduccion_pct'], 1,         'reducción de H2 aire vs tanques [%]')
val('AhorroHdos',    o['ahorro_kg'], 1,             'H2 ahorrado por entrega [kg]')
val('AguaTanques',   o['agua_tanques_kg'], 0,       'agua por entrega, O2 de tanques [kg]')
val('AguaDescenso',  o['agua_aire_kg'], 0,          'agua por entrega, O2 del aire [kg]')
val('EgenTanques',   o['E_gen_tanques_kWh'], 0,     'energía recuperable, O2 de tanques [kWh]')
val('EgenDescenso',  o['E_gen_aire_kWh'], 0,        'energía recuperable por entrega [kWh]')
for nom, x in o['regeneracion_solar'].items():
    c = x['corto']
    val('RegenKgDia' + c, x['regen_kg_dia'], 0, f'regeneración solar [kg H2/día] — {nom}')
    val('RegenDias' + c,  x['dias'], 2,         f'días para regenerar 1 entrega — {nom}')
    val('RegenHoras' + c, x['horas'], 1,        f'horas para regenerar 1 entrega — {nom}')

# ── Condensación (condensacion.py) ───────────────────────────────────────────
k = rutas.leer_datos('condensacion')
mismo(k['m_H2_descenso_kg'], o['m_H2_aire_kg'], 'H2 de descenso (condensacion vs fuente_O2)')
seccion('Condensación (condensacion.json)')
val('QcondMin', k['Q_cond_min_kW'], 0, 'potencia mínima de condensación [kW]')

# ── Esquema estándar (t2_gain_scheduled.py) ──────────────────────────────────
e = rutas.leer_datos('estandar')
seccion('Esquema estándar, pulso fijo (estandar.json)')
val('EgenEstandar',      e['base_1t']['E_kWh'], 0,     'energía, 1 entrega [kWh]')
val('DhEstandar',        e['base_1t']['dh_max_m'], 2,  '|Δh| máx, 1 entrega [m]')
val('EgenEstandarMulti', e['multi_4x1t']['E_kWh'], 0,  'energía, 4x1 t [kWh]')
val('DhEstandarMulti',   e['multi_4x1t']['dh_max_m'], 2, '|Δh| máx, 4x1 t [m]')

# ── Minimax por plataforma (drop_flujo_optimo_perdrop.py) ────────────────────
m = rutas.leer_datos('minimax_plataformas')
b = rutas.leer_datos('minimax_base')
pf = m['single_1t'][p.nombre_dirigible]
mismo(b['t_anticipo'], pf['t_anticipo'], 't_anticip del caso base (minimax_base vs plataformas)')
mismo(b['dm_opt'], pf['dm_opt'], 'ṁ* del caso base (minimax_base vs plataformas)')
import minimax as _mmx
_rx = _mmx.resolver(p.nombre_dirigible, 1000.0, p.h_max_precision, p.D_lim_oper)
mismo(pf['alpha_star'], _rx['alpha'], 'alpha* EXACTO del caso base (simulación vs minimax.resolver)')
mismo(pf['dm_opt'], _rx['dm_star'], 'ṁ* EXACTO del caso base (simulación vs minimax.resolver)')
mismo(d['alpha_star'], _rx['alpha_aprox'], 'alpha* de la LEY GENERAL (dlim_flujo vs minimax.resolver)')
seccion('Quema limpia por plataforma, 1 entrega de 1 t (minimax_plataformas.json) — Tabla 3')
_ax = [x['alpha_star'] for x in m['single_1t'].values()]
rango('AlphaStarExacto', min(_ax), max(_ax), 3, 'alpha* exacto, rango de las 4 plataformas')
_dm = [(x['dm_opt'], _mmx.resolver(n, 1000.0, p.h_max_precision, p.D_lim_oper)['dm_star_aprox'])
       for n, x in m['single_1t'].items()]
_dif = [abs(e / a_ - 1.0) for e, a_ in _dm]
_gr = [abs(e / a_ - 1.0) for (e, a_), n in zip(_dm, m['single_1t']) if p.PLATAFORMAS[n]['corto'] != 'ZNT']
val('ExactoDifGrandes', 100 * max(_gr), 1, 'diferencia exacto vs ley general, naves grandes [%]')
val('ExactoDifMax', 100 * max(_dif), 0, 'diferencia exacto vs ley general, peor caso (Zeppelin NT) [%]')
for nom, x in m['single_1t'].items():
    c = x['corto']
    mismo(x['dm_opt'], d['plataformas'][nom]['dm_star_kgps'], f'ṁ* de {nom} (simulación vs analítico)')
    val('TAnticip' + c, x['t_anticipo'], 2,  f't_anticip [s] — {nom}')
    val('Egen' + c,     x['Energia_kWh'], 0, f'energía generada [kWh] — {nom}')
    val('DhPico' + c,   x['dh_peak_m'], 2,   f'|Δh| pico [m] — {nom}')
    val('HdosNeto' + c, x['m_H2_neto'], 1,   f'H2 neto consumido [kg] — {nom}')
for _pref in ('TAnticip', 'Egen', 'DhPico', 'HdosNeto'):
    alias_base(_pref)
seccion('Minimax multi-drop, 4 entregas (Zeppelin NT: 4 x 250 kg)')
for nom, x in m['multidrop'].items():
    c = x['corto']
    val('Ecum' + c,          x['E_cum_kWh'], 0, f'energía acumulada [kWh] — {nom}')
    val('HdosNetoMulti' + c, x['m_H2_neto'], 1, f'H2 neto [kg] — {nom}')
    val('DhPicoMulti' + c,   x['dh_peak_m'], 2, f'|Δh| pico [m] — {nom}')

# ── Barrido reactivo (t3_barrido_reactivo.py) ────────────────────────────────
r = rutas.leer_datos('barrido_reactivo')
seccion('Barrido reactivo puro (barrido_reactivo.json) — Tabla 3')
for nom, x in r.items():
    c = x['corto']
    val('MMin' + c,       x['dm_min_kgps'], 0,   f'ṁ_min (|Δh| <= 5 m) [kg/s] — {nom}')
    val('MFine' + c,      x['dm_fine_kgps'], 0,  f'ṁ_fine (|Δh| <= 1 m) [kg/s] — {nom}')
    val('DhMin' + c,      x['dh_en_min_m'], 2,   f'|Δh| a ṁ_min [m] — {nom}')
    val('DhFine' + c,     x['dh_en_fine_m'], 3,  f'|Δh| a ṁ_fine [m] — {nom}')
    val('TEstabFine' + c, x['t_estab_fine_s'], 1, f't_estab a ṁ_fine [s] — {nom}')

# ── Mapa de requisitos verificado (mapa_requisitos.py, Paso 1) ───────────────
q = rutas.leer_datos('mapa_requisitos')
seccion(f"Mapa de requisitos verificado por simulación (mapa_requisitos.json), h_max = {q['h_max_tabla_m']:g} m")
val('AirePorKgHdos',          q['aire_por_kg_H2'], 1,           'kg de aire estequiométrico por kg de H2')
val('SobrepasoDlimHmaxMin',   q['sobrepaso_D_igual_h_min'], 1,  'excursión / h_max con D_lim = h_max (mín.)')
val('SobrepasoDlimHmaxMax',   q['sobrepaso_D_igual_h_max'], 1,  'excursión / h_max con D_lim = h_max (máx.)')
val('MapaConsInsuficiente',   q['conservador_forma_cerrada_insuficiente'], 0,
    'casos en que la forma cerrada (D_lim 0.10 m) se queda corta')
val('MapaConsFactorMax',      q['conservador_factor_max'], 2,   'ṁ verificado / forma cerrada, máx. (D_lim 0.10 m)')
for nom, x in q['en_capacidad_h5'].items():
    c = x['corto']
    val('MapaCarga' + c,     x['carga_t'], 0 if x['carga_t'] == int(x['carga_t']) else 1,
        f'carga evaluada (su capacidad, máx. 30 t) [t] — {nom}')
    val('MapaMCons' + c,     x['cons_dm_verificado_kgps'], 1,  f'ṁ verificado, D_lim 0.10 m [kg/s] — {nom}')
    val('MapaMMejor' + c,    x['mejor_dm_verificado_kgps'], 1, f'ṁ verificado, mejor D_lim [kg/s] — {nom}')
    val('MapaDlimMejor' + c, x['mejor_D_lim_m'], 2,            f'mejor D_lim verificado [m] — {nom}')
    val('MapaPth' + c,       x['mejor_potencia_termica_MW'], 0, f'potencia térmica, mejor caso [MW] — {nom}')
    val('MapaAire' + c,      x['mejor_aire_min_kgps'], 0,      f'aire mínimo, mejor caso [kg/s] — {nom}')
    val('MapaTEstab' + c,    x['mejor_t_estab_s'], 0,          f't_estab, mejor caso [s] — {nom}')

# ── EXP-02: quema constante centrada, sin reinyección (exp02_quema_centrada.py) ──
x = rutas.leer_datos('exp02_quema_centrada')
seccion('EXP-02: quema constante centrada, sin reinyección, de reposo a reposo (exp02_quema_centrada.json)')
val('KMasaAnadida',      x['k_masa_anadida'], 2,      'k: masa añadida vertical / aire desplazado')
val('LeffEscape',        x['L_phi0'], 1,              'L si el escape sale con el vapor (phi = 0)')
val('FactorMasaAnadida', x['factor_masa_anadida'], 2, 'factor del flujo con masa añadida, 1/sqrt(1+k)')
val('FactorEscape',      x['factor_escape'], 2,       'factor del flujo si el agua no queda a bordo')
q1, q2, q3 = x['base_1t_1m'], x['FW_20t_30m'], x['FW_20t_5m']
val('CentradaBaseIdeal',   q1['ideal_dm_kgps'], 2,            'ṁ quema centrada, caso base 1 t, H = 1 m, ideal [kg/s]')
val('CentradaBaseMa',      q1['con_masa_anadida_dm_kgps'], 2, 'ṁ quema centrada, caso base 1 t, H = 1 m, + masa añadida')
val('CentradaFWIdeal',   q2['ideal_dm_kgps'], 1,            'ṁ quema centrada, FW 20 t, H = 30 m, ideal [kg/s]')
val('CentradaFWMa',      q2['con_masa_anadida_dm_kgps'], 1, 'ṁ quema centrada, FW 20 t, H = 30 m, + masa añadida')
val('CentradaFWEscape',  q2['escape_sale_dm_kgps'], 1,      'ṁ quema centrada, FW 20 t, H = 30 m, escape sale')
val('CentradaFWArrastre', q2['con_arrastre_dm_kgps'], 1,   'ṁ quema centrada, FW 20 t, H = 30 m, + masa añadida + arrastre')
val('CentradaFWKappa',   q2['con_arrastre_kappa'], 2,     'kappa = c H / M_e (validez del primer orden), FW 20 t, 30 m')
val('CentradaFWAlphaReposo', q2['con_arrastre_alpha_reposo'], 3, 'alfa que termina en reposo con arrastre, FW 20 t, 30 m')
val('CentradaFWT',       q2['con_masa_anadida_T_s'], 0,     'duración de la quema, FW 20 t, H = 30 m [s]')
val('CentradaFWVmax',    q2['con_masa_anadida_v_max_ms'], 1, 'v máx. de bajada, FW 20 t, H = 30 m [m/s]')
val('CentradaFWCincoIdeal', q3['ideal_dm_kgps'], 0,         'ṁ quema centrada, FW 20 t, H = 5 m, ideal [kg/s]')

# ── EXP-05: condensación en tiempo real, brecha de diseño (exp05_condensacion_tiempo_real.py) ──
z = rutas.leer_datos('exp05_condensacion_tiempo_real')
mismo(z['articulo_Qcond_min_kW'], k['Q_cond_min_kW'], 'Q_cond,min del artículo (exp05 vs condensacion)')
mismo(z['L_phi0'], x['L_phi0'], 'L(phi = 0) (exp05 vs exp02)')
seccion('EXP-05: condensación en tiempo real como brecha de diseño (exp05_condensacion_tiempo_real.json)')
val('TRocioEscape',      z['T_rocio_C'], 0,          'punto de rocío del escape estequiométrico [°C]')
val('XHdosOEscape',      100 * z['x_H2O'], 0,        'fracción molar de H2O en el escape [%]')
val('TSalidaCond',       z['T_salida_C'], 0,         'temperatura de salida del condensador [°C]')
val('PhiMaxSalida',      z['phi_max_salida'], 2,     'fracción de agua que condensa a T_salida')
val('PhiMaxTuno',        z['phi_max_T1'], 2,         'fracción de agua que condensa a 15 °C')
val('LeffPhiSalida',     z['L_phi'], 1,              'L(phi) con phi = phi_max(T_salida)')
val('QRech',             z['q_rech_MJkg'], 0,        'calor a rechazar por kg H2, eta = 0 [MJ/kg]')
val('QRechGen',          z['q_rech_gen_MJkg'], 0,    'calor a rechazar por kg H2 si rend_gen sale como trabajo [MJ/kg]')
val('QLatente',          z['q_latente_MJkg'], 1,     'solo latente por kg H2 (lo que cuenta el artículo) [MJ/kg]')
val('FactorVsLatente',   z['factor_vs_latente'], 1,  'calor total / solo latente')
val('PrtArticuloBase',     z['articulo_base_P_rt_MW'], 0, 'potencia en tiempo real en el punto del artículo (caso base 1 t, ṁ*) [MW]')
pf5, fw5 = z['base_1t'], z['FW_20t']
val('QUnaTon',           pf5['Q_GJ'], 1,             'calor de una liberación de 1 t [GJ]')
val('PDiezMinUnaTon',    pf5['P_10min_MW'], 1,       'potencia para rechazarlo en 10 min, 1 t [MW]')
val('PCascoMinBase',       pf5['P_casco_min_MW'], 1,   'casco caso base como intercambiador, aire en calma [MW]')
val('PCascoMaxBase',       pf5['P_casco_max_MW'], 1,   'casco caso base como intercambiador, aire forzado [MW]')
val('TMinCascoMaxBase',    pf5['t_min_casco_max_min'], 0, 'tiempo mínimo para rechazar 1 t con el casco caso base (forzado) [min]')
val('TMinCascoMinBase',    pf5['t_min_casco_min_min'], 0, 'idem, aire en calma [min]')
val('HUnaLibCascoBase',    pf5['H_una_liberacion_casco_max_m'], 0, 'bajada de una liberación de 1 t limitada por el casco caso base [m]')
val('PrtCincoBase',        pf5['P_rt_H5_MW'], 0,       'potencia en tiempo real, caso base 1 t, H = 5 m [MW]')
val('QVeinteFW',         fw5['Q_GJ'], 0,             'calor de una liberación de 20 t [GJ]')
val('PCascoMaxFW',       fw5['P_casco_max_MW'], 1,   'casco FW como intercambiador, aire forzado [MW]')
val('PrtTreintaFW',      fw5['P_rt_H30_MW'], 0,      'potencia en tiempo real, FW 20 t, H = 30 m [MW]')
val('PrtCincoFW',        fw5['P_rt_H5_MW'], 0,       'potencia en tiempo real, FW 20 t, H = 5 m [MW]')
val('TMinCascoMaxFW',    fw5['t_min_casco_max_min'], 0, 'tiempo mínimo para rechazar 20 t con el casco FW (forzado) [min]')
val('PDiezMinFW',        fw5['P_10min_MW'], 0,       'potencia para rechazar 20 t en 10 min [MW]')
val('BrechaRtBase',        pf5['P_rt_H5_MW'] / pf5['P_casco_max_MW'], 0, 'brecha: tiempo real (H = 5 m) / casco forzado, caso base 1 t')
val('BrechaRtFW',        fw5['P_rt_H5_MW'] / fw5['P_casco_max_MW'], 0, 'brecha: tiempo real (H = 5 m) / casco forzado, FW 20 t')
val('BrechaDiezBase',      pf5['P_10min_MW'] / pf5['P_casco_max_MW'], 1, 'brecha: 10 min / casco forzado, caso base 1 t')
val('BrechaDiezFW',      fw5['P_10min_MW'] / fw5['P_casco_max_MW'], 0, 'brecha: 10 min / casco forzado, FW 20 t')
val('PasoFraccionado',   z['paso_fraccionado_kg'], 0, 'tamaño de paso de la liberación fraccionada [kg]')
val('NPasosFW',          fw5['n_500kg'], 0,          'pasos de 500 kg para 20 t')
val('PPasoFW',           fw5['P_500kg_MW'], 0,       'potencia por paso de 500 kg, FW 20 t, H total 5 m [MW]')
val('DmPasoFW',          fw5['dm_500kg_kgps'], 1,    'flujo por paso de 500 kg, FW 20 t, H total 5 m [kg/s]')
val('DuracionPasosFW',   fw5['duracion_500kg_min'], 0, 'duración total en pasos de 500 kg, FW 20 t [min]')
val('NMinCascoFW',       fw5['n_min_casco_max'], 0,  'pasos mínimos para caber en el casco FW (forzado)')
val('PasoMinCascoFW',    fw5['paso_min_kg'], 0,      'tamaño de paso correspondiente [kg]')
val('HdosUnaTonPhi',     pf5['H2_kg'], 1,            'H2 total por t de carga con phi = phi_max(T_salida) [kg]')
pr = {int(r['P_atm']): r for r in z['presion']}
for n_atm, nom in ((3, 'Tres'), (10, 'Diez')):
    val('TRocio' + nom + 'Atm', pr[n_atm]['T_rocio_C'], 0, f'rocío del escape a {n_atm} atm [°C]')
    val('Phi' + nom + 'Atm',    pr[n_atm]['phi_40C'], 2,   f'fracción condensada a 40 °C y {n_atm} atm')
    val('WIso' + nom + 'Atm',   pr[n_atm]['W_iso_MJkg'], 1, f'trabajo isotérmico mínimo para comprimir el escape a {n_atm} atm [MJ/kg H2]')
val('DTNdosLatente',     z['dT_N2_latente_K'], 0,    'calentamiento del N2 si absorbiera todo el latente [K]')

# ── EXP-06: requisitos mínimos, minimax + condensar todo (exp06_requisitos_minimos.py) ──
r6 = rutas.leer_datos('exp06_requisitos_minimos')
mismo(r6['base_1t_h1']['dm_star_kgps'], d['plataformas'][p.nombre_dirigible]['dm_star_kgps'], 'ṁ* exacto caso base 1 t (exp06 vs dlim_flujo)')
seccion('EXP-06: requisitos mínimos, minimax sin reinyección y condensando todo (exp06_requisitos_minimos.json)')
val('PhiObjetivo',       r6['phi_objetivo'], 2,     'fracción de agua a condensar ("todo")')
val('LeffObjetivo',      r6['L_phi'], 1,            'L con phi_objetivo')
val('QTodo',             r6['q_MJkg'], 0,           'calor a rechazar por kg H2 condensando todo [MJ/kg]')
val('PReqAtm',           r6['P_req_atm'], 0,        'presión del escape para condensar phi_objetivo a T_salida [atm]')
val('TUnaAtm',           r6['T_1atm_C'], 0,         'temperatura para condensar phi_objetivo a 1 atm [°C]')
val('WComp',             r6['W_comp_MJkg'], 1,      'trabajo mínimo (isotérmico) de compresión [MJ/kg H2]')
val('WCompPct',          100 * r6['W_comp_frac_PCI'], 0, 'trabajo de compresión como % del PCI')
# M-5 (auditoría 3): la cota ideal junta phi = 1 (quema) con el calor de condensar phi objetivo a 40 C, que a 1 atm
# pide 14 atm; a 1 atm (phi_max) el calor por tonelada sería mayor, y el compresor no está en la masa
from exp02_quema_centrada import L_ef as _Lm5
from exp05_condensacion_tiempo_real import q_rech as _qm5, phi_max as _phm5
_Qt = lambda ph: 1000.0 / _Lm5(ph) * _qm5(0.0, ph)
val('CalorUnaAtmPct', 100.0 * (_Qt(_phm5(p.T_salida_cond)) / _Qt(p.phi_objetivo) - 1.0), 0,
    'calor por tonelada a 1 atm (phi_max) frente al de phi objetivo [% más]')
# R-04 (revisión 2026-10-01): caso coherente a 1 atm (phi_max en la tasa y en el calor)
mismo(r6['base_1t_1atm']['phi'], _phm5(p.T_salida_cond), 'phi_max a 1 atm (exp06 vs exp05)')
val('ReqDmUnaAtm',    r6['base_1t_1atm']['dm_star_kgps'], 2, 'ṁ* caso base con phi_max a 1 atm en la tasa [kg/s]')
val('ReqPcondUnaAtm', r6['base_1t_1atm']['P_cond_MW'], 0,    'P_cond caso base coherente a 1 atm [MW]')
val('ReqLUnaAtm',     r6['base_1t_1atm']['L'], 2,            'L(phi_max) a 1 atm')
val('ReqQUnaAtm',     r6['base_1t_1atm']['q_MJkg'], 1,       'q(0, phi_max) a 1 atm [MJ/kg]')
val('CompMinPotBase', r6['base_1t_h1']['dm_star_kgps'] * r6['W_comp_MJkg'], 0,
    'potencia mínima del compresor del escape a la presión de phi objetivo, caso base [MW]')
val('AireEnfr',          p.dT_aire_enfr, 0,         'calentamiento del aire de enfriamiento [K]')
val('DTMediaCond',       p.dT_media_cond, 0,        'diferencia media del condensador [K]')
for clave, suf in (('base_1t_h1', 'BaseUno'), ('base_1t_h5', 'BaseCinco'), ('FW_20t_h1', 'FWUno'), ('FW_20t_h5', 'FWCinco')):
    x6 = r6[clave]
    val('ReqDm' + suf,     x6['dm_star_kgps'], 1 if x6['dm_star_kgps'] >= 100 else 2, f'ṁ* minimax [kg/s] — {clave}')
    val('ReqPth' + suf,    x6['P_termica_MW'], 0,  f'potencia térmica [MW] — {clave}')
    val('ReqAire' + suf,   x6['aire_comb_kgps'], 0, f'aire de combustión [kg/s] — {clave}')
    val('ReqPcond' + suf,  x6['P_cond_MW'], 0,     f'potencia de condensación en tiempo real [MW] — {clave}')
    val('ReqEnfr' + suf,   x6['aire_enfr_tps'], 1, f'aire de enfriamiento [t/s] — {clave}')
    val('ReqUA' + suf,     x6['UA_MWK'], 1,        f'UA del condensador [MW/K] — {clave}')
    val('ReqEsp' + suf,    x6['P_especifica_kWkg'], 0, f'potencia de condensación por kg de carga [kW/kg] — {clave}')
    val('ReqDmDh' + suf,   x6['dm_star_Dh3_kgps'], 1 if x6['dm_star_Dh3_kgps'] >= 100 else 2, f'ṁ* con D_lim = h/3 [kg/s] — {clave}')
    val('ReqPcondDh' + suf, x6['P_cond_Dh3_MW'], 0, f'P_cond con D_lim = h/3 [MW] — {clave}')
val('RatioEnfrComb',     r6['base_1t_h1']['ratio_aire_enfr_comb'], 0, 'aire de enfriamiento / aire de combustión')
val('FactorMaReq',       r6['factor_masa_anadida'], 2, 'factor de masa añadida sobre ṁ* y P_cond')


def paso6(nombre, h, P_mw):
    return next(x for x in r6['pasos'] if x['plataforma'] == nombre and x['h_m'] == h and x['P_instalada_MW'] == P_mw)


# escalón máximo con 10 MW en las tres naves grandes (antes el texto daba solo caso base-FW, 62-68 kg, y dejaba fuera
# al Hindenburg; auditoría R2-05)
_pg = [paso6(n_, p.h_max_seguridad, 10.0)['paso_max_kg'] for n_, d_ in p.PLATAFORMAS.items() if d_['corto'] != 'ZNT']
rango('PasoMaxGrandesDiez', min(_pg), max(_pg), 0, 'paso máximo con 10 MW, naves grandes (Hind, FW, ZR) [kg]')
for (nombre, P_mw), suf in (((p.nombre_dirigible, 10.0), 'BaseDiez'), ((p.nombre_dirigible, 100.0), 'BaseCien'),
                            (('Flying Whales LCA60T', 10.0), 'FWDiez'), (('Flying Whales LCA60T', 100.0), 'FWCien')):
    x6 = paso6(nombre, p.h_max_seguridad, P_mw)
    val('PasoMax' + suf,  x6['paso_max_kg'], 0, f'paso máximo, h = 5 m, {P_mw:g} MW instalados [kg] — {nombre}')
    val('NPasos' + suf,   x6['n_pasos'], 0,     f'número de pasos para la capacidad, {P_mw:g} MW — {nombre}')
    val('TEntrega' + suf, x6['t_total_min'], 0, f'duración de la entrega completa, {P_mw:g} MW [min] — {nombre}')
    val('TEntregaHoras' + suf, x6['t_total_min'] / 60.0, 1, f'duración de la entrega completa, {P_mw:g} MW [h] — {nombre}')
val('CargaTotalBase',      paso6(p.nombre_dirigible, p.h_max_seguridad, 10.0)['carga_t'], 0, 'carga total caso base (capacidad) [t]')
val('BufferAgua',        r6['buffer_agua_t_por_t'], 0, 'agua a bordo para absorber el calor de 1 t con +50 K [t]')
val('CargaTotalFW',      paso6('Flying Whales LCA60T', p.h_max_seguridad, 10.0)['carga_t'], 0, 'carga total FW evaluada [t]')

v2 = rutas.leer_datos('ver02_numerica')
seccion('VER-02: verificación numérica de la simulación del artículo (ver02_numerica.json)')
ex2 = {a['config']: a for a in v2['exacta']}
val('VerDosExactaZ', 100 * ex2['paso /2, rtol 1e-8']['err_rel_z'], 3, 'error vs solución exacta, altura, tolerancia estricta [%]')
val('VerDosExactaAnt', 100 * ex2['anterior (rtol 1e-3)']['err_rel_z'], 1, 'error vs solución exacta con la tolerancia anterior [%]')
un = [max(x['u_rel'].values()) for filas in v2['convergencia'].values() for x in filas if x['config'].startswith('adoptada')]
ua = [max(x['u_rel'].values()) for filas in v2['convergencia'].values() for x in filas if x['config'].startswith('anterior')]
def _cota_sup(x):
    """Cota superior a una cifra significativa (el texto dice "below ..."): 0.0132 -> (0.02, 2)."""
    e = math.floor(math.log10(x))
    return math.ceil(round(x / 10.0 ** e, 9)) * 10.0 ** e, max(0, -e)
val('VerDosUnum', *_cota_sup(100 * max(un)), 'cota superior (1 cifra) de la incertidumbre numérica de la simulación adoptada (máx. excursión, H2, agua) [%]')
val('VerDosUnumAnt', 100 * max(ua), 1, 'incertidumbre numérica con la tolerancia anterior [%]')

# ── VAL-03: cotas y sensibilidades de la PIRT (val03_cotas.py) ────────────────
v3 = rutas.leer_datos('val03_cotas')
seccion('VAL-03: cotas y sensibilidades de la PIRT (val03_cotas.json)')
t3 = v3['tau']
rango('TresTauFly', min(t3['tau_FlyZero_s']), max(t3['tau_FlyZero_s']), 1, 'tau = c_Al/k del radiador de FlyZero (2035-2026) [s]')
rango('TresTauGhatos', t3['tau_Ghatos_forzada_s'][0], t3['tau_Ghatos_forzada_s'][1], 0, 'tau medido, intercambiador de laboratorio, aire forzado (Ghatos 2023) [s]')
val('TresCAl', t3['c_Al_J_kgK'], 0, 'Cp del aluminio a 298 K (JANAF) [J/(kg K)]')
val('TresCdBase', 100 * max(abs(f['d_excursion_rel']) for f in v3['cd']['base']), 1, 'cambio máx. de la excursión con Cd 0-1, cota ideal [%]')
val('TresCdCond', 100 * max(abs(f['d_excursion_rel']) for f in v3['cd']['condensador_100MW_tau5']), 1, 'cambio máx. de la excursión con Cd 0-1, 100 MW y tau 5 s [%]')
te3 = v3['termica']
val('TresHint', te3['h_W_m2K'], 1, 'convección natural interna del H2 (Churchill & Chu) [W/(m2 K)]')
val('TresTauTerm', te3['tau_th_s'][1] / 60, 0, 'tiempo térmico del gas del caso base [min]')
rango('TresTauTermRango', te3['tau_th_s'][0] / 60, te3['tau_th_s'][2] / 60, 0, 'tiempo térmico con h x2 a x0.5 [min]')
val('TresDurMinimax', te3['dur_minimax_s'], 0, 'duración de la quema minimax (1-99 % del H2) [s]')
val('TresDurPasos', te3['dur_pasos_s'] / 60, 1, 'duración de 1 t en escalones con 10 MW [min]')
rango('TresDerivaMinimax', te3['deriva']['5K_minimax']['dm_kg'][1], te3['deriva']['20K_minimax']['dm_kg'][1], 0, 'deriva en la maniobra minimax, dT_eq 5-20 K [kg]')
rango('TresDerivaPasos', te3['deriva']['5K_pasos_10MW_1t']['dm_kg'][1] / 1e3, te3['deriva']['10K_pasos_10MW_1t']['dm_kg'][1] / 1e3, 1, 'deriva en 1 t por escalones a 10 MW, dT_eq 5-10 K [t]')
rango('TresDerivaEntrega', te3['deriva']['5K_entrega']['dm_kg'][1] / 1e3, te3['deriva']['10K_entrega']['dm_kg'][1] / 1e3, 1,
      'deriva en toda la entrega de referencia (LCA60T cargado al salir, sueltas una tras otra), dT_eq 5-10 K [t]')
val('TresTauTermEntrega', te3['entrega']['tau_th_s'][1] / 60, 0, 'tiempo térmico del gas del LCA60T cargado al salir de la entrega [min]')
rango('TresDTMinimax', te3['deriva']['5K_minimax']['dT_K'], te3['deriva']['20K_minimax']['dT_K'], 0, 'cambio de sobrecalentamiento solar supuesto, quema minimax [K]')
rango('TresDTPasos', te3['deriva']['5K_pasos_10MW_1t']['dT_K'], te3['deriva']['10K_pasos_10MW_1t']['dT_K'], 0, 'cambio de sobrecalentamiento solar supuesto, escalones [K]')
rango('TresCdRango', min(f['Cd'] for f in v3['cd']['condensador_100MW_tau5']), max(f['Cd'] for f in v3['cd']['condensador_100MW_tau5']), 0, 'rango de Cd de la sensibilidad')
val('TresTauCota', v3['tau']['tau_cota_s'], 0, 'cota alta del retardo del condensador [s]')
inc = v3['incertidumbre']
for esc, nom in (('ideal', 'Ideal'), ('real', 'Real')):
    for k, abrev, dec in (('dm_star', 'Dm', 2), ('P_MW', 'P', 0), ('H2_t', 'Hdos', 1), ('Q_GJ_t', 'Qt', 2)):
        m = inc[esc]['rango'][k]
        rango('TresInc' + nom + abrev, m['min'], m['max'], dec, f'incertidumbre {esc}: {k} mín.-máx. en los extremos medidos')
val('TresPhiMax', inc['phi_max'], 2, 'phi_max a la salida del condensador (escenario real)')

# ── CRED-01: evaluación de credibilidad, NASA-STD-7009B (credibilidad.py; revisión 2026-10-01) ──
cr1 = rutas.leer_datos('credibilidad')
seccion('CRED-01: evaluación de credibilidad con los factores de NASA-STD-7009B (credibilidad.json)')
for _f in cr1['capacidad'] + cr1['resultados']:
    val('CredNivel' + _f['clave'], _f['nivel'], 0, f"nivel {_f['factor']} (juicio de los autores, credibilidad.py)")
val('CredUmbral', cr1['umbral'], 0, 'umbral de nivel fijado por este estudio [M&S 43e]')
val('CredCumplenCap', sum(f['nivel'] >= cr1['umbral'] for f in cr1['capacidad']), 0, 'factores de capacidad que llegan al umbral')
val('CredCumplenRes', sum(f['nivel'] >= cr1['umbral'] for f in cr1['resultados']), 0, 'factores de resultados que llegan al umbral')
val('CredCumplen', sum(f['nivel'] >= cr1['umbral'] for f in cr1['capacidad'] + cr1['resultados']), 0, 'factores que llegan al umbral, de 11')
_ea = cr1['elasticidades_base']['analitica']
val('CredElasCarga', _ea['m'], 2, 'elasticidad de ṁ* respecto de la carga, caso base (ecuación 7)')
val('CredElasInercia', -_ea['M_I'], 2, 'elasticidad de ṁ* respecto de la inercia, caso base, EN VALOR ABSOLUTO (es negativa)')
_en = cr1['elasticidades_naves']
val('CredElasCargaMax', max(x['m'] for x in _en.values()), 2, 'elasticidad respecto de la carga, máx. entre naves (Zeppelin NT)')
val('CredElasInerciaMax', -min(x['M_I'] for x in _en.values()), 2, 'elasticidad respecto de la inercia, máx. |.| entre naves')
val('CredElasPhi', -cr1['elasticidad_phi']['ideal'], 2, '|d ln ṁ*/d ln phi| = r phi/L en phi = 1 (también del H2 por tonelada)')
val('CredElasPhiAtm', -cr1['elasticidad_phi']['atm'], 2, '|d ln ṁ*/d ln phi| en phi_max a 1 atm')
_fl = {x['clave']: x['requisitos'] for x in cr1['robustez']['filas']}
_b1 = cr1['robustez']['base']['dm_star']
rango('CredPhiDmPct', 100 * (min(r['dm_star'] for r in _fl['phi'][1:]) / _b1 - 1), 100 * (max(r['dm_star'] for r in _fl['phi'][1:]) / _b1 - 1), 0,
      'cuánto sube ṁ* con el agua retenida medida (phi_max x recolección de Zakhvatkin) [%]')
rango('CredPhiRet', min(cr1['robustez']['phi_max'] * x for x in cr1['robustez']['fphi']),
      max(cr1['robustez']['phi_max'] * x for x in cr1['robustez']['fphi']), 2, 'fracción del agua retenida con la recolección medida')
import tablas_articulo as _ta
_v1 = {x['magnitud']: x for x in rutas.leer_datos('val01_componentes')['filas']}
val('CredRefDesvMax', 100 * max(abs(_v1[c]['desvio_rel']) for _, c, _t, _f, _d in _ta.SEL_VALIDACION
                                if c != 'densidad del aire en h0'), 2,
    'desvío máx. de los datos de entrada frente a su referencia, sin la densidad del aire (tabla S3) [%]')
val('CredFuentesTotal', cr1['fuentes_paper1']['total'], 0, 'fuentes del paper 1 en el registro (fuentes.py)')
for _g in 'ABC':
    val('CredFuentes' + _g, cr1['fuentes_paper1'][_g], 0, f'fuentes del paper 1 con pedigrí {_g}')

# ── Control reactivo puro vs quema limpia (reactivo_vs_limpia.py; revisión 2026-10-01, sec. 4) ──
rvl = rutas.leer_datos('reactivo_vs_limpia')
seccion('Control reactivo puro frente a la quema limpia, mismo L (reactivo_vs_limpia.json)')
_rv = rvl['naves']
_rvb = _rv[p.nombre_dirigible]
_rvg = [x for x in _rv.values() if x['corto'] != 'ZNT']
val('ReacUnoBase', _rvb['reactivo_sim']['precision'], 1, 'reactivo, tasa para no subir más de 1 m, caso base [kg/s]')
val('ReacIgualBase', _rvb['reactivo_sim']['igual'], 0, 'reactivo, tasa para no pasar la excursión de la quema limpia, caso base [kg/s]')
rango('ReacRazonUno', min(x['razon_sim']['precision'] for x in _rv.values()),
      max(x['razon_sim']['precision'] for x in _rv.values()), 1, 'reactivo simulado / quema limpia, 1 m, 4 naves')
rango('ReacRazonIgual', min(x['razon_sim']['igual'] for x in _rv.values()),
      max(x['razon_sim']['igual'] for x in _rv.values()), 0, 'reactivo simulado / quema limpia, misma excursión, 4 naves')
rango('ReacRazonCinco', min(x['razon_sim']['seguridad'] for x in _rv.values()),
      max(x['razon_sim']['seguridad'] for x in _rv.values()), 2, 'reactivo simulado / quema limpia, 5 m, 4 naves')
val('ReacIdealUno', _rvb['razon_ideal']['precision'], 1, 'reactivo ideal / quema limpia, 1 m, caso base')
val('ReacIdealIgual', _rvb['razon_ideal']['igual'], 0, 'reactivo ideal / quema limpia, misma excursión, caso base')
val('ReacIdealCinco', _rvb['razon_ideal']['seguridad'], 1, 'reactivo ideal / quema limpia, 5 m, caso base')
rango('ReacMismoUno', min(x['razon_mismo']['precision'] for x in _rv.values()),
      max(x['razon_mismo']['precision'] for x in _rv.values()), 1,
      'reactivo / quema limpia con el MISMO desplazamiento, 1 m (la limpia con D_lim = 1 m / (1+s)), 4 naves')
rango('ReacMismoCinco', min(x['razon_mismo']['seguridad'] for x in _rv.values()),
      max(x['razon_mismo']['seguridad'] for x in _rv.values()), 1, 'lo mismo, 5 m, 4 naves')
rango('ReacIgualTodas', min(x['reactivo_sim']['igual'] for x in _rv.values()),
      max(x['reactivo_sim']['igual'] for x in _rv.values()), 0, 'reactivo, tasa para la misma excursión, 4 naves [kg/s]')

# ── Revisión del 2026-10-01: f físico (R-03), humedad (R-07), hélices existentes (discusión) ──
seccion('Revisión 2026-10-01: ventiladores por caída de presión, humedad del aire, hélices existentes')
_rq = rutas.leer_datos('req_masa_rechazo')['requisitos']
rango('DpRechazo', p.dp_rechazo_Pa[0], p.dp_rechazo_Pa[1], 0, 'caída de presión del lado del aire supuesta [Pa]')
val('EtaFan', p.eta_ventilador, 1, 'rendimiento supuesto del soplado de la unidad de rechazo')
rango('FFisico', _rq['f_fisico'][0], _rq['f_fisico'][1], 3, 'f = dp/(eta rho cp dT) con esa caída de presión')
rango('ReqKWkgFfis', _rq['kWkg_pide_ffis'][0], _rq['kWkg_pide_ffis'][1], 0,
      'requisito del rechazo con batería y el f por caída de presión [kW/kg]')
val('AireEnfrBase', round(_rq['aire_enfr_m3s'], -2), 0, 'aire de enfriamiento del caso base, al 100 m3/s [m3/s]')
val('AreaEntradaBase', round(_rq['area_entrada_m2'], -2), 0, 'área de entrada de ese aire a v_aire_entrada, al 100 m2 [m2]')
val('VAireEntrada', p.v_aire_entrada, 0, 'velocidad del aire en la entrada supuesta [m/s]')
val('SeccionCascoBase', round(_rq['seccion_casco_m2'], -2), 0, 'sección transversal del casco del caso base, al 100 m2 [m2]')
_hu = v3['humedad']['humedo']
rango('HRAire', 100 * p.HR_aire[0], 100 * p.HR_aire[1], 0, 'humedad relativa del día de diseño [%]')
rango('HumAgua', _hu[0]['agua_extra_kg_por_kgH2'], _hu[-1]['agua_extra_kg_por_kgH2'], 2,
      'agua que trae el aire de combustión húmedo y condensa [kg por kg de H2]')
rango('HumHdos', -100 * _hu[0]['dH2_t_rel'], -100 * _hu[-1]['dH2_t_rel'], 1, 'cuánto BAJA el H2 por tonelada con aire húmedo [%]')
rango('HumQt', -100 * _hu[0]['dQ_t_rel'], -100 * _hu[-1]['dQ_t_rel'], 1, 'cuánto BAJA el calor por tonelada con aire húmedo [%]')
_fs = rutas.leer_datos('helices_vectoriales')
_fr, _fb = _fs['referencia'], _fs['bajo']
_lim = lambda xs, k: (min(x[k] for x in xs), max(x[k] for x in xs))
val('FrenoTmin', _fr['t_quema_s'] / 60, 1, 'hélices existentes: quema de 1 t al rechazo de referencia [min]')
val('FrenoF', _fs['F_N'] / 1e3, 1, 'hélices existentes: fuerza inicial [kN]')
rango('FrenoA', p.freno_A_discos[0], p.freno_A_discos[1], 0, 'área total de discos supuesta [m2]')
val('FrenoFM', p.freno_FM, 1, 'figura de mérito supuesta')
rango('FrenoP', *[x / 1e3 for x in _lim(_fr['helices'], 'P_max_W')], 0, 'potencia de pico de las hélices [kW]')
rango('FrenoE', *[x / 3.6e6 for x in _lim(_fr['helices'], 'E_J')], 0, 'energía de las hélices por tonelada, 10 MW [kWh]')
rango('FrenoEPct', *[100 * x for x in _lim(_fr['helices'], 'E_frac_quim')], 1,
      'esa energía como % de la química del H2 que se quema por tonelada')
val('FrenoEQuim', _fr['E_quim_J'] / 3.6e6, 0, 'energía química (PCI) del H2 que se quema por tonelada [kWh]')
rango('FrenoHXAvion', _fr['rechazo']['avion_t'][0], _fr['rechazo']['avion_t'][1], 1,
      'unidad de rechazo de 10 MW con intercambiadores de avión (80-160 °C) [t]')
val('FrenoHXACC', _fr['rechazo']['acc_haces_t'], 0, 'unidad de rechazo de 10 MW con haces de ACC (temperatura del ciclo) [t]')
_fr_haces_t_por_MW = _fr['rechazo']['acc_haces_t'] / (p.P_ref_masa / 1e6)    # lineal en la potencia [t/MW]
val('FrenoPotBajo', p.freno_P_bajo / 1e6, 0, 'rechazo menor de referencia [MW]')
val('FrenoTminBajo', _fb['t_quema_s'] / 60, 0, 'quema de 1 t con el rechazo menor [min]')
val('FrenoHXACCBajo', _fb['rechazo']['acc_haces_t'], 1, 'haces de ACC del rechazo menor [t]')
rango('FrenoEBajo', *[x / 3.6e6 for x in _lim(_fb['helices'], 'E_J')], 0, 'energía de las hélices por tonelada, rechazo menor [kWh]')
val('FrenoEscGrande', _fs['escalon_grande']['m_kg'] / 1e3, 0, 'escalón grande de referencia [t]')
rango('FrenoPGrande', min(_fs['escalon_grande']['P_W']) / 1e6, max(_fs['escalon_grande']['P_W']) / 1e6, 1,
      'potencia de las hélices para ese escalón [MW]')
rango('FrenoEGrande', min(_fs['escalon_grande']['E_J']) / 3.6e9, max(_fs['escalon_grande']['E_J']) / 3.6e9, 1,
      'energía de las hélices para ese escalón, rechazo de referencia [MWh]')
_rp = _fs['respuesta']
val('FrenoTauFW', _rp['Flying Whales LCA60T']['tau_techo_s'], 1, 'cualquier freno: techo de respuesta, 1 t en el LCA60T [s]')
val('FrenoTauNT', _rp['Zeppelin NT']['tau_techo_s'], 1, 'cualquier freno: techo de respuesta, 1 t en el Zeppelin NT [s]')
# ── El freno como 4.ª caja negra (freno.py; sección nueva del paper 1, revisión 2026-10-06) ──
_fn = rutas.leer_datos('freno')
_ffn = lambda c, mt, P: next(x for x in _fn['filas'] if x['corto'] == c and x['m_kg'] == mt * 1000.0
                              and abs(x['P_W'] - P * 1e6) < 1.0)
_fwn = [x for x in _fn['filas'] if x['corto'] == 'FW']
_fwm = lambda mt: [x for x in _fwn if x['m_kg'] == mt * 1000.0]
_pmm = lambda xs, k, esc: (min(min(x[k]) for x in xs) / esc, max(max(x[k]) for x in xs) / esc)
rango('FrenoPGrilla', min(x['P_W'] for x in _fn['filas']) / 1e6, max(x['P_W'] for x in _fn['filas']) / 1e6, 0,
      'rechazo instalado de la tabla del freno; el mayor es el límite de la carga útil del LCA60T [MW]')
rango('FrenoMeFW', min(x['m_e_kg'] for x in _fwn), max(x['m_e_kg'] for x in _fwn), 0,
      'escalón que el ciclo compensa en tiempo real, LCA60T, grilla de la tabla [kg]')
val('FrenoMeDiezFW', _ffn('FW', 1, 10)['m_e_kg'], 0, 'escalón en tiempo real con 10 MW, LCA60T, 1 t [kg]')
val('FrenoMeDiezNT', _ffn('ZNT', 1, 10)['m_e_kg'], 0, 'escalón en tiempo real con 10 MW, Zeppelin NT, 1 t [kg]')
rango('FrenoFrac', 100 * min(x['frac_freno'] for x in _fn['filas']), 100 * max(x['frac_freno'] for x in _fn['filas']), 1,
      'parte del escalón que toma el freno, en la tabla [%]')
rango('FrenoTrecUno', min(x['t_rec_s'] for x in _fwm(1)) / 60, max(x['t_rec_s'] for x in _fwm(1)) / 60, 1,
      'tiempo de recuperación de 1 t, grilla de la tabla [min]')
val('FrenoTrecDiezT', _ffn('FW', 10, 10)['t_rec_s'] / 60, 0, 'tiempo de recuperación de 10 t con 10 MW [min]')
val('FrenoTrecTreintaT', _ffn('FW', 30, 10)['t_rec_s'] / 3600, 1, 'tiempo de recuperación de 30 t con 10 MW [h]')
rango('FrenoPHelUno', *_pmm(_fwm(1), 'P_hel_W', 1e6), 2, 'potencia de las hélices, 1 t, LCA60T, grilla de la tabla [MW]')
rango('FrenoEHelUno', *_pmm(_fwm(1), 'E_hel_J', 3.6e6), 0, 'energía de las hélices, 1 t, LCA60T, grilla de la tabla [kWh]')
rango('FrenoPHelDiezT', *_pmm(_fwm(10), 'P_hel_W', 1e6), 1, 'potencia de las hélices, 10 t, LCA60T [MW]')
rango('FrenoEHelDiezT', *_pmm(_fwm(10), 'E_hel_J', 3.6e9), 1, 'energía de las hélices, 10 t, LCA60T, grilla de la tabla [MWh]')
rango('FrenoPHelTreintaT', *_pmm(_fwm(30), 'P_hel_W', 1e6), 0, 'potencia de las hélices, 30 t, LCA60T [MW]')
rango('FrenoEHelTreintaT', *_pmm(_fwm(30), 'E_hel_J', 3.6e9), 0, 'energía de las hélices, 30 t, LCA60T, grilla de la tabla [MWh]')
rango('FrenoRegenDiezT', 10 * min(_fn['regen_kWh_por_t']) / 1e3, 10 * max(_fn['regen_kWh_por_t']) / 1e3, 0,
      'regenerar el H2 de 10 t (electrólisis y compresión, 1 atm) [MWh]')
# auditoría 2026-10-06 (M1, m12): el límite de la carga útil en el NT y el equipo del LCA60T con 10 MW
_lc = _fn['limite_carga']['Zeppelin NT']
val('FrenoPLimNT', _lc['P_lim_W'] / 1e6, 2, 'NT: rechazo con el que los haces de ACC pesan su carga útil [MW]')
val('FrenoTrecLimNT', _lc['t_rec_1t_s'] / 3600, 1, 'NT: recuperación de 1 t con ese rechazo [h]')
_bf10 = [b['m_bat_kg'] for b in rutas.leer_datos('exp11_balance_energia')['buffer'] if b['P_instalada_MW'] == 10.0]
_eq10 = _ffn('FW', 1, 10)['masa_acc_haces_t'] + _rq['fijas_t']
rango('FrenoEquipoDiez', _eq10 + min(_bf10) / 1e3, _eq10 + max(_bf10) / 1e3, 0,
      'LCA60T, 10 MW: haces + buffer de sus ventiladores + cajas fijas [t]')
val('FrenoPLimHaces', _fn['figura']['P_lim_haces_W'] / 1e6, 0, 'rechazo con el que los haces de ACC pesan la carga útil del LCA60T [MW]')
val('FrenoMargen', 100 * _fn['margen'], 0, 'margen del autor con el que se dimensiona la fuerza del freno [%]')
val('FrenoMu', _fn['mu'], 2, 'exceso de fuerza del freno, mu = 1/(1 - margen) - 1')
val('FrenoU', _fn['respuesta']['Flying Whales LCA60T']['u'], 3, 'raíz u de (1 + mu)(1 - e^-u) = mu u')
val('FrenoG', _fn['respuesta']['Flying Whales LCA60T']['G'], 3, 'G(mu) = u (1 - mu u / 2): subida mínima / (a tau^2)')
val('FrenoTauPrimerFW', _fn['respuesta']['Flying Whales LCA60T']['tau_primer_orden_s'], 1,
    'plazo del freno con fuerza (1+mu) m g y el mejor control, 1 t, LCA60T [s]')
val('FrenoTauPrimerNT', _fn['respuesta']['Zeppelin NT']['tau_primer_orden_s'], 1,
    'plazo del freno con fuerza (1+mu) m g y el mejor control, 1 t, Zeppelin NT [s]')
val('FrenoNTLiv', _fn['NT']['liviandad_cert_kg'], 0, 'liviandad estática certificada del Zeppelin NT (FAA TCDS AS1CE) [kg]')
val('FrenoNTCarga', _fn['NT']['carga_kg'] / 1e3, 0, 'carga útil del Zeppelin NT [t]')
val('FrenoTnatFW', next(x for x in _fwm(1))['t_nat_s'], 0, 'el ciclo solo recupera 1 t en el LCA60T en [s]')
val('FrenoTnatNT', _ffn('ZNT', 1, 10)['t_nat_s'], 1, 'el ciclo solo recupera 1 t en el Zeppelin NT en [s]')
_sfn = {x['nave']: x for x in _fn['simulacion']}
val('FrenoSimSubeNT', _sfn['Zeppelin NT']['sube_m'], 2, 'simulación con freno, NT, 1 t con 10 MW: subida [m]')
val('FrenoSimSubeFW', _sfn['Flying Whales LCA60T']['sube_m'], 2, 'simulación con freno, LCA60T, 1 t con 10 MW: subida [m]')
rango('FrenoSimERazon', min(x['razon_E'] for x in _fn['simulacion']), max(x['razon_E'] for x in _fn['simulacion']), 2,
      'energía del freno, simulación / forma cerrada')
rango('FrenoSimHdos', min(x['H2_kg'] for x in _fn['simulacion']), max(x['H2_kg'] for x in _fn['simulacion']), 1,
      'simulación con freno, 1 t con 10 MW: H2 quemado [kg]')
val('FrenoSimHdosEsp', _fn['simulacion'][0]['H2_esperado_kg'], 1, 'H2 de equilibrio, m / L_eff, para 1 t [kg]')
# para validar a mano en las derivaciones (1 t, LCA60T, 10 MW)
_f1 = _ffn('FW', 1, 10)
val('FrenoPrtFW', _f1['P_rt_W'] / 1e6, 0, 'rechazo en tiempo real de 1 t en el LCA60T, P_rt = dm* q [MW]')
val('FrenoAFW', _fn['respuesta']['Flying Whales LCA60T']['a_mps2'], 4, 'aceleración al soltar 1 t, LCA60T [m/s2]')
val('FrenoANT', _fn['respuesta']['Zeppelin NT']['a_mps2'], 3, 'aceleración al soltar 1 t, Zeppelin NT [m/s2]')
val('FrenoFDiezFW', _f1['F_N'] / 1e3, 2, 'fuerza del freno, 1 t con 10 MW, LCA60T [kN]')
rango('FrenoPHelDiezFW', min(_f1['P_hel_W']) / 1e3, max(_f1['P_hel_W']) / 1e3, 0,
      'potencia de las hélices, 1 t con 10 MW, LCA60T [kW]')
rango('FrenoEHelDiezFW', min(_f1['E_hel_J']) / 3.6e6, max(_f1['E_hel_J']) / 3.6e6, 1,
      'energía de las hélices, 1 t con 10 MW, LCA60T [kWh]')
val('FrenoTrecDiezFWs', _f1['t_rec_s'], 0, 'tiempo de recuperación de 1 t con 10 MW [s]')
val('FrenoKWkgHaces', _fn['figura']['kWkg_haces'], 3, 'haces de ACC: potencia específica usada por freno.py [kW/kg]')
# ronda 2 (2026-10-06): anticipar con el freno y prescindir de la quema previa
_an = _fn['anticipacion']
val('FrenoRampaFW', _an['Flying Whales LCA60T']['rampa_Nps'] / 1e3, 2, 'freno anticipado: rampa L g dm*, 1 t, LCA60T [kN/s]')
val('FrenoRampaNT', _an['Zeppelin NT']['rampa_Nps'] / 1e3, 1, 'freno anticipado: rampa L g dm*, 1 t, Zeppelin NT [kN/s]')
val('FrenoRampaPreviaFW', _an['Flying Whales LCA60T']['t_previa_s'], 1, 'freno anticipado: rampa antes de soltar, LCA60T [s]')
val('FrenoRampaPreviaNT', _an['Zeppelin NT']['t_previa_s'], 1, 'freno anticipado: rampa antes de soltar, Zeppelin NT [s]')
rango('FrenoSinPreviaPct', 100 * (min(x['factor_P_sin_previa'] for x in _fwm(1)) - 1),
      100 * (max(x['factor_P_sin_previa'] for x in _fwm(1)) - 1), 0,
      'sin quema previa: aumento de la potencia de las hélices, 1 t, LCA60T, grilla de la tabla [%]')
val('FrenoSinPreviaDiez', 100 * (_ffn('FW', 1, 10)['factor_P_sin_previa'] - 1), 0,
    'sin quema previa: aumento de la potencia de las hélices, 1 t, LCA60T, 10 MW [%]')
_tz = _fn['tau_nt']
val('FrenoTauEsperaNT', _tz['espera_s'], 1, 'NT, 1 t, 10 MW: retardo máximo del freno que espera, |dh| <= h_max [s]')
val('FrenoTauAnticipaNT', _tz['anticipa_s'], 1, 'NT, 1 t, 10 MW: retardo máximo del freno que anticipa, |dh| <= h_max [s]')
val('FrenoEsperaSubeNT', max(_tz['espera_con_tau_ant']), 1, 'NT: con ese retardo, el freno que espera deja subir [m]')
val('FrenoSinPreviaSubeNT', max(_tz['sin_previa_tau_base']), 2, 'NT, sin quema previa, retardo base: |dh| máximo [m]')
val('FrenoHelTau', p.hel_tau, 1, 'retardo de primer orden de la fuerza del freno en la simulación [s]')
val('FrenoHelWn', p.hel_wn, 1, 'frecuencia natural del lazo del freno en la simulación [rad/s]')
val('FrenoHelZeta', p.hel_zeta, 1, 'amortiguamiento del lazo del freno en la simulación')

r9 = rutas.leer_datos('exp09_condensacion_dinamica')
seccion('EXP-09: el minimax del artículo con dinámica del condensador (exp09_condensacion_dinamica.json)')
c9 = {x['caso']: x for x in r9['casos_base']}
val('NuevePhi',    r9['phi_max'], 2,        'fracción máxima que condensa (equilibrio, 40 °C, aire justo)')
val('NueveTau',    r9['tau_base_s'], 0,     'retardo térmico base del condensador [s]')
val('NuevePulmonT', r9['T_pulmon_C'], 0,    'temperatura del vapor retenido como gas [°C]')
# el pulmón sale del paper (M3): a 1 atm y bajo el rocío el vapor solo es gas mezclado con el inerte del escape
from exp05_condensacion_tiempo_real import AIRE_ESTEQ as _AE9
val('NuevePulmonInerte', (_AE9 - p.ratio_O2_H2) / (1.0 + p.ratio_O2_H2), 1,
    'N2 + Ar del escape por kg de vapor: lo que habría que guardar con él [kg/kg]')
for clave, suf in (('ideal', 'Ideal'), ('ventear_infMW', 'VentInf'), ('ventear_100MW', 'VentCien'),
                   ('ventear_10MW', 'VentDiez'), ('pulmon_infMW', 'PulInf'), ('pulmon_100MW', 'PulCien'),
                   ('pulmon_10MW', 'PulDiez'), ('tau_0s', 'TauCero'), ('tau_30s', 'TauTreinta'),
                   ('sin_prearranque_inf', 'SinPre')):
    x = c9[clave]
    val('Nueve' + suf + 'Ret',   100 * x['frac_retenida'], 0, f'agua retenida [%] — {clave}')
    val('Nueve' + suf + 'Sube',  x['dh_arriba_m'], 2,         f'subida máxima [m] — {clave}')
    val('Nueve' + suf + 'Hdos',  x['H2_quemado_kg'], 0,       f'H2 quemado [kg] — {clave}')
    # (autor, 2026-10-11: la figura 3 sale del paper; el agua retenida y venteada pasa al texto, en kg)
    val('Nueve' + suf + 'AguaRet',  x['agua_retenida_kg'], 0, f'agua retenida a bordo [kg] — {clave}')
    val('Nueve' + suf + 'AguaVent', x['agua_venteada_kg'], 0, f'agua venteada con el nitrógeno [kg] — {clave}')
val('NueveHmax', c9['ideal']['h_max_m'], 0, 'umbral de la maniobra del artículo [m]')
pl9 = [y for y in r9['plataformas'] if y['modo'] != 'ideal']
rango('NuevePlatSube', min(y['dh_arriba_m'] for y in pl9), max(y['dh_arriba_m'] for y in pl9), 1,
      'subida con 100 MW (ventear), rango de las 4 plataformas [m]')
rango('NuevePlatRet', 100 * min(y['frac_retenida'] for y in pl9), 100 * max(y['frac_retenida'] for y in pl9), 0,
      'agua retenida con 100 MW (ventear), rango de las 4 plataformas [%]')

# ── EXP-10: liberación por escalones con retardo del condensador (exp10_escalones_condensador.py) ──
r10 = rutas.leer_datos('exp10_escalones_condensador')
seccion('EXP-10: escalones de EXP-06 con retardo del condensador (exp10_escalones_condensador.json)')
c10 = r10['casos']
grandes = [x for x in c10 if x['corto'] != 'ZNT']


def sel(casos, caso):
    return [x for x in casos if x['caso'] == caso]


val('DiezKgtIdeal', r10['H2_por_t_ideal_kg'], 1, 'H2 por tonelada con condensación instantánea [kg/t]')
val('DiezKgtVenteo', r10['H2_por_t_venteo_kg'], 1, 'H2 por tonelada venteando H2 directamente [kg/t]')
for caso, suf in (('ideal', 'Ideal'), ('ventear_tau0s', 'TauCero'), ('ventear_tau5s', 'TauCinco'),
                  ('ventear_tau30s', 'TauTreinta'), ('pulmon_tau5s', 'Pulmon')):
    xs = sel(c10, caso)
    rango('Diez' + suf + 'Ret', 100 * min(x['frac_retenida'] for x in xs), 100 * max(x['frac_retenida'] for x in xs),
          0, f'agua retenida [%], 4 plataformas x 10/100 MW — {caso}')
    rango('Diez' + suf + 'Kgt', min(x['H2_por_t_kg'] for x in xs), max(x['H2_por_t_kg'] for x in xs),
          0, f'H2 quemado por tonelada [kg/t], 4 plataformas x 10/100 MW — {caso}')
    val('Diez' + suf + 'Sube', max(x['dh_arriba_m'] for x in xs), 2, f'subida máxima [m], peor caso — {caso}')
    val('Diez' + suf + 'Fallan', sum(not x['cumple_h_max'] for x in xs), 0, f'casos que no cumplen h_max (de 8) — {caso}')
    xg = sel(grandes, caso)
    rango('Diez' + suf + 'KgtGrandes', min(x['H2_por_t_kg'] for x in xg), max(x['H2_por_t_kg'] for x in xg),
          0, f'H2 por tonelada [kg/t], sin el Zeppelin NT — {caso}')
rango('DiezQuema', min(x['dur_quema_s'] for x in sel(c10, 'ideal')), max(x['dur_quema_s'] for x in sel(c10, 'ideal')),
      0, 'duración de la quema por escalón, ideal [s]')
rango('DiezEscalon', min(x['escalon_kg'] for x in sel(c10, 'ideal')), max(x['escalon_kg'] for x in sel(c10, 'ideal')),
      0, 'tamaño del escalón (EXP-06, h = 1 m, 10 y 100 MW) [kg]')

# ── EXP-11: cajas negras, balance de energía y regeneración (exp11_balance_energia.py) ──
r11 = rutas.leer_datos('exp11_balance_energia')
seccion('EXP-11: cajas negras, balance de energía y regeneración (exp11_balance_energia.json)')
aa = r11['autoabasto']
for caso, suf in (('1atm', 'Atm'), ('presurizado', 'Pres')):
    xs = [a for a in aa if a['caso'] == caso]
    rango('OnceEtaMin' + suf, min(a['eta_min'] for a in xs), max(a['eta_min'] for a in xs), 2,
          f'eta mínima para autoabastecer el enfriamiento — {caso}')
    rango('OnceRazon' + suf, min(a['razon'] for a in xs), max(a['razon'] for a in xs), 1,
          f'trabajo generado / pedido con eta de referencia — {caso}')
    rango('OncePide' + suf, min(a['E_pide_kWh_t'] for a in xs), max(a['E_pide_kWh_t'] for a in xs), 0,
          f'trabajo que pide el enfriamiento [kWh por tonelada] — {caso}')
    val('OnceGenera' + suf, xs[0]['E_gen_kWh_t'], 0, f'trabajo generado [kWh por tonelada] — {caso}')
    val('OncePhi' + suf, xs[0]['phi'], 2, f'fracción que condensa — {caso}')
val('OnceEtaRef', r11['eta_ref'], 2, 'eta de referencia (el grupo electrógeno del artículo)')
val('OncePresAtm', next(a['P_atm'] for a in aa if a['caso'] == 'presurizado'), 0, 'presión del caso presurizado [atm]')
bf = r11['buffer']
grandes_b = [b for b in bf if b['corto'] != 'ZNT']
rango('OnceKgKW', min(b['kg_por_kW'] for b in bf), max(b['kg_por_kW'] for b in bf), 2,
      'radiador + buffer por kW de condensación instalada [kg/kW]')
for P, suf in ((10.0, 'Diez'), (100.0, 'Cien')):
    xs = [b for b in bf if b['P_instalada_MW'] == P]
    rango('OnceBuf' + suf, min(b['m_bat_kg'] for b in xs) / 1e3, max(b['m_bat_kg'] for b in xs) / 1e3,
          1 if P < 50 else 0, f'masa del buffer [t] — {P:g} MW')
    rango('OnceTasaC' + suf, min(b['tasa_C'] for b in xs), max(b['tasa_C'] for b in xs), 0,
          f'tasa de descarga del buffer [C] — {P:g} MW')
    xg = [b for b in grandes_b if b['P_instalada_MW'] == P]
    rango('OnceTasaCGrandes' + suf, min(b['tasa_C'] for b in xg), max(b['tasa_C'] for b in xg), 0,
          f'tasa de descarga de un buffer dimensionado por energía [C], sin el Zeppelin NT — {P:g} MW')
    rango('OnceTasaCPaqGrandes' + suf, min(b['tasa_C_paquete'] for b in xg), max(b['tasa_C_paquete'] for b in xg), 0,
          f'tasa de descarga del buffer dimensionado por potencia [C], sin el Zeppelin NT — {P:g} MW')
    rango('OnceCarga' + suf, 100 * min(b['frac_carga'] for b in xg), 100 * max(b['frac_carga'] for b in xg), 0,
          f'radiador + buffer [% de la carga útil], sin el Zeppelin NT — {P:g} MW')
rango('OnceRecarga', min(b['recarga'] for b in bf), max(b['recarga'] for b in bf), 1,
      'energía que genera la quema del escalón / energía del buffer')
val('OnceBufManda', sum(b['manda'] == 'potencia' for b in bf), 0, 'casos en que la masa del buffer la fija la potencia')
val('OnceBufCasos', len(bf), 0, 'casos de buffer evaluados')
rg = r11['regeneracion']
gr = [g for g in rg if g['corto'] != 'ZNT']
val('OnceVenteo', r11['H2_por_t_venteo_kg'], 1, 'H2 perdido por tonelada venteando H2 [kg/t]')
for caso, suf in (('ventear_tau0s', 'TauCero'), ('ventear_tau5s', 'TauCinco')):
    xs = [g for g in gr if g['caso'] == caso]
    rango('OncePerdido' + suf, min(g['perdido_kg_t'] for g in xs), max(g['perdido_kg_t'] for g in xs), 0,
          f'H2 perdido por tonelada, plataformas grandes [kg/t] — {caso}')
    rango('OnceAhorro' + suf, 100 * (1 - max(g['perdido_kg_t'] for g in xs) / r11['H2_por_t_venteo_kg']),
          100 * (1 - min(g['perdido_kg_t'] for g in xs) / r11['H2_por_t_venteo_kg']), 0,
          f'pérdida evitada frente a ventear H2 [%] — {caso}')
    rango('OnceDias' + suf, min(g['dias_por_t'] for g in xs), max(g['dias_por_t'] for g in xs), 2,
          f'días de regeneración solar por tonelada, plataformas grandes — {caso}')
zn = [g for g in rg if g['corto'] == 'ZNT' and g['caso'] == 'ventear_tau5s']
rango('OncePerdidoZNT', min(g['perdido_kg_t'] for g in zn), max(g['perdido_kg_t'] for g in zn), 0,
      'H2 perdido por tonelada, Zeppelin NT, retardo 5 s [kg/t]')
sl = r11['sumidero_LN2']
val('OnceQPorT',        sl['Q_por_t_GJ'], 1, 'calor a rechazar por tonelada liberada, condensando todo [GJ/t]')
val('OnceLNdos',        sl['LN2_t_por_t'], 0, 'LN2 que habría que evaporar por tonelada liberada [t/t]')
val('OnceLNdosLatente', sl['LN2_latente_t_por_t'], 1, 'LN2 por tonelada contando solo el latente del agua [t/t]')

# ── EXP-12: comparación con helio (exp12_helio.py) ──
r12 = rutas.leer_datos('exp12_helio')
seccion('EXP-12: comparación con helio a nivel de caja negra (exp12_helio.json)')
ru = {x['ruta'].split(' (')[0].split(',')[0] + (' 1atm' if '1 atm' in x['ruta'] else ''): x for x in r12['rutas']}
rh, rhv, rhe, rhev = ru['H2 quemado'], ru['H2 venteado'], ru['He comprimido a bordo'], ru['He venteado']
mismo(rh['L'], rutas.leer_datos('exp06_requisitos_minimos')['L_phi'], 'L(phi) H2 (EXP-12 vs EXP-06)')
val('DoceRatioRho',  r12['rho_He_sobre_rho_H2'], 2, 'rho_He / rho_H2 a igual T y P')
val('DoceLHdos',     rh['L'], 1,       'pesadez por kg de H2 quemado con agua retenida (phi objetivo) [kg/kg]')
val('DoceLHdosVent', rhv['L'], 1,      'pesadez por kg de H2 venteado [kg/kg]')
val('DoceLHe',       rhe['L'], 2,      'pesadez por kg de He comprimido a bordo [kg/kg]')
val('DoceLHeVent',   rhev['L'], 2,     'pesadez por kg de He venteado [kg/kg]')
val('DoceKgHdos',    rh['kg_por_t'], 0,  'H2 quemado por tonelada liberada [kg/t]')
val('DoceKgHe',      rhe['kg_por_t'], 0, 'He comprimido por tonelada liberada [kg/t]')
val('DoceKgHeVent',  rhev['kg_por_t'], 0, 'He perdido por tonelada liberada si se ventea [kg/t]')
val('DoceSuccionHdos', rh['m3_por_t'], 0, 'gas procesado por tonelada, H2 [m3/t a 1 atm]')
val('DoceSuccionHe', rhe['m3_por_t'], 0, 'gas procesado por tonelada, He comprimido [m3/t a 1 atm]')
val('DoceCalorHdos', rh['calor_GJ_por_t'], 1, 'calor a rechazar por tonelada, H2 con agua retenida [GJ/t]')
alm = {x['P_MPa']: x for x in r12['almacen']}
for P, suf in ((5.0, 'Cinco'), (20.0, 'Veinte'), (30.0, 'Treinta')):
    val('DoceW' + suf, alm[P]['W_iso_MJkg'], 2, f'trabajo mínimo isotérmico para comprimir He a {P:g} MPa [MJ/kg]')
    val('DoceCalorHe' + suf, alm[P]['calor_GJ_por_t'], 2, f'calor (= trabajo) por tonelada, He a {P:g} MPa [GJ/t]')
    val('DoceAlmacen' + suf, alm[P]['m3_almacen_por_t'], 1, f'volumen de He comprimido a {P:g} MPa por tonelada [m3/t]')
val('DoceFactorCalor', r12['factor_calor_20MPa'], 0, 'calor por tonelada H2 / He a 20 MPa')
bh, be = r12['base']['H2'], r12['base']['He']
val('DoceBaseDmHdos', bh['dm_kgps'], 2, 'caso base caso base 1 t: flujo minimax H2 (forma cerrada exacta) [kg/s]')
val('DoceBaseDmHe',   be['dm_kgps'], 1, 'caso base caso base 1 t: flujo minimax He comprimido [kg/s]')
val('DoceBaseAlphaHe', be['alpha'], 3, 'caso base caso base 1 t: alpha* con He')
val('DoceFactorFlujo', r12['factor_flujo_base'], 2, 'flujo He / H2 con la misma forma cerrada')
val('DoceFactorL',   r12['factor_flujo_L'], 2, 'flujo He / H2 por regla de tres con L (no es el minimax)')
val('DoceBaseSuccionHdos', bh['succion_m3ps'], 0, 'caso base: succión H2 [m3/s a 1 atm]')
val('DoceBaseSuccionHe',   be['succion_m3ps'], 0, 'caso base: succión He [m3/s a 1 atm]')
val('DoceBasePotHdos', bh['potencia_MW'], 0, 'caso base: calor a rechazar en tiempo real, H2 [MW]')
val('DoceBasePotHe',   be['potencia_MW']['20.0'], 0, 'caso base: trabajo mínimo de compresión He a 20 MPa [MW]')
# revisión 2026-10-06: el tiempo de la maniobra limpia no depende del gas; potencia en el tiempo de recuperación
import minimax as _mmD
from exp05_condensacion_tiempo_real import q_rech as _qD
from exp02_quema_centrada import L_ef as _LD
_tnatH = 1000.0 / _mmD.L_eff / bh['dm_kgps']
_tnatHe = 1000.0 / r12['rutas'][3]['L'] / be['dm_kgps']
val('DoceTnatHdos', _tnatH, 1, 'caso base: duración de la maniobra limpia, H2 [s]')
val('DoceTnatHe', _tnatHe, 1, 'caso base: duración de la maniobra limpia, He (la misma: no depende de L) [s]')
# (2026-10-14, objeción 1.9) la quema se dimensiona con phi = 1 y la unidad retiene phi_objetivo: con la tasa
# proporcional a 1 / L (ec. de la quema limpia), cuánto más alta saldría
val('PhiUnoDifPct', 100.0 * (_LD(1.0) / _LD(p.phi_objetivo) - 1.0), 1,
    'tasa instalada con phi_objetivo en lugar de phi = 1: cuánto más alta [%]')
_trecD = (1000.0 / _LD(p.phi_objetivo)) * _qD(0.0, p.phi_objetivo) / p.P_ref_masa
val('DocePotHeRec', r12['rutas'][3]['calor_GJ_por_t'] * 1e9 / _trecD / 1e6, 1,
    'COSH: potencia de compresión para 1 t en el tiempo de recuperación del ciclo con 10 MW [MW]')
# comparación en dos etapas (2026-10-06): almacenamiento de lo que se quita y recuperación
_alm20 = next(x for x in r12['almacen'] if abs(x['P_MPa'] - r12['P_ref_MPa']) < 1e-9)
val('DoceTanqueHe', _alm20['tanque_compuesto_t_por_t'], 2, 'COSH: compuesto del tanque por tonelada liberada, cota inferior [t]')
val('DoceTanqueHeCarga', _alm20['tanque_compuesto_t_por_t'] * p.NAVE_BASE['carga_max'] / 1e3, 0,
    'COSH: compuesto de los tanques para toda la carga útil del LCA60T, cota inferior [t]')
val('DoceAguaM', r12['rutas'][0]['agua_m3_por_t'], 2, 'agua retenida por tonelada liberada (phi = 0.99) [m3]')
_tq = next(c for c in rutas.leer_datos('balance_masa')['cajas'] if c['caja'] == 'tanques_por_t')
rango('DoceTanqueCiclo', _tq['min_kg'] / 1e3, _tq['max_kg'] / 1e3, 2, 'tanques de H2 del ciclo, una tonelada regenerada [t]')
from exp16_balance_neto import e_comp as _ecD
_Rg = r12['rutas'][0]['kg_por_t'] * (1.0 - r12['rutas'][0]['perdida_frac'])
rango('DoceRecupHdos', _Rg * (p.E_H2_electrolisis + min(_ecD())) / 1e3, _Rg * (p.E_H2_electrolisis + max(_ecD())) / 1e3, 1,
      'energía para devolver la sustentación de 1 t: electrólisis y compresión del H2 regenerable (phi = 0.99) [MWh]')
# paso [5] (corrección del autor): las dos rutas devuelven la sustentación sacando gas de sus tanques a bordo; el ciclo
# saca los B kg que quemó y la electrólisis devuelve el resto (el O2 del agua)
_h0 = r12['rutas'][0]
val('DoceInstantePct', _h0['al_instante_kg'] / 10.0, 0, 'ciclo: sustentación devuelta al instante con los B kg de sus tanques [% de 1 t]')
val('DoceRestoKg', _h0['resto_kg'], 0, 'ciclo: resto que vuelve con la electrólisis, B (r phi - 1) [kg por t]')
val('DoceTanqueCicloComp', _h0['reserva_compuesto_t'], 2, 'ciclo: compuesto de los tanques de B kg a 35 MPa, cota inferior [t por t]')
val('DocePTanqueMPa', _h0['reserva_P_MPa'], 0, 'presión de los tanques del ciclo en la comparación [MPa]')
val('DoceTodoKgHdos', r12['reserva_todo']['H2_kg'], 1, 'H2 de los tanques para devolver 1 t entera de golpe [kg]')
val('DoceTanqueTodo', r12['reserva_todo']['H2_compuesto_t'], 2, 'compuesto para devolver 1 t entera de golpe con H2 a 35 MPa [t]')
for e in r12['escalones']:
    suf = {10.0: 'Diez', 100.0: 'Cien'}[e['P_MW']]
    val('DoceMinHdos' + suf, e['H2_min_por_t'], 1, f"minutos por tonelada con {e['P_MW']:g} MW, H2")
    val('DoceMinHe' + suf, e['He_min_por_t']['20.0'], 2, f"minutos por tonelada con {e['P_MW']:g} MW, He a 20 MPa")

# ── EXP-13: masa añadida en la dinámica (inercia.py; datos de exp06_requisitos_minimos.json) ──
seccion('EXP-13: masa añadida en la dinámica (inercia.py)')
t13 = r6['base_1t_h1']
M13 = t13['M_total_kg']
_aire13 = p.rho_aire_std * p.NAVE_BASE['V_casco']            # aire del volumen del casco
_I13 = _mmx.resolver(p.nombre_dirigible, p.m_carga, p.h_max_precision, p.D_lim_oper)['M_inercial']
mismo(_I13, (1.0 + p.NAVE_BASE['k']) * _aire13, 'inercia del caso base = (1 + k) rho V_casco')
# M-12 (auditoría 3): esbelteces de las naves (el texto daba "3-6", el ZR-3 tiene 7.3) y cuánto se aparta la bajada
# pico simulada de la forma cerrada (1 + s) D_lim (el apéndice A decía "a few millimetres" sin número)
rango('EsbeltezNaves', min(d_['esbeltez'] for d_ in p.PLATAFORMAS.values()),
      max(d_['esbeltez'] for d_ in p.PLATAFORMAS.values()), 1, 'esbelteces L/D de las naves')
import minimax as _mmc
_mpc = rutas.leer_datos('minimax_plataformas')['single_1t']
val('CotaSimFormulaMm', 1e3 * max(abs(_mpc[n_]['dh_peak_m'] - (1.0 + _mmc.resolver(n_, p.m_carga, p.h_max_precision,
                                                                                     p.D_lim_oper)['s']) * p.D_lim_oper)
                                  for n_ in _mpc), 0,
    'bajada pico simulada frente a (1 + s) D_lim, la mayor de las cuatro naves [mm]')
rango('KNaves', min(d_['k'] for d_ in p.PLATAFORMAS.values()), max(d_['k'] for d_ in p.PLATAFORMAS.values()), 2,
      'k de masa añadida de las naves, cada una con su esbeltez (Lamb)')
val('TreceUnoMasK',     1.0 + p.NAVE_BASE['k'], 2, 'inercia / aire del volumen del casco en el equilibrio, 1 + k (caso base)')
val('TreceMasaBase',      M13 / 1000.0, 0, 'masa total caso base con 1 t [t]')
val('TreceMasaAireBase',  _aire13 / 1000.0, 0, 'aire del volumen del casco del caso base, rho V_casco [t]')
val('TreceVcascoBase',    p.NAVE_BASE['V_casco'], 0, 'volumen total del casco del caso base [m3]')
# M-2 (auditoría 3): el LCA60T usa su volumen de helio publicado (180 000 m3) como casco, cota inferior; con los
# 200 000 m3 de la FAQ (fuentes.FWH26) el flujo instalado baja
import minimax as _mmv
_dmv = _mmv.resolver(p.NAVE_BASE, p.m_carga, p.h_max_precision, p.D_lim_oper)['dm_star']
_dmv2 = _mmv.resolver(dict(p.NAVE_BASE, V_casco=p.V_casco_FW_alternativo), p.m_carga, p.h_max_precision, p.D_lim_oper)['dm_star']
val('BaseVolAltPct', 100.0 * (1.0 - _dmv2 / _dmv), 0, 'flujo instalado con 200 000 m3 en vez de 180 000: cuánto baja [%]')
val('TreceInerciaBase',   _I13 / 1000.0, 0, 'inercia vertical caso base con 1 t [t]')
val('TreceDmConBase',     t13['dm_star_kgps'], 2, 'ṁ* caso base 1 t, h = 1 m, con masa añadida [kg/s]')
val('TreceDmSinBase',     t13['dm_star_sin_ma_kgps'], 2, 'ṁ* caso base 1 t, h = 1 m, sin masa añadida [kg/s]')
val('TrecePcondConBase',  t13['P_cond_MW'], 0, 'P_cond en tiempo real caso base 1 t, con masa añadida [MW]')
val('TrecePcondSinBase',  t13['P_cond_sin_ma_MW'], 0, 'P_cond en tiempo real caso base 1 t, sin masa añadida [MW]')
val('TreceZetaAnterior', 1.35 * math.sqrt(   # 1.35: zeta de aquel diseño (2026-10-09: el actual es 1.30)
    (p.NAVE_BASE['m_vacio'] + p.m_carga) / _I13), 2,
    'zeta real que tenía el diseño anterior (Kd sin masa añadida ni H2), caso base 1 t')

# ── VAL-01: validación por componente contra datos publicados (val01_componentes.py) ──
v1 = rutas.leer_datos('val01_componentes')
seccion('VAL-01: validación por componente contra datos publicados (val01_componentes.json)')
vr, va, vn = v1['referencias_derivadas'], v1['cifras_antes'], v1['cifras_ahora']
vf = {f['magnitud']: f for f in v1['filas']}
vs = {s['supuesto'].split(' (')[0]: s for s in v1['supuestos']}
mismo(vn['H2_por_t_kg'], 1000.0 / (p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1)) + p.ratio_O2_H2),
      'H2 por tonelada (VAL-01 vs parametros)')
val('ValidFilas',   v1['n_filas'], 0,  'comparaciones modelo vs dato publicado')
val('ValidDentro',  v1['n_dentro'], 0, 'comparaciones dentro de su criterio')
val('ValidOdos',    vr['r_O2_H2'], 3,  'O2 por kg de H2, estequiometría exacta (IUPAC)')
mismo(p.ratio_O2_H2, vr['r_O2_H2'], 'O2 por kg de H2 (parametros vs VAL-01)', tol=1e-3)
val('OdosLastre',   p.ratio_O2_H2, 1,  'lastre que dejaría guardar el O2 a bordo, por kg de H2 regenerado [kg/kg]')
val('ValidAgua',    vr['r_H2O_H2'], 3, 'H2O por kg de H2, estequiometría exacta (IUPAC)')
# (2026-10-14, autor) introducción: el agua guarda el H2 quemado en ~1/1300 de su volumen de gas. Gas a T1 y 1 atm (la
# densidad del modelo); agua líquida saturada (IAPWS) a la salida de la unidad de rechazo, donde se guarda (como en
# balance_masa). Redondeada a centenas: el texto dice "about".
from exp05_condensacion_tiempo_real import rho_liq_sat as _rlsV
_razV = _rlsV(p.T_salida_cond) / (vr['r_H2O_H2'] * p.P_atmos / (p.R_H2 * p.T1))
mismo(_razV, (p.R_univ * p.T1 / (p.M_H2 * p.P_atmos)) / ((p.M_H2O / p.M_H2) / _rlsV(p.T_salida_cond)),
      'volumen del H2 / volumen de su agua (VAL-01 vs masas molares de parametros)', tol=1e-3)
val('RazonVolHdosAgua', round(_razV, -2), 0,
    'volumen del H2 quemado (gas a T1, 1 atm) / volumen de su agua líquida a la salida de la unidad de rechazo [veces]')
val('ValidPCI',     vr['LHV_Jkg'] / 1e6, 2, 'PCI del H2, CODATA [MJ/kg]')
val('ValidPCS',     vr['HHV_Jkg'] / 1e6, 2, 'PCS del H2, CODATA [MJ/kg]')
val('ValidPCSModelo', vf['PCS = PCI + (M_H2O/M_H2) L_v']['modelo'], 2, 'PCS reconstruido con PCI + (M_H2O/M_H2) L_v [MJ/kg]')
val('ValidHfg',     vr['h_fg_25C_Jkg'] / 1e3, 0, 'calor latente del agua a 25 °C, NIST [kJ/kg]')
val('ValidPsatDesv', 100 * max(abs(f['desvio_rel']) for f in v1['filas'] if f['magnitud'].startswith('p_sat')), 3,
    '|desvío| máx. de p_sat (IF97) frente a IAPWS-95 de NIST [%]')
val('ValidRhoHDesv', 100 * vf['rho_H2 a 15 °C, 1 atm (gas ideal)']['desvio_rel'], 2,
    'desvío de rho_H2 gas ideal frente a NIST [%]')
val('ValidHsISA',   vr['H_escala_isa_h0_m'], 0, 'altura de escala local de la densidad ISA en h0 [m]')
val('ValidNinerte',  vf['N inerte por mol de H2O (escape)']['modelo'], 2, 'mol de N2 + Ar por mol de H2O en el escape estequiométrico')
val('ValidHs',      p.H_escala, 0, 'altura de escala de la densidad del modelo [m]')
val('ValidHsKm',    p.H_escala / 1000.0, 1, 'la misma altura de escala, en km (texto: tres cifras significativas)')
val('ValidEfIEA',   100 * vr['efic_iea_PCI'], 0, 'eficiencia de electrólisis IEA 2024 [% PCI]')
val('ValidElecIEA', vr['E_el_iea_kWhkg'], 1, 'energía de electrólisis, promedio IEA 2024 [kWh/kg]')
val('ValidElecTermo', vr['E_el_termoneutral_kWhkg'], 1, 'energía termoneutral de electrólisis (PCS) [kWh/kg]')
vma = v1['masa_anadida']
val('ValidMaRestringido', vma['medida_sobre_analitica']['restringido'], 2,
    'masa vertical medida / analítica, ensayo restringido (Yamasaki & Goto 2003)')
val('ValidMaVueloLibre', vma['medida_sobre_analitica']['vuelo_libre'], 2,
    'masa vertical medida / analítica, vuelo libre (Yamasaki & Goto 2003)')
val('ValidMaFlujoAlto', 100 * (vma['efecto_flujo']['restringido'] - 1), 0,
    'flujo dinámico si la inercia fuera la del ensayo restringido [+%]')
val('ValidMaFlujoBajo', 100 * (1 - vma['efecto_flujo']['vuelo_libre']), 0,
    'flujo dinámico si la inercia fuera la del vuelo libre [-%]')
# supuestos que quedan y su efecto
val('ValidGasReal', 100 * abs(vs['H2 como gas ideal']['cambio_rel']), 2, '|efecto| del gas ideal en H2 por tonelada [%]')
val('ValidRhoHcero', 100 * abs(vs['aire a 1.225 kg/m3 en h0']['cambio_rel']), 1,
    '|diferencia| de densidad entre ISA nivel del mar e ISA en h0 [%]')
val('ValidQQuince', 100 * vs['aire exterior a %g °C' % p.T_aire_exterior]['cambio_rel'], 1,
    'cambio de q con aire exterior a 15 °C [%]')
val('ValidQLatente', 100 * vs['q con calor sensible']['cambio_rel'], 1, 'cambio de q contando solo el latente [%]')
val('ValidElecOpt', 100 * abs(vs['electrólisis al 63 % del PCI']['cambio_rel']), 0,
    'reducción del tiempo de regeneración con 45 kWh/kg [%]')
val('ValidHsEfecto', abs(vs['H_s local ISA']['cambio_abs_mm']), 1,
    '|cambio| de la excursión máx. de la simulación del artículo con H_s = 8500 m [mm]')
# cifras con los valores anteriores y con los actuales
for k, nom, dec in (('L_eff', 'Leff', 2), ('H2_por_t_kg', 'Hdos', 1), ('q_MJkg', 'Q', 1),
                    ('phi_max_40C_1atm', 'Phi', 3), ('T_rocio_1atm_C', 'Rocio', 1),
                    ('P_phi_objetivo_atm', 'Pres', 1), ('E_electrolisis_kWhkg', 'Elec', 1)):
    val('ValidAntes' + nom, va[k], dec, f'{k} con los valores anteriores (redondeados)')
    val('ValidAhora' + nom, vn[k], dec, f'{k} con los datos publicados')

# ── EXP-14: calor por nivel de temperatura (exp14_calor_por_temperatura.py) ──
r14 = rutas.leer_datos('exp14_calor_por_temperatura')
seccion('EXP-14: calor por nivel de temperatura (exp14_calor_por_temperatura.json)')
f14 = {(x['P_atm'], round(x['eta'], 2), x['tec']): x for x in r14['filas']}
b0 = f14[(1.0, 0.0, '2026')]
val('CatorceRocio', b0['T_rocio_C'], 0, 'punto de rocío del escape a 1 atm [°C]')
val('CatorceFria', b0['q_fria_MJkg'], 1, 'calor de la banda fría (rocío -> 40 °C), 1 atm [MJ/kg H2]')
val('CatorceFracFriaCero', 100 * b0['frac_fria'], 0, 'fracción fría del calor, eta = 0 [%]')
val('CatorceFracFriaEta', 100 * f14[(1.0, round(p.rend_gen, 2), '2026')]['frac_fria'], 0, 'fracción fría, eta de referencia [%]')
# (2026-10-13) la banda fría de los dos casos de retención (auditoría, objeción 5.1): 1 atm y la presión que retiene 0.99
_bd = r14['bandas']
mismo(_bd['1atm']['frac'], b0['frac_fria'], 'banda fría a 1 atm: la fracción de bandas = la de la fila de 1 atm (W = 0)')
val('CatorceFriaObj', _bd['objetivo']['Q_banda_MJkg'], 1, 'calor de la banda fría a la presión que retiene 0.99, '
    'rocío -> 40 °C [MJ/kg H2]')
val('CatorceFracFriaObj', 100 * _bd['objetivo']['frac'], 0, 'fracción de q(0, 0.99) que sale en esa banda [%]')
rango('CatorceKWkgFria', b0['kWkg_fria_efectivo'], f14[(1.0, 0.0, '2035')]['kWkg_fria_efectivo'], 1,
      'kW/kg efectivo del radiador en la banda fría, 1 atm, 2026--2035')
rango('CatorceKgKW', b0['kg_por_kW'], f14[(1.0, 0.0, '2035')]['kg_por_kW'], 3,
      'masa de rechazo con dos bandas, 1 atm, eta = 0, 2026--2035 [kg/kW]')
rango('CatorceRazonArt', b0['razon_vs_articulo'], f14[(1.0, 0.0, '2035')]['razon_vs_articulo'], 2,
      'dos bandas / cifra única del artículo, 1 atm, 2026--2035')
for P, suf in ((5.0, 'Cinco'), (20.0, 'Veinte')):
    x = f14[(P, 0.0, '2026')]
    val('CatorceRocio' + suf, x['T_rocio_C'], 0, f'rocío a {P:g} atm [°C]')
    val('CatorceFin' + suf, x['T_fin_C'], 0, f'fin de la condensación a {P:g} atm [°C]')
    val('CatorceKWkgFria' + suf, x['kWkg_fria_efectivo'], 1, f'kW/kg banda fría a {P:g} atm, 2026')
    val('CatorceRazon' + suf, x['kg_por_kW'] / b0['kg_por_kW'], 2, f'masa de rechazo a {P:g} atm / a 1 atm, 2026')
    val('CatorceW' + suf, x['W_comp_MJkg'], 1, f'compresión isotérmica mínima a {P:g} atm [MJ/kg H2]')
t14 = r14['tina']
val('CatorceTinaDiez', t14['agua_t_por_t']['10.0'] if '10.0' in t14['agua_t_por_t'] else t14['agua_t_por_t'][10.0], 0,
    'agua de la tina para guardar el calor frío de 1 t subiendo 10 °C [t]')
val('CatorceLastreFrac', 100 * t14['frac_lastre_10K'], 1, 'calor frío que absorbe el propio lastre subiendo 10 °C [%]')

# ── Irradiancia sobre los paneles del casco (irradiancia_casco.py) ───────────
_ir = rutas.leer_datos('irradiancia_casco')
seccion('Irradiancia sobre el casco: SENAMHI (Lima) x factor de forma (irradiancia_casco.json)')
_EST = (('invierno', 'Inv'), ('verano', 'Ver'), ('anual', 'Anual'))
val('IrrTablaAnual', _ir['anual_tabla'], 2, 'SENAMHI, energía solar diaria anual, tabla 2 [kWh/m2/día]')
for est, se in _EST:
    e = _ir['estaciones'][est]
    val('IrrGhi' + se, e['GHI'], 1, f'irradiación global horizontal medida, {est} [kWh/m2/día]')
    val('IrrForma' + se, e['F'], 2, f'factor de forma de la mitad superior del casco, {est}')
    val('IrrCasco' + se, e['efectiva'], 1, f'irradiación sobre los paneles del casco, {est} [kWh/m2/día]')
    val('IrrFd' + se, 100 * e['fd'], 0, f'fracción difusa (Page 1964), {est} [%]')
rango('IrrFormaRumbo', min(e['F_min'] for e in _ir['estaciones'].values()),
      max(e['F_max'] for e in _ir['estaciones'].values()), 2, 'factor de forma según el rumbo del eje (N-S a E-W)')
rango('IrrFormaFd', min(e['F_fd0'] for e in _ir['estaciones'].values()),
      max(e['F_fd1'] for e in _ir['estaciones'].values()), 2, 'factor de forma con toda la luz directa o toda difusa')
val('IrrDifusa', _ir['F_difusa'], 2, 'difusa isotrópica sobre la mitad superior: 1/2 + 1/pi')
val('IrrAlbedo', _ir['albedo'], 1, 'reflectancia del suelo (supuesto estándar)')
val('IrrDifAnual', 100 * (_ir['anual_reconstruido'] / _ir['anual_tabla'] - 1), 1,
    'anual rearmado con los anexos 6 y 7 vs tabla 2 [%]')

# ── Balance de masa por caja negra (balance_masa.py) ─────────────────────────
bm = rutas.leer_datos('balance_masa')
seccion('Balance de masa por caja negra (balance_masa.json)')
val('MasaPref', bm['P_ref_W'] / 1e6, 0, 'potencia de rechazo de referencia [MW]')
val('MasaQuema', bm['mdot_kgps'], 3, 'quema que produce esa potencia de calor [kg/s]')
val('MasaEje', bm['P_eje_W'] / 1e6, 1, 'potencia de eje de la unidad de conversión [MW]')
rango('MasaKWkgTurb', min(bm['kWkg_turbina']), max(bm['kWkg_turbina']), 1, 'turbina CT7 (EASA TCDS) [kW/kg]')
_kr = [w / 1e3 / m for w, m in p.MOTOR_ROTAX]
rango('MasaKWkgRotax', min(_kr), max(_kr), 1, 'Rotax 916 iS, continua instalado / despegue [kW/kg] (proxy de aviación)')
val('MasaKWkgDeutz', p.MOTOR_DEUTZ[0][0] / 1e3 / p.MOTOR_DEUTZ[0][1], 2, 'DEUTZ TCG 7.8 H2 [kW/kg]')
bi = bm['ideal']
seccion('Tabla 6: cajas negras con el flujo IDEAL del caso base (caso base, 1 t de una vez)')
mismo(bi['mdot_kgps'], rutas.leer_datos('minimax_plataformas')['single_1t'][p.nombre_dirigible]['dm_opt'],
      'ṁ* de la tabla 6 = ṁ* del caso base simulado')
val('MasaIdealQuema', bi['mdot_kgps'], 2, 'quema ideal del caso base [kg/s]')
val('MasaIdealEje', bi['P_eje_W'] / 1e6, 0, 'potencia de eje de la conversión con el flujo ideal [MW]')
val('MasaIdealPrej', bi['P_rej_W'] / 1e6, 0, 'calor a rechazar con eta de referencia y phi objetivo [MW]')
val('MasaIdealPrejCero', bi['P_rej_eta0_W'] / 1e6, 0, 'calor a rechazar sin exportar trabajo (= tabla 5, caso base 1 t) [MW]')
val('MasaIdealTbuf', bi['t_buf_s'], 0, 'tiempo del buffer: 3 tau + quema [s]')
_cj = {c['caja']: c for c in bi['cajas']}
rango('MasaIdealConv', _cj['conversion']['min_kg'] / 1e3, _cj['conversion']['max_kg'] / 1e3, 0, 'unidad de conversión [t]')
rango('MasaIdealRech', _cj['rechazo']['min_kg'] / 1e3, _cj['rechazo']['max_kg'] / 1e3, 0, 'unidad de rechazo [t]')
rango('MasaIdealBuf', _cj['buffer']['min_kg'] / 1e3, _cj['buffer']['max_kg'] / 1e3, 0, 'buffer [t]')
rango('MasaIdealTotal', bi['total_kg'][0] / 1e3, bi['total_kg'][1] / 1e3, 0, 'total de las cajas negras, flujo ideal [t]')
rango('MasaIdealVeces', bi['veces_carga'][0], bi['veces_carga'][1], 1, 'total / carga útil del caso base [veces]')
# extremos sueltos para el abstract ("between 223 and 813 tonnes", "by 22 to 81 times")
val('MasaIdealTotalMin', bi['total_kg'][0] / 1e3, 0, 'total de las cajas negras, flujo ideal, optimista [t]')
val('MasaIdealTotalMax', bi['total_kg'][1] / 1e3, 0, 'total de las cajas negras, flujo ideal, pesimista [t]')
val('MasaIdealVecesMin', bi['veces_carga'][0], 1, 'total / carga útil del caso base, optimista [veces]')
val('MasaIdealVecesMax', bi['veces_carga'][1], 1, 'total / carga útil del caso base, pesimista [veces]')
val('CargaMaxBase', p.PLATAFORMAS[p.nombre_dirigible]['carga_max'] / 1e3, 0, 'carga útil máxima del caso base [t]')
# arreglo con quemador (eta = 0) y rango entre arreglos (auditoría R2-02, autor 2026-09-27)
_bq = bi['quemador']
rango('MasaEjeSobra', bi['eje_sobra_W'][0] / 1e6, bi['eje_sobra_W'][1] / 1e6, 0,
      'eje del motor que sobra frente a los ventiladores durante la quema [MW]')
rango('MasaQuemadorTotal', _bq['total_kg'][0] / 1e3, _bq['total_kg'][1] / 1e3, 0,
      'cajas con quemador (eta 0) y buffer para los ventiladores, sin el quemador [t]')
rango('MasaQuemadorVeces', _bq['veces_carga'][0], _bq['veces_carga'][1], 1, 'lo mismo / carga útil del caso base [veces]')
rango('MasaRangoTotal', _bq['total_kg'][0] / 1e3, bi['total_kg'][1] / 1e3, 0, 'cajas, de quemador a motor [t]')
rango('MasaRangoVeces', _bq['veces_carga'][0], bi['veces_carga'][1], 1, 'cajas, de quemador a motor / carga del caso base [veces]')
val('MasaRangoMin', _bq['total_kg'][0] / 1e3, 0, 'cajas con quemador, extremo liviano [t] (abstract)')
val('MasaRangoVecesMin', _bq['veces_carga'][0], 1, 'cajas con quemador, extremo liviano / carga del caso base [veces] (abstract)')
rango('MasaIdealTurb', bi['turbina_kg'][0] / 1e3, bi['turbina_kg'][1] / 1e3, 0, 'la conversión con turbinas CT7 [t]')
# brecha de masa como FLUJO que cabe en la carga útil (autor, 2026-09-28), "al menos para la configuración propuesta"
_cb = bi['cabe_en_carga']
_cl, _cp = _cb['quemador_liviano'], _cb['motor_pesado']
rango('CabeFlujo', _cp['mdot_kgps'], _cl['mdot_kgps'], 2,
      'flujo que cabe en la carga útil del caso base: motor pesado -- quemador liviano [kg/s]')
rango('CabeVeces', _cl['veces_menos'], _cp['veces_menos'], 1, 'veces menos que el flujo requerido: quemador -- motor')
rango('CabeEscalon', _cp['escalon_kg'] / 1e3, _cl['escalon_kg'] / 1e3, 1,
      'escalón único que ese flujo sostiene (ley m^{3/2}) [t]: motor -- quemador')

# ── Requisito de masa del rechazo y brecha (req_masa_rechazo.py; ruta nueva del autor, 2026-09-28) ──
seccion('Requisito de masa del rechazo por kg/s de H2 y brecha con lo que existe (req_masa_rechazo.py)')
_rq = rutas.leer_datos('req_masa_rechazo')
_r = _rq['requisitos']
_c = {c['clave']: c for c in _rq['candidatos']}
val('ReqPresupuesto', _r['presupuesto_t_por_kgps'], 1, 'presupuesto de las cajas que escalan con el flujo [t por kg/s de H2]')
# ANTES (revisión 2026-10-01, P-04: la tabla 6 ya da 1.45 t; la macro con un decimal no coincidía):
# val('ReqFijas', _r['fijas_t'], 1, 'cajas fijas (regeneración, tanques, compresor), extremo liviano [t]')
val('ReqFijas', _r['fijas_t'], 2, 'cajas fijas (regeneración, tanques, compresor), extremo liviano [t]')
rango('ReqConv', *_r['conversion_t_por_kgps'], 0, 'motor alternativo de H2 [t por kg/s]')
rango('ReqBuf', *_r['buffer_quemador_t_por_kgps'], 0, 'batería que mueve los ventiladores del quemador [t por kg/s]')
val('ReqQueda', _r['queda_rechazo_t_por_kgps'], 1, 'lo que deja el presupuesto al rechazo (quemador, batería liviana) [t por kg/s]')
rango('ReqTurb', *_r['turbina_t_por_kgps'], 1, 'conversión con turbinas CT7 [t por kg/s]')
val('ReqKWkgSinVent', _r['kWkg_pide_f0'], 1, 'kW/kg que pide el rechazo si sus ventiladores no gastaran nada')
val('ReqKWkgFZ', _r['kWkg_pide_f2035'], 0, 'kW/kg que pide el rechazo con la f de FlyZero 2035')
mismo(_r['f_2035'], p.perdida_rad[1], 'f de la curva en el punto de FlyZero 2035 = valFRadFuturo')
val('ReqFMax', _r['f_max'], 2, 'f sobre la cual la batería sola llena el presupuesto [fracción del calor]')
# ventiladores movidos por una turbina de la clase CT7 que solo los mueve (I-2, auditoría 3): el extremo liviano
_tv = _r['turbina_ventiladores']
rango('ReqTurbVentKWkg', *_tv['kWkg'], 1, 'turbinas CT7 (TCDS de EASA), potencia de eje por masa [kW/kg]')
rango('ReqTurbVentT', *_tv['t_por_kgps_f2035'], 1, 'esa turbina para los ventiladores con f de 2035 [t por kg/s]')
rango('ReqKWkgTurb', *_tv['kWkg_pide_f2035'], 1, 'kW/kg que pide el rechazo con esa turbina y f de 2035')
val('ReqKWkgTurbMin', _tv['kWkg_pide_f2035'][0], 1, 'lo mismo, con la turbina más liviana (extremo bajo del requisito)')
val('ReqTurbHdosExtra', _tv['H2_extra_kgps_por_kgps'], 2, 'H2 que quema esa turbina por kg/s del quemador, f de 2035 [kg/s]')
rango('ReqTempVeces', *_r['brecha_temp_ciclo'], 0, 'brecha de las clases que condensan a la temperatura del ciclo [veces]')
# (2026-10-13) el mismo requisito a 1 atm, sin el compresor del escape (auditoría, objeción 5.2): con el calor de 0.99
# el escape va a ~14 atm y su compresor solo ya no cabe (PilaCompEscSinFW); a 1 atm se retiene phi_max
_u1 = _r['una_atm']
mismo(_u1['P_MW'], r6['base_1t_1atm']['P_cond_MW'], 'requisito a 1 atm: la potencia = ReqPcondUnaAtm (EXP-06)', tol=1e-3)
rango('ReqKWkgFfisUnaAtm', *_u1['kWkg_pide_ffis'], 0, 'requisito del rechazo a 1 atm (phi_max, sin compresor del '
      'escape), con batería y el f por caída de presión [kW/kg]')
rango('ReqTempVecesUnaAtm', *_u1['brecha_temp_ciclo'], 0, 'brecha de las clases a la temperatura del ciclo, 1 atm [veces]')
val('ReqPresupuestoUnaAtm', _u1['presupuesto_t_por_kgps'], 1, 'presupuesto para lo que escala, 1 atm [t por kg/s]')
# (2026-10-14, suplementario A2) las tablas S7 y S8 y la fila del PIRT van a 1 atm: lo que falta del mismo requisito
val('ReqKWkgSinVentUnaAtm', _u1['kWkg_pide_f0'], 1, 'kW/kg que pide el rechazo sin ventiladores, 1 atm')
val('ReqKWkgFZUnaAtm', _u1['kWkg_pide_f2035'], 0, 'kW/kg que pide el rechazo con la f de FlyZero 2035, 1 atm')
val('ReqFMaxUnaAtm', _u1['f_max'], 2, 'f sobre la cual la batería sola llena el presupuesto, 1 atm [fracción del calor]')
val('ReqKWkgFZBatMejorUnaAtm', _u1['kWkg_pide_bat_mejor_f2035'], 0, 'kW/kg que pide el rechazo con la mejor batería de '
    '30 s y f de 2035, 1 atm')
val('ReqQAguaUnaAtm', _u1['q_agua_MJkg'], 1, 'calor del ciclo por kg de agua retenida a 1 atm (quemador, phi_max) [MJ/kg]')
mismo(min(min(fl['brecha_rango']) for fl, c in zip(_u1['filas'], _rq['candidatos']) if c['temp_ciclo']),
      _u1['brecha_temp_ciclo'][0], 'tabla S8 a 1 atm: el extremo bajo de la brecha de las clases que condensan = '
      'ReqTempVecesUnaAtm')
if _tv['cumplen'] != ['FZ_sistema_2035', 'FZ_radiador_2026', 'FZ_radiador_2035']:
    sys.exit(f"turbina para los ventiladores: cambian los que cumplen ({_tv['cumplen']}): revisar sec. 5.4 y la tabla 7")
# (2026-10-14) a 1 atm, el caso de la tabla S8 y de la sección 5.1: alcanzan solo las clases de avión de 2035 (100 y
# 160 C); la frase de 5.1 dice "whose heat-rejection temperature of 100--160 C"
_alc1 = [c['clave'] for fl, c in zip(_u1['filas'], _rq['candidatos']) if fl['brecha_rango'][0] < 1.0]
if _alc1 != ['FZ_sistema_2035', 'FZ_radiador_2035']:
    sys.exit(f"a 1 atm cambian las clases que alcanzan el requisito ({_alc1}): revisar la frase de 5.1 (100-160 C)")
if any(c['cumple_turbina'] for c in _rq['candidatos'] if c['temp_ciclo']):
    sys.exit('una clase que condensa a la temperatura del ciclo cumple con la turbina: revisar el resumen y la conclusión 4')
val('ReqACCQ', p.ACC_Q_MWth, 0, 'calor del condensador enfriado por aire de referencia [MW] (Zhang et al. 2020)')
val('ReqACCMasa', _rq['acc_masa_t'], 0, 'masa de ese condensador [t] (tabla S7)')
val('ReqACCT', p.ACC_T_aire_C + p.ACC_ITD_C, 0, 'temperatura a la que condensa [°C]')
val('ReqACCf', _rq['f_acc'], 2, 'potencia de sus ventiladores / calor [fracción]')
# ANTES (revisión 2026-10-01, P-03: la tabla 7 da tres cifras bajo 1 kW/kg; el texto quedaba en 0.14, 0.65 y
# 0.27-0.50 frente a 0.135, 0.653 y 0.267-0.497 de la tabla):
# val('ReqACCKWkg', _c['ACC_sistema']['kWkg'][0], 2, 'kW/kg del sistema completo')
# val('ReqACCHacesKWkg', _c['ACC_haces']['kWkg'][0], 2, 'kW/kg de los haces de tubos solos')
val('ReqACCKWkg', _c['ACC_sistema']['kWkg'][0], 3, 'kW/kg del sistema completo')
val('ReqACCHacesKWkg', _c['ACC_haces']['kWkg'][0], 3, 'kW/kg de los haces de tubos solos')
val('ReqACCPide', _c['ACC_sistema']['pide'], 1, 'kW/kg que pide el rechazo con la f de ese condensador')
# ANTES (P-03, como arriba): rango('ReqRigKWkg', *_c['rigidos']['kWkg'], 2, ...)
rango('ReqRigKWkg', *_c['rigidos']['kWkg'], 3, 'recuperación de agua de los rígidos, equivalente [kW/kg]')
val('ReqRigFav', _c['rigidos']['kWkg_favorable'], 1, 'lo mismo, a potencia máxima y con el agua del metano [kW/kg]')
val('ReqTRocioPres', _r['T_rocio_99_C'], 0, 'rocío del escape a la presión que retiene el 99 % [°C]')
val('ReqQAgua', _r['q_agua_MJkg'], 1, 'calor de este ciclo por kg de agua retenida (quemador) [MJ/kg]')
val('ReqVConsumo', _r['v_consumo_kn'], 0, 'velocidad del consumo de los rígidos (Goodyear 1975, nota ***) [kn]')
val('ReqCandidatos', len(_c), 0, 'candidatos con fuente en la tabla de brecha')
# sensibilidad a la batería (ítem 8): la mejor celda encontrada para 30 s (Amprius SA88, fuentes.AMP26)
_sb = _r['sens_bateria']
val('ReqBatMejorCelda', _sb['kWkg_celda'], 1, 'mejor celda para 30 s (Amprius SA88, pulso x voltaje nominal) [kW/kg]')
val('ReqBatMejorPaquete', _sb['kWkg_paquete'], 1, 'lo mismo a nivel paquete (x 0.66 del X-57) [kW/kg]')
val('ReqKWkgFZBatMejor', _sb['kWkg_pide_f2035'], 0, 'kW/kg que pide el rechazo con esa batería y f de 2035')
if _sb['cumplen'] != ['FZ_radiador_2035']:
    sys.exit(f"sensibilidad de la batería: cambian los que cumplen ({_sb['cumplen']}): revisar la frase de la sec. 5.4")
mismo(_r['presupuesto_t_por_kgps'] * _r['mdot_kgps'] + _r['fijas_t'] + _r['suelta_t'], p.NAVE_BASE['carga_max'] / 1e3,
      'presupuesto x flujo + fijas + la tonelada que se suelta = carga útil (objeción 1.8)')
val('MasaPem', bm['P_pem_W'] / 1e3, 0, 'electrolizador para recuperar 1 t en 24 h [kW]')
rango('MasaTotal', bm['total_kg'][0] / 1e3, bm['total_kg'][1] / 1e3, 1, 'masa total de las cajas negras [t]')
for x in bm['plataformas']:
    c = x['corto']
    rango('MasaCarga' + c, 100 * x['frac_carga'][0], 100 * x['frac_carga'][1], 0, f'masa total / carga útil [%] — {c}')
    for tec, suf in (('a-Si', 'ASi'), ('CIGS', 'Cigs')):
        val('MasaPv' + suf + c, x['pv'][tec]['masa_kg'] / 1e3, 1, f'masa de paneles {tec} [t] — {c}')
        for est, se in _EST:
            val('DiasPv' + suf + c + se, x['pv'][tec]['dias_por_t'][est], 2, f'días por tonelada con paneles {tec}, {est} — {c}')
gr = [x for x in bm['plataformas'] if x['corto'] != 'ZNT']
rango('MasaCargaGrandes', 100 * min(x['frac_carga'][0] for x in gr), 100 * max(x['frac_carga'][1] for x in gr), 0,
      'masa total / carga útil, naves grandes [%]')
for tec, suf in (('a-Si', 'ASi'), ('CIGS', 'Cigs')):
    rango('MasaPv' + suf, min(x['pv'][tec]['masa_kg'] for x in gr) / 1e3, max(x['pv'][tec]['masa_kg'] for x in gr) / 1e3,
          1, f'masa de paneles {tec}, naves grandes [t]')
    for est, se in _EST:
        rango('DiasPv' + suf + se, min(x['pv'][tec]['dias_por_t'][est] for x in gr),
              max(x['pv'][tec]['dias_por_t'][est] for x in gr), 1, f'días por tonelada con paneles {tec}, {est}, naves grandes')
for tec, suf in (('a-Si', 'ASi'), ('CIGS', 'Cigs')):
    ef, kgm2 = p.PV_TEC[tec]
    val('PvEf' + suf, 100 * ef, 0 if abs(100 * ef - round(100 * ef)) < 1e-9 else 1, f'eficiencia {tec} sobre el área del módulo [%]')
    val('PvKgm' + suf, kgm2, 2 if kgm2 < 1 else 1, f'masa por área {tec} [kg/m2]')
rango('ElecCubre', 100 * min(bm['cubre']), 100 * max(bm['cubre']), 0, 'electrólisis cubierta por la propia quema [%]')
rango('ElecCubreAvion', 100 * min(bm['cubre_avion']), 100 * max(bm['cubre_avion']), 0,
      'electrólisis cubierta por la propia quema, con la f de radiadores de avión [%]')
val('ElecPorT', bm['E_el_por_t_kWh'], 0, 'electricidad de electrólisis por tonelada [kWh]')
val('WtCopv', 100 * p.wt_COPV, 1, 'H2 / (H2 + sistema), Tipo IV, estado 2023 del DOE llevado a 350 bar [%]')
val('WtCopvSetecientos', p.kWhkg_tanque_sistema * 3.6e6 / p.Q_comb_H2 * 100, 1, 'lo mismo a 700 bar, como lo da el DOE [%]')
val('PTanqueBar', p.P_tanque_bar, 0, 'presión de los tanques del modelo [bar]')
rango('KWkgPem', p.kWkg_PEM[0], p.kWkg_PEM[1], 1, 'electrolizador PEM [kW/kg] (proxy de pilas PEM)')

# ── Texto nuevo del paper 1 (revisión 2026-09-25) ────────────────────────────
seccion('Paper 1, revisión del autor: potencia por escalón, maniobra, FlyZero, motores')
import minimax as _mmx
from exp05_condensacion_tiempo_real import q_rech as _qr
_q99 = _qr(0.0, p.phi_objetivo)
_grandes = [(n_, d_) for n_, d_ in p.PLATAFORMAS.items() if d_['corto'] != 'ZNT']
_P1 = [_mmx.resolver(d_, 1000.0, p.h_max_seguridad, p.D_lim_oper)['dm_star'] * _q99 / 1e6 for n_, d_ in _grandes]
_P10 = [_mmx.resolver(d_, 10000.0, p.h_max_seguridad, p.D_lim_oper)['dm_star'] * _q99 / 1e9
        for n_, d_ in _grandes if d_['carga_max'] >= 10000.0]
rango('EscPotUno', min(_P1), max(_P1), 0, 'P_rej para soltar 1 t (5 m), naves grandes [MW]')
# D_lim: el mayor que todavía cumple el umbral de 1 m es h_max/(1+s) (la quema limpia baja (1+s) D_lim);
# las cifras de la cota ideal escalan con D_lim^(-1/2) (auditoría R2-06)
_pf = p.PLATAFORMAS[p.nombre_dirigible]
_s = _mmx.resolver(_pf, 1000.0, p.h_max_precision, p.D_lim_oper)['s']
_dmx = p.h_max_precision / (1.0 + _s)
_rmx = _mmx.resolver(_pf, 1000.0, p.h_max_precision, _dmx)
mismo(_mmx.resolver(_pf, 1000.0, p.h_max_precision, p.D_lim_oper)['dm_star'] * _q99 / 1e6,
      rutas.leer_datos('balance_masa')['ideal']['P_rej_eta0_W'] / 1e6, 'P_rej del caso base = la de la tabla 6 con eta 0',
      tol=1e-6)
val('DlimMaxBase', _dmx, 2, 'mayor D_lim que cumple h_max = 1 m, h_max/(1+s) [m] — caso base')
val('DlimMaxQuemaBase', _rmx['dm_star'], 2, 'ṁ* con D_lim = h_max/(1+s) [kg/s] — caso base')
val('DlimMaxPrejBase', _rmx['dm_star'] * _q99 / 1e6, 0, 'P_rej con D_lim = h_max/(1+s) [MW] — caso base')
_rop = _mmx.resolver(_pf, 1000.0, p.h_max_precision, p.D_lim_oper)
val('DlimVceroBase', 3.0 * p.D_lim_oper / _rop['t_anticip'], 3, 'velocidad de bajada al soltar, |v0| = 3 D_lim/t_anticip, D_lim operacional [m/s]')
val('DlimMaxVceroBase', 3.0 * _dmx / _rmx['t_anticip'], 3, 'velocidad de bajada al soltar con D_lim = h_max/(1+s) [m/s]')
rango('EscPotDiez', min(_P10), max(_P10), 1, 'P_rej para soltar 10 t (5 m), naves grandes [GW]')
_op = rutas.leer_datos('figura_operacion')
val('OpArranque', _op['t_release_s'] - _op['t_arranque_s'], 0, 'la unidad de rechazo arranca antes de soltar [s]')
mismo(_op['t_release_s'] - _op['t_arranque_s'], p.n_tau_arranque * p.tau_condensador + _op['t_anticip_s'],
      'arranque previo = n_tau tau + anticipación (figura_operacion vs parametros)', tol=1e-6)
val('OpNTau', p.n_tau_arranque, 0, 'tiempos de respuesta entre la orden de arranque y la quema previa [-]')
val('OpCapArranque', 100.0 * (1.0 - math.exp(-p.n_tau_arranque)), 0,
    'capacidad alcanzada al comenzar la quema previa (primer orden) [%]')
val('OpAnticip', _op['t_anticip_s'], 1, 'quema previa antes de soltar [s]')
val('OpReactiva', _op['t_fin_reactiva_s'] - _op['t_release_s'], 0, 'duración de la quema reactiva [s]')
_r14 = rutas.leer_datos('exp14_calor_por_temperatura')
val('CatorceTFlyHoy', _r14['T_rechazo_flyzero_C'][0], 0, 'temperatura de rechazo de FlyZero, 2026 [°C]')
val('CatorceTFlyFuturo', _r14['T_rechazo_flyzero_C'][1], 0, 'temperatura de rechazo de FlyZero, 2035 [°C]')
rango('EtaMotor', p.eta_motor_H2[0], p.eta_motor_H2[1], 2, 'eficiencia de motores de H2 publicada')
val('CatorceTinaDT', 10.0, 0, 'calentamiento admitido del agua de la tina [°C]')
rango('KWkgCaliente', _r14['kWkg_caliente'][0], _r14['kWkg_caliente'][1], 0, 'intercambiadores de alta temperatura, FlyZero 2026 (en 2035 solo ">20") [kW/kg]')

# ── EXP-15: operación en altura (exp15_altura.py) ────────────────────────────
_a15 = rutas.leer_datos('exp15_altura')
seccion('EXP-15: operación en altura (exp15_altura.json)')
_f0 = _a15['filas'][0]
_fr = next(f for f in _a15['filas'] if f['h_m'] == _a15['H_ref_m'])
_c15 = _a15['cambio_ref']
val('QuinceHRef', _a15['H_ref_m'] / 1e3, 0, 'altura de comparación [km]')
val('QuinceTAireRef', _fr['T_aire_C'], 1, 'aire a esa altura, día caluroso + gradiente ISA [°C]')
val('QuincePhiCero', 100 * _f0['phi'], 0, 'agua retenible a nivel del mar [%]')
val('QuincePhiRef', 100 * _fr['phi'], 0, 'agua retenible a la altura de comparación [%]')
val('QuinceHdosRef', -100 * _c15['H2_por_t'], 1, 'baja del H2 por tonelada a esa altura [%]')
rango('QuinceMasaRef', min(-100 * x for x in _c15['kg_por_kW']), max(-100 * x for x in _c15['kg_por_kW']), 0,
      'baja de la masa de rechazo por kW a esa altura, 2026--2035 [%]')
val('QuinceSustRef', -100 * _c15['sustentacion'], 0, 'pérdida de sustentación de la misma nave a esa altura [%]')

# ── EXP-16: balance neto por tonelada, materia y energía (exp16_balance_neto.py) ──
_b16 = rutas.leer_datos('exp16_balance_neto')
seccion('EXP-16: balance neto por tonelada (exp16_balance_neto.json)')
_c16 = {f['caso'] + (f['corto'] or ''): f for f in _b16['filas']}
_C16 = [f for f in _b16['filas'] if f['caso'] == 'C']
rango('DieciseisEComp', _b16['e_comp_kWhkg'][0], _b16['e_comp_kWhkg'][1], 1, 'electricidad de compresión 20 -> 350 bar, típica DOE -- medida [kWh/kg]')
for _k, _suf in (('A', 'Ideal'), ('B', 'Atm')):
    _f = _c16[_k]
    val('DieciseisQuim' + _suf, _f['E_quimica'], 0, f'energía química quemada por tonelada, caso {_k} [kWh/t]')
    val('DieciseisTrabajo' + _suf, _f['W_generado'], 0, f'trabajo generado por tonelada, caso {_k} [kWh/t]')
    rango('DieciseisVent' + _suf, min(_f['E_ventiladores']), max(_f['E_ventiladores']), 0, f'ventiladores por tonelada, caso {_k} [kWh/t]')
    val('DieciseisElec' + _suf, _f['E_electrolisis'], 0, f'electrólisis por tonelada, caso {_k} [kWh/t]')
    rango('DieciseisComp' + _suf, min(_f['E_compresion']), max(_f['E_compresion']), 0, f'compresión por tonelada, caso {_k} [kWh/t]')
    rango('DieciseisNeta' + _suf, _f['E_externa_neta'][0], _f['E_externa_neta'][1], 0, f'electricidad externa neta por tonelada, caso {_k} [kWh/t]')
val('DieciseisReposAtm', _c16['B']['materia']['H2_reposicion'], 1, 'H2 de reposición por tonelada, 1 atm [kg/t]')
# con quemador (ventiladores con la batería, recargada desde afuera) y días de paneles con la NETA (hallazgo M7):
# del motor (extremo bajo) al quemador (extremo alto), 1 atm, naves grandes
# R-03 (revisión 2026-10-01): las de arriba van con la f por caída de presión; con la de radiadores de avión, cota alta
rango('DieciseisNetaAtmAvion', *_c16['B']['avion']['E_externa_neta'], 0,
      'electricidad externa neta por tonelada, 1 atm, con la f de radiadores de avión (cota alta) [kWh/t]')
rango('DieciseisNetaQuemAtm', *_c16['B']['E_externa_neta_quemador'], 0,
      'electricidad externa neta por tonelada con quemador, 1 atm [kWh/t]')
for _tec, _sf in (('a-Si', 'ASi'), ('CIGS', 'Cigs')):
    for _est, _se in (('invierno', 'Inv'), ('verano', 'Ver')):
        _k16 = [k for k in _c16['B']['dias_paneles'] if k.endswith(f'|{_tec}|{_est}')]
        rango('DiasNeta' + _sf + _se, min(_c16['B']['dias_paneles'][k][0] for k in _k16),
              max(_c16['B']['dias_paneles_quemador'][k][1] for k in _k16), 1,
              f'días por tonelada con paneles {_tec}, {_est}, electricidad neta de motor a quemador, 1 atm, naves grandes')
# ── M-11 (auditoría 3): entradas de escenario que el .tex escribía a mano; cada macro sale de la constante que usa
# el script y se comprueba contra los demás que usan el mismo nivel ──
import exp16_balance_neto as _e16m
import exp12_helio as _e12m
import exp14_calor_por_temperatura as _e14m
_caps9 = {x['P_cap_MW'] for x in rutas.leer_datos('exp09_condensacion_dinamica')['casos_base'] if x.get('P_cap_MW')}
_caps11 = {b['P_instalada_MW'] for b in rutas.leer_datos('exp11_balance_energia')['buffer']}
_pdiez, _pcien = _e16m.P_ESCALON_MW, 100.0
for _P in (_pdiez, _pcien):
    if _P not in _caps9 or _P not in _caps11:
        sys.exit(f'M-11: el nivel de {_P:g} MW no está en EXP-09 o EXP-11: revisar las macros PotDiez/PotCien')
mismo(bm['P_ref_W'] / 1e6, _pdiez, 'potencia de referencia de balance_masa = escalones de EXP-16 (10 MW)')
seccion('Entradas de escenario citadas en el texto (M-11, auditoría 3)')
val('PotDiez', _pdiez, 0, 'potencia instalada de referencia: escalones (EXP-16), buffer (EXP-11), EXP-09 [MW]')
val('PotCien', _pcien, 0, 'potencia instalada alta de EXP-09 y EXP-11 [MW]')
# sección 5 reescrita (autor, 2026-10-12): lo que pesaría el caso de 100 MW de la Fig. 3 con los haces de ACC
val('FrenoHXCien', _fr_haces_t_por_MW * _pcien, 0, 'haces de ACC para los 100 MW de la Fig. 3 [t]')
val('TRegenH', p.t_regeneracion_h, 0, 'tiempo para regenerar el H2 de una tonelada [h]')
if 20.0 not in _e12m.P_ALMACEN_MPA:
    sys.exit('M-11: EXP-12 ya no evalúa 20 MPa (macros DoceWVeinte, DoceFactorCalor)')
val('DocePRefMPa', 20.0, 0, 'presión de almacenamiento del He de referencia (EXP-12, macros ...Veinte) [MPa]')
if 5.0 not in _e14m.PRESIONES_ATM:
    sys.exit('M-11: EXP-14 ya no evalúa 5 atm (macros ...Cinco)')
val('CatorcePCinco', 5.0, 0, 'presión del escape de la discusión (EXP-14, macros ...Cinco) [atm]')
rango('DieciseisNetaEsc', min(f['E_externa_neta'][0] for f in _C16), max(f['E_externa_neta'][1] for f in _C16), 0,
      'electricidad externa neta por tonelada, escalones de 10 MW, naves grandes [kWh/t]')
rango('DieciseisVentEsc', min(min(f['E_ventiladores']) for f in _C16), max(max(f['E_ventiladores']) for f in _C16), 0,
      'ventiladores por tonelada con arranque previo, escalones de 10 MW [kWh/t]')
rango('DieciseisReposEsc', min(f['materia']['H2_reposicion'] for f in _C16), max(f['materia']['H2_reposicion'] for f in _C16), 1,
      'H2 de reposición por tonelada, escalones de 10 MW [kg/t]')
rango('DieciseisTierraEsc', min(f['E_tierra_reposicion'] for f in _C16), max(f['E_tierra_reposicion'] for f in _C16), 0,
      'energía en tierra del H2 de reposición, escalones de 10 MW [kWh/t]')
val('DieciseisCompFlujo', _b16['compresor']['flujo_kgph'], 2, 'flujo del compresor (= electrolizador) para regenerar 1 t por día [kg/h]')
rango('DieciseisCompPot', _b16['compresor']['P_kW'][0], _b16['compresor']['P_kW'][1], 1, 'potencia del compresor [kW]')
rango('DieciseisCompMasa', _b16['compresor']['masa_cota_sup_kg'][0] / 1e3, _b16['compresor']['masa_cota_sup_kg'][1] / 1e3, 2,
      'masa del compresor eléctrico de referencia (Hydro-Pac C06-15), cota superior [t]')
rango('DieciseisMasaConComp', _b16['masa_con_compresor_t'][0], _b16['masa_con_compresor_t'][1], 1,
      'masa de las cajas negras con el compresor (cota superior) [t]')

# ── VAL-02: condensación contra experimentos H2-aire (val02_condensacion.py) ──
v2 = rutas.leer_datos('val02_condensacion')
seccion('VAL-02: Ec. 5 contra Zakhvatkin et al. 2022 (val02_condensacion.json)')
val('ValDosCurvas', 100 * v2['dif_max_curvas'], 2, 'diferencia máxima con su figura 4 [puntos de %]')
_c = v2['experimentos'][:5]
rango('ValDosDesvio', 100 * min(abs(x['esperada_g'] / x['esperada_zak_recalc_g'] - 1) for x in _c),
      100 * max(abs(x['esperada_g'] / x['esperada_zak_recalc_g'] - 1) for x in _c), 1,
      '|agua esperada: este modelo vs su apéndice A| (experimentos 1-5; el modelo queda por debajo) [%]')
# NUEVA (revisión 2026-09-30, I-01): cuánto difiere SU tabla 1 de SU apéndice A recalculado con sus flujos (mismos
# experimentos 1-5): explica por qué la referencia del texto y de la figura es el recálculo y no la tabla
rango('ValDosTablaApendice', 100 * min(abs(x['esperada_zak_g'] / x['esperada_zak_recalc_g'] - 1) for x in _c),
      100 * max(abs(x['esperada_zak_g'] / x['esperada_zak_recalc_g'] - 1) for x in _c), 0,
      '|su tabla 1 vs su apéndice A recalculado con sus flujos| (experimentos 1-5) [%]')
rango('ValDosMedido', 100 * min(x['medido_sobre_modelo'] for x in v2['experimentos']),
      100 * max(x['medido_sobre_modelo'] for x in v2['experimentos']), 0, 'agua recogida / cota de equilibrio [%]')
rango('ValDosT', min(x['T_amb_C'] for x in v2['experimentos']), max(x['T_amb_C'] for x in v2['experimentos']), 0,
      'temperatura ambiente de los experimentos [°C]')
val('ValDosN', len(v2['experimentos']), 0, 'experimentos comparados')
import fuentes as _fz
_zp = _fz.FUENTES['ZAK22']['datos']['recoleccion_sobre_su_modelo'][0]
rango('ValDosZakPropio', 100 * _zp[0], 100 * _zp[1], 0, 'agua recogida / SU modelo, según los autores [%] (fuentes.ZAK22)')
from exp05_condensacion_tiempo_real import p_sat as _ps2
_e5, _e6 = v2['experimentos'][4], v2['experimentos'][5]
val('ValDosCincoRecalc', _e5['esperada_zak_recalc_g'], 0, 'exp. 5: su apéndice A recalculado [g]')
val('ValDosCincoTabla', _e5['esperada_zak_g'], 0, 'exp. 5: valor esperado de su tabla [g]')
val('ValDosSeisT', _e6['T_amb_C'], 1, 'exp. 6: temperatura ambiente [°C]')
val('ValDosSeisAH', _ps2(_e6['T_amb_C']) / (p.R_vapor * (_e6['T_amb_C'] + 273.15)) * 1e3, 1,
    'exp. 6: humedad de saturación que corresponde a su temperatura (IF97) [g/m3]')

# ── Admisión de aire al casco (req_aire_casco.py) ─────────────────────────────
ra = rutas.leer_datos('req_aire_casco')
seccion('Admisión de aire al casco (req_aire_casco.json)')
_fa = {x['corto']: x for x in ra['filas']}
_ga = [x for x in ra['filas'] if x['corto'] != 'ZNT']
val('AireCascoBase', _fa[BASE]['aire_1m_m3s'], 0, 'aire al casco, caso base (caso base, 1 t, 1 m) [m3/s]')
rango('AireCascoUno', min(x['aire_1m_m3s'] for x in _ga), max(x['aire_1m_m3s'] for x in _ga), 0,
      'aire al casco, 1 t, 1 m, naves grandes [m3/s]')
rango('AireCascoCinco', min(x['aire_5m_m3s'] for x in _ga), max(x['aire_5m_m3s'] for x in _ga), 0,
      'aire al casco, 1 t, 5 m, naves grandes [m3/s]')
val('AireCascoZNT', _fa['ZNT']['aire_1m_m3s'], 0, 'aire al casco, Zeppelin NT, 1 t, 1 m [m3/s]')
val('AireCascoEscalon', _fa[BASE]['aire_escalon_m3s'], 2, 'aire al casco en los escalones de 10 MW [m3/s]')
val('AireCombBase', _fa[BASE]['aire_comb_1m_m3s'], 0, 'aire de combustión, caso base [m3/s]')
val('AirePorT', ra['volumen_por_t_m3'], 0, 'aire que entra al casco por tonelada liberada [m3]')

# ── Recuperación de agua histórica (fuentes.GY75, Goodyear 1975, p. 20, ítem 18) ──
import fuentes as _fu
seccion('Recuperación de agua del escape en los rígidos históricos (fuentes.GY75)')
_ag, _br = _fu.FUENTES['GY75']['datos']['agua_lb'][0], _fu.FUENTES['GY75']['datos']['bruto_lb'][0]
rango('AguaHistMasa', min(_ag.values()) * p.LB / 1e3, max(_ag.values()) * p.LB / 1e3, 1,
      'masa de los equipos de recuperación de agua (Los Angeles, Graf Zeppelin, Macon, Hindenburg) [t]')
rango('AguaHistPct', min(100 * _ag[k] / _br[k] for k in _ag), max(100 * _ag[k] / _br[k] for k in _ag), 1,
      'masa de esos equipos / sustentación bruta [%]')

# ── Derivaciones del paper 1 (derivaciones.py) ────────────────────────────────
seccion('Derivaciones del paper 1 (derivaciones.json): entradas y pasos intermedios')
for nombre, valor, dec, com in rutas.leer_datos('derivaciones')['macros']:
    if isinstance(valor, (list, tuple)):
        rango(nombre, valor[0], valor[1], dec, com)
    else:
        val(nombre, valor, dec, com)
alias_base('DerCasco')

# ── VER-01: verificación interna del subsistema (verificacion.py) ──
ver = rutas.leer_datos('verificacion')
if ver['n_fallas']:
    sys.exit(f'verificacion.py tiene {ver["n_fallas"]} chequeos en falla: no se exporta numeros.tex')
# la verificación debe ser de ESTOS datos (M-6, auditoría 3): si un script cambió resultados/datos/ después, o
# verificacion.py se cayó antes de escribir, el sello no coincide y no se exporta
if ver.get('sello_datos') != rutas.sello_datos():
    sys.exit('verificacion.json no corresponde a los datos actuales: correr verificacion.py antes de exportar')
seccion('VER-01: verificación interna del subsistema (verificacion.json)')
val('VerTotal', ver['n_total'], 0, 'chequeos de verificación')
val('VerOk',    ver['n_ok'], 0,    'chequeos que pasan')
# R-09 (revisión 2026-10-01): los chequeos por tipo, no solo el total
val('VerCond',   ver['n_condiciones'], 0, 'chequeos que son condiciones (desigualdades, signos, cotas)')
val('VerNula',   ver['n_diferencia_nula'], 0, 'chequeos de diferencia nula (identidades entre scripts, regresiones)')
val('VerTol',    ver['n_con_tolerancia'], 0, 'chequeos con tolerancia (forma cerrada, solución exacta, convergencia)')
ch = ver['chequeos']

# ── Auditoría 2 de la copia del 2026-10-06 ──────────
seccion('Auditoría 2 de la copia del 2026-10-06: H2 comprimido, haces, devolución, hielo, ventear y detalles')
from exp05_condensacion_tiempo_real import q_rech as _qA
from exp16_balance_neto import e_comp as _ecA
_hc = r12['h2_comprimido']
# M1: H2 comprimido a 35 MPa (los tanques del ciclo)
val('DoceHcW', _hc['W_iso_MJkg'], 1, 'H2 comprimido a 35 MPa: trabajo mínimo isotérmico, gas real (Leachman) [MJ/kg]')
val('DoceHcGJ', _hc['calor_GJ_por_t'], 2, 'H2 comprimido a 35 MPa: trabajo (= calor) por tonelada liberada [GJ/t]')
val('DoceHcPot', _hc['potencia_MW'], 0, 'H2 comprimido: potencia de compresión de golpe, caso base [MW]')
val('DoceHcPotRec', _hc['calor_GJ_por_t'] * 1e9 / _trecD / 1e6, 1, 'H2 comprimido: potencia en el tiempo de recuperación del ciclo con 10 MW [MW]')
val('DoceHcLiber', _hc['m3_almacen_por_t'], 1, 'H2 comprimido a 35 MPa: volumen guardado por tonelada [m3]')
# (2026-10-14, objeción 3.1) el mismo balance en una entrega: comprimir guarda la sustentación de cada tonelada
# entregada, así que el compuesto crece con toda la carga (la entrega de referencia de EXP-22: n sueltas de s t)
_En22 = rutas.leer_datos('exp22_energia_potencia')['E']
val('DoceHcCompEntrega', _hc['reserva_compuesto_t'] * _En22['n'] * _En22['s_t'], 0,
    'H2 comprimido: compuesto de los tanques para la entrega de referencia del LCA60T (36 sueltas de 1 t), cota '
    'inferior [t]')
# (2026-10-14, objeción 3.21) el compuesto, proporcional a P V, casi no depende de la presión de almacenamiento
_alm = r12['almacen']
rango('DoceHeCompRango', min(a['tanque_compuesto_t_por_t'] for a in _alm),
      max(a['tanque_compuesto_t_por_t'] for a in _alm), 2, 'He: compuesto de los tanques por tonelada en las presiones '
      'evaluadas [t por t]')
rango('DoceAlmacenPRango', min(a['P_MPa'] for a in _alm), max(a['P_MPa'] for a in _alm), 0,
      'He: presiones de almacenamiento evaluadas [MPa]')
# M3: tiempo para devolver el resto, al ritmo del electrolizador de la tabla 6; energía a 1 atm
val('DoceTRestoH', r12['rutas'][0]['t_resto_h'], 0, 'ciclo (phi = 0.99): horas para devolver el resto con el electrolizador de la tabla 6 [h]')
mismo(r12['rutas'][1]['t_resto_h'], p.t_regeneracion_h, 'el electrolizador regenera la tonelada a 1 atm en t_regeneracion_h')
mismo(r12['electrolizador_kgph'], _b16['compresor']['flujo_kgph'], 'ritmo del electrolizador = flujo del compresor (EXP-16)', tol=1e-6)
_g1 = r12['rutas'][1]
_Rg1 = _g1['kg_por_t'] * (1.0 - _g1['perdida_frac'])
rango('DoceRecupHdosUnaAtm', _Rg1 * (p.E_H2_electrolisis + min(_ecA())) / 1e3, _Rg1 * (p.E_H2_electrolisis + max(_ecA())) / 1e3, 1,
      'energía para devolver la sustentación de 1 t a 1 atm: electrólisis y compresión [MWh]')
# M2: el sistema ACC entero frente a los haces
_cA = {c['clave']: c for c in _rq['candidatos']}
_kwACC = _cA['ACC_sistema']['kWkg'][0]
val('FrenoACCEnteroDiez', p.P_ref_masa / 1e3 / _kwACC / 1e3, 0, 'sistema ACC entero para 10 MW [t]')
val('FrenoPLimEntero', p.NAVE_BASE['carga_max'] * _kwACC / 1e3, 1, 'potencia con la que el ACC entero pesa la carga del LCA60T [MW]')
val('XNdosEscape', 100 * (1.0 - z['x_H2O']), 0, 'fracción molar de N2 y Ar en el escape [%]')
val('CatorceFracCaliente', 100 * (1.0 - b0['frac_fria']), 0, 'fracción del calor por encima del rocío, eta = 0 [%]')
# m10: unidad en dos tramos para la cota ideal del caso base
_PrA = r12['base']['H2']['potencia_MW']
_Pband = _PrA * b0['frac_fria']
_Pcal = _PrA - _Pband
val('TramoBandaMW', _Pband, 0, 'cota ideal del caso base: calor de la banda [MW]')
val('TramoCalienteMW', _Pcal, 0, 'cota ideal del caso base: calor del gas caliente, sobre el rocío [MW]')
val('TramoBandaT', _Pband * 1e3 / _cA['ACC_haces']['kWkg'][0] / 1e3, 0, 'banda con haces de ACC [t]')
rango('TramoCalienteT', _Pcal * 1e3 / max(p.kWkg_radiador) / 1e3, _Pcal * 1e3 / min(p.kWkg_radiador) / 1e3, 0,
      'gas caliente con radiadores de avión 2035-2026 [t]')
rango('TramoTotalT', (_Pband / _cA['ACC_haces']['kWkg'][0] + _Pcal / max(p.kWkg_radiador)),
      (_Pband / _cA['ACC_haces']['kWkg'][0] + _Pcal / min(p.kWkg_radiador)), 0, 'unidad en dos tramos, total [t]')
# M6: hielo para la banda
_Qb = r14['tina']['Q_frio_por_t_GJ'] * 1e9
_hielo = (_Qb / (p.L_fusion_hielo + p.cp_agua_liq * p.T_aire_exterior) / 1e3, _Qb / p.L_fusion_hielo / 1e3)
val('HieloBandaGJ', _Qb / 1e9, 2, 'calor de la banda por tonelada liberada [GJ/t]')
rango('HieloT', _hielo[0], _hielo[1], 1, 'hielo que absorbe la banda por tonelada: fundido y calentado al ambiente -- solo fundido [t/t]')
rango('HieloCargaT', _hielo[0] * p.NAVE_BASE['carga_max'] / 1e3, _hielo[1] * p.NAVE_BASE['carga_max'] / 1e3, 0,
      'hielo para toda la carga del LCA60T [t]')
# M7: frente a ventear
_vent = r12['rutas'][2]['perdido_kg_por_t']
_Event = _vent * p.E_H2_electrolisis / 1e3
val('VentAhorroKg', _vent - r12['rutas'][1]['perdido_kg_por_t'], 1, 'H2 que el ciclo ahorra frente a ventear, 1 atm [kg/t]')
val('VentReponerMWh', _Event, 2, 'electrólisis para reponer el H2 venteado por tonelada [MWh]')
# (2026-10-13) CicloExtMWh y VentAhorroPct salen de EXP-22 (la fuente de los números de energía; antes, EXP-16): se
# definen después del bloque de EXP-21/22
# m4: motor al ritmo de 10 MW
_dm10 = p.P_ref_masa / _qA(0.0, p.phi_objetivo)
_eje = [e_ * p.Q_comb_H2 * _dm10 for e_ in p.eta_motor_H2]
rango('MotorEjeDiez', min(_eje) / 1e6, max(_eje) / 1e6, 1, 'eje del motor al ritmo de 10 MW de rechazo [MW]')
rango('MotorMasaDiez', min(_eje) / 1e3 / max(bm['kWkg_motor']) / 1e3, max(_eje) / 1e3 / min(bm['kWkg_motor']) / 1e3, 0,
      'masa del motor al ritmo de 10 MW [t]')
# m2: la inercia depende de la densidad absoluta; día de diseño de 30 °C
_rho30 = p.P_atmos / (p.R_aire * (p.T_aire_exterior + 273.15))
val('InerciaTreintaPct', 100 * (1 - _rho30 / p.rho_aire_std), 0, 'inercia vertical más baja el día de 30 °C [%]')
val('TasaTreintaPct', 100 * ((p.rho_aire_std / _rho30) ** 0.5 - 1), 1, 'tasa instalada más alta el día de 30 °C [%]')
# m3: un esferoide de 200 x 50 m
_Vesf = math.pi / 6.0 * _pf['L'] * _pf['D'] ** 2
val('VolEsferoideFW', round(_Vesf, -3), 0, 'volumen del esferoide de las dimensiones publicadas del LCA60T [m3]')
val('TasaEsferoidePct', 100 * (1 - (_pf['V_casco'] / _Vesf) ** 0.5), 0, 'tasa instalada más baja con ese volumen [%]')
# m13: abortar después de empezar la quema previa
_kr = p.rho_aire_std / (p.P_atmos / (p.R_H2 * p.T1))
_pesada = _mmx.L_eff * _rop['dm_star'] * _rop['t_anticip']
val('AbortoKg', round(_pesada, -1), 0, 'nave pesada si la liberación falla tras la quema previa, caso base [kg]')
val('AbortoHdosKg', _pesada / _kr, 0, 'H2 a devolver desde los tanques [kg]')
val('AbortoS', _pesada / _kr / _rop['dm_star'], 0, 'tiempo para devolverlo al ritmo instalado [s]')
# M4: la maniobra limpia usa una fracción del umbral de precisión
# pregunta del autor (2026-10-07): las naves de helio como variantes nocionales de H2 del mismo casco. Con la misma masa,
# el gas ocupa el volumen que la flota: V_H2/V_He = (rho_a - rho_He)/(rho_a - rho_H2); el resto del casco es aire
_ra = p.rho_aire_std
val('HeHdosVolRazon', (_ra - r12['rho_He']) / (_ra - r12['rho_H2']), 2, 'volumen de H2 / volumen de He para la misma sustentación')
val('HeHdosVolPct', 100 * (1 - (_ra - r12['rho_He']) / (_ra - r12['rho_H2'])), 0, 'gas que el H2 ahorra frente al He con la misma sustentación [%]')
val('DlimFraccionPct', 100 * (1.0 + _s) * p.D_lim_oper / p.h_max_precision, 0, '(1+s) D_lim frente a h_max de precisión [%]')

# ── Figura de los esferoides (2026-10-08) ──
seccion('Figura de los esferoides: cascos como esferoides prolatos y k de Lamb (fig_esferoides.py)')
_fe = rutas.leer_datos('fig_esferoides')
rango('EsfRazonOtras', *_fe['razon_V_otras'], 2, 'volumen del esferoide / V_casco, LZ 129, ZR-3 y NT')
val('EsfRazonFW', _fe['naves']['FW']['razon_V'], 2, 'volumen del esferoide / V_casco, LCA60T')
# (2026-10-15, autor) el LCA60T se toma con su volumen de helio publicado, cota inferior de su casco (sección 3.5): el
# esferoide de sus medidas encierra EsfRazonFW veces más, y con la elasticidad de ṁ* respecto de la inercia sus
# requisitos dinámicos bajarían este porcentaje
_bajaFW = 100.0 * (1.0 - _fe['naves']['FW']['razon_V'] ** _ea['M_I'])
mismo(_bajaFW, 100.0 * (1.0 - _fe['naves']['FW']['razon_V'] ** -0.5),
      'baja de los requisitos dinámicos del LCA60T con el esferoide (elasticidad vs -1/2)', tol=0.02)
val('EsfFWBajaPct', _bajaFW, 0, 'cuánto bajarían los requisitos dinámicos del LCA60T con el esferoide de sus medidas [%]')
val('EsfFWEqL', _fe['naves']['FW']['L_equiv_m'], 0, 'largo del esferoide del LCA60T con su V_casco y su esbeltez [m]')
val('EsfFWEqD', _fe['naves']['FW']['D_equiv_m'], 0, 'diámetro de ese esferoide [m]')

# ── Respuestas a críticas externas (autor, 2026-10-11): mismas fuentes y funciones que la tabla 7 (tablas_articulo.helio)
seccion('Críticas externas: 1 atm sin compresor, compresor y motor a 10 MW, potencia neta con el freno, inercia')
from exp05_condensacion_tiempo_real import q_rech as _qrR
from balance_masa import f_caida_presion as _fcpR
_c99, _c1 = r12['rutas'][0], r12['rutas'][1]
if '1 atm' not in _c1['ruta']:
    sys.exit('rutas de exp12: la segunda ya no es la de 1 atm')
_tr99 = 1000.0 / _c99['L'] * _c99['q_Jkg'] / p.P_ref_masa          # recuperación de 1 t con 10 MW, phi = 0.99 [s]
_tr1 = 1000.0 / _c1['L'] * _c1['q_Jkg'] / p.P_ref_masa              # lo mismo a 1 atm, phi_max, sin compresor [s]
val('FrenoTminUnaAtm', _tr1 / 60, 1, 'recuperación de 1 t con 10 MW a 1 atm, sin compresor (phi_max) [min]')
val('FrenoTminUnaAtmPct', 100 * (_tr1 / _tr99 - 1), 0, 'cuánto más tarda a 1 atm que con phi = 0.99 [%]')
val('DocePerdidoCiclo', _c99['perdido_kg_por_t'], 1, 'H2 perdido por tonelada, phi = 0.99 [kg]')
val('DocePerdidoUnaAtm', _c1['perdido_kg_por_t'], 1, 'H2 perdido por tonelada a 1 atm, phi_max [kg]')
_dm10 = _c99['kg_por_t'] / _tr99                                     # ritmo de quema con 10 MW [kg/s]
val('FrenoDmDiez', _dm10, 3, 'ritmo de quema con 10 MW de rechazo y el dispositivo de freno [kg/s]')
val('FrenoFlotDiez', _dm10 * _c99['L'], 1, 'con 10 MW, flotabilidad que la quema quita por segundo: L x tasa [kg/s]')
val('FrenoCompDiez', _dm10 * r6['W_comp_MJkg'], 1, 'compresor del escape a la presión de phi objetivo, al ritmo de 10 MW [MW]')
# masa de ese compresor con el paquete del SGT-800 (fuentes.SGT800): todo el paquete se atribuye a su compresor y se
# escala por el trabajo isentrópico de una etapa, sin enfriamiento intermedio (el mayor); la eficiencia se cancela.
# Entra al compresor el inerte con el vapor saturado a T_salida y 1 atm; el SGT-800 comprime aire ISO (escape - combustible)
from exp05_condensacion_tiempo_real import p_sat as _psC, N_INERTE as _NIC, AIRE_ESTEQ as _AEC
_nC = _NIC / p.M_H2                                                     # inerte [mol por kg de H2]
_nvC = _nC * _psC(p.T_salida_cond) / (p.P_atmos - _psC(p.T_salida_cond))    # vapor saturado [mol por kg de H2]
_CC = (_AEC - p.ratio_O2_H2) * p.cp_N2 + _nvC * p.M_H2O * p.cp_vapor    # capacidad calorífica [J/K por kg de H2]
_kC = (_nC + _nvC) * p.R_univ / _CC                                     # R / cp de la mezcla
_wsC = _CC * (p.T_salida_cond + 273.15) * (r6['P_req_atm'] ** _kC - 1.0)    # trabajo isentrópico [J/kg H2]
_s8, _ka = p.SGT800, (p.gamma_aire - 1.0) / p.gamma_aire
_wS8 = [(me - hr * P / 3.6 / p.PCI_CH4 * 1e3) * p.cp_aire * (p.T_ISO_C + 273.15) * (pr ** _ka - 1.0)
        for P, pr, me, hr in zip(_s8['P_MWe'], _s8['PR'], _s8['m_escape_kgps'], _s8['heat_rate_kJkWh'])]   # [W]
_mC = [_s8['masa_kg'] * _dm10 * _wsC / w / 1e3 for w in _wS8]          # [t]
_mCb = [_s8['masa_kg'] * r6['base_1t_h1']['dm_star_kgps'] * _wsC / w / 1e3 for w in _wS8]
val('SGTMasa', _s8['masa_kg'] / 1e3, 0, 'paquete clásico del SGT-800, con el generador [t]')
rango('SGTPR', min(_s8['PR']), max(_s8['PR']), 1, 'relación de presión del compresor del SGT-800')
rango('FrenoCompMasaDiez', min(_mC), max(_mC), 1,
      'compresor del escape al ritmo de 10 MW, escalado por trabajo de compresión con el paquete del SGT-800 [t]')
rango('FrenoCompMasaDobleDiez', 2 * min(_mC), 2 * max(_mC), 1, 'lo mismo con el doble de kg/kW (máquina chica) [t]')
rango('CompMasaBase', min(_mCb), max(_mCb), 0, 'lo mismo para el compresor de la cota ideal, caso base [t]')
rango('FrenoMotorDiez', min(_r['conversion_t_por_kgps']) * _dm10, max(_r['conversion_t_por_kgps']) * _dm10, 1,
      'motor alternativo de H2 al ritmo de 10 MW [t]')
# (2026-10-12) a 14 atm el motor también paga el compresor del escape, y su trabajo sale como calor: neto por kg =
# eta PCI - f (q0 - eta PCI + W) - W (ANTES sin W: sobrestimaba el neto en ~0.5 MW a 10 MW)
_WcN = r6['W_comp_MJkg'] * 1e6
_netaR = [eta * p.Q_comb_H2 * _dm10 - fv * (_c99['q_Jkg'] - eta * p.Q_comb_H2 + _WcN) * _dm10 - _WcN * _dm10
          for eta in p.eta_motor_H2 for fv in (min(_fcpR()), max(p.perdida_rad))]
rango('DoceNetaCiclo', min(_netaR) / 1e6, max(_netaR) / 1e6, 1,
      'potencia neta del ciclo con motor hacia la planta, con el freno y 10 MW, ya descontados los ventiladores y el '
      'compresor del escape [MW]')
rango('ZetaInercia', p.zeta / math.sqrt(vma['medida_sobre_analitica']['vuelo_libre']),
      p.zeta / math.sqrt(vma['medida_sobre_analitica']['restringido']), 2,
      'amortiguamiento del lazo con la masa vertical medida (zeta / raiz(f), f de Yamasaki & Goto 2003)')

# ── Presupuesto de masa apilado (autor, 2026-10-11): la figura de la sección 5.2 (fig_presupuesto_masa.py) ──
seccion('Presupuesto de masa sin y con freno, las cuatro naves (fig_presupuesto_masa.py; la tabla 5 pasó al suplemento)')
_ap = rutas.leer_datos('fig_presupuesto_masa')['apilado']
_apn = {nv['corto']: nv for nv in _ap['naves']}
rango('PilaSinNaves', min(sum(nv['sin']['t']) for nv in _ap['naves']), max(sum(nv['sin']['t']) for nv in _ap['naves']), 0,
      'sin freno, equipo y la tonelada (clase más liviana, haces de ACC), de la nave más liviana a la más pesada [t]')
val('PilaSinFW', sum(_apn['FW']['sin']['t']), 0, 'lo mismo en el LCA60T [t]')
val('PilaSinHXFW', _apn['FW']['sin']['t'][0], 0, 'de eso, el rechazo con los haces de ACC en el LCA60T [t]')
# (2026-10-13) a 14 atm, sin freno, el compresor del escape solo ya pesa más que la carga útil (auditoría, objeción 5.2)
val('PilaCompEscSinFW', _apn['FW']['sin']['t'][1], 0, 'sin freno, LCA60T: compresor del escape (14 atm) a la tasa de '
    'tiempo real [t]')
val('PilaCompEscSinVeces', _apn['FW']['sin']['t'][1] / _apn['FW']['carga_t'], 1, 'lo mismo / la carga útil [veces]')
val('PilaRazonCondMotor', (_apn['FW']['sin']['t'][0] + _apn['FW']['sin']['t'][1]) / _apn['FW']['sin']['t'][2], 0,
    'sin freno, LCA60T: haces + compresor del escape (lo que existe para condensar) / motor (igual con cualquier tasa)')
_fjk = rutas.leer_datos('req_masa_rechazo')['requisitos']['fijas_kg']
rango('FijaPem', _fjk['regeneracion'][0] / 1e3, _fjk['regeneracion'][1] / 1e3, 2,
      'caja fija: electrolizador PEM para regenerar 1 t en 24 h [t por t]')
val('PilaConDiez', sum(_apn['FW']['con']['t']), 0, 'con freno y 10 MW, equipo y la tonelada (igual en las cuatro naves) [t]')
val('PilaTporMW', _ap['t_por_MW'], 1, 'con freno: masa de lo que escala con el rechazo (haces, motor, compresor) [t por MW]')
for _c, _n in (('FW', 'FW'), ('Hind', 'Hind'), ('ZR', 'ZR')):
    val('PilaMax' + _n, _apn[_c]['max']['P_MW'], 1, f'con freno, el rechazo más grande que cabe en la carga útil, {_c} [MW]')
    val('PilaMax' + _n + 'Min', _apn[_c]['max']['t_rec_min'], 1 if _apn[_c]['max']['t_rec_min'] < 10 else 0,
        f'con él, la recuperación de 1 t, {_c} [min]')
val('PilaMaxACCFW', _apn['FW']['max_acc']['P_MW'], 1, 'lo mismo con el sistema ACC entero, LCA60T [MW]')
val('PilaMaxACCFWMin', _apn['FW']['max_acc']['t_rec_min'], 0, 'con él, la recuperación de 1 t, LCA60T [min]')
# tercera barra (autor, 2026-10-12): HDM a 1 atm, sin el compresor del escape, dentro de la misma carga útil
val('PilaMaxFWUnaAtm', _apn['FW']['max_una_atm']['P_MW'], 1, 'a 1 atm, sin compresor del escape, el rechazo más grande '
    'que cabe en la carga útil, LCA60T [MW]')
# la figura del paper (2026-10-12): el Zeppelin NT suelta 0.2 t
_apf = {nv['corto']: nv for nv in rutas.leer_datos('fig_presupuesto_masa')['apilado_fig']['naves']}
rango('PilaSinGrandes', min(sum(_apn[c]['sin']['t']) for c in ('Hind', 'FW', 'ZR')),
      max(sum(_apn[c]['sin']['t']) for c in ('Hind', 'FW', 'ZR')), 0, 'sin freno, 1 t: LZ 129, LCA60T y ZR-3 [t]')
val('PilaNTSuelta', _apf['ZNT']['masa_t'], 1, 'suelta del Zeppelin NT en la figura del desglose [t]')
val('PilaNTSin', sum(_apf['ZNT']['sin']['t']), 0, 'sin freno, esa suelta en el Zeppelin NT: equipo y carga [t]')
val('PilaNTMax', _apf['ZNT']['max']['P_MW'], 2, 'con freno, el rechazo más grande que cabe para esa suelta, NT [MW]')
val('PilaNTMaxMin', _apf['ZNT']['max']['t_rec_min'], 0, 'con él, la recuperación de esa suelta, NT [min]')
val('PilaNTMaxH', _apf['ZNT']['max']['t_rec_min'] / 60.0, 0, 'lo mismo en horas (2026-10-13: con una unidad entera del '
    'compresor de H2) [h]')
if _apn['ZNT']['max'] is not None or _apn['ZNT']['max_una_atm'] is not None:
    sys.exit('fig_presupuesto_masa: en el Zeppelin NT ahora cabe algo; revisar el texto de la sección 5.2')

# ── EXP-18 (autor, 2026-10-11): alturas de las cuatro naves con la dinámica del condensador (figura de la sección 5.1)
seccion('EXP-18: alturas de las cuatro naves, 1 t, con la dinámica del condensador (exp18_por_nave.py)')
_e18 = rutas.leer_datos('exp18_por_nave')['agua']
_sube = [_e18[c][k]['dh_arriba_m'] for c in _e18 for k in ('100', '10')]
rango('AlturasSube', min(_sube), max(_sube), 1, 'subida máxima con 100 o 10 MW, las cuatro naves [m]')
# carta de diseño (EXP-18, parte 4): m_max / C = 1 / (1 + b + a Q / t), la misma fracción en las cuatro naves
_dd = rutas.leer_datos('exp18_por_nave')['diseno']
_a, _b, _Qd = _dd['a_t_por_MW'], _dd['b_t_por_t'], _dd['Q_MJ_por_t']
val('DisenoA', _a, 1, 'con freno: masa de lo que escala con el rechazo (haces, motor, compresor del escape) [t por MW]')
val('DisenoB', _b, 2, 'cajas fijas por tonelada liberada [t por t]')
val('DisenoQ', _Qd / 1e3, 2, 'calor por tonelada liberada, phi objetivo [GJ por t]')
# (2026-10-13) el compresor de H2 en unidades de catálogo: la fracción ya no es la misma en toda nave; se da para las
# tres grandes (en el NT la unidad pesa más: DisenoTechoNT)
_GR = ('Hind', 'FW', 'ZR')
_pct = lambda c, x: 100.0 * x / _dd['naves'][c]['carga_t']
rango('DisenoTechoPct', min(_pct(c, _dd['naves'][c]['limite_t']) for c in _GR),
      max(_pct(c, _dd['naves'][c]['limite_t']) for c in _GR), 0,
      'techo de la fracción de la carga útil que se puede soltar (t largo), LZ 129, LCA60T y ZR-3 [%]')
rango('DisenoUnaHoraPct', min(_pct(c, _dd['naves'][c]['tabla']['3600']) for c in _GR),
      max(_pct(c, _dd['naves'][c]['tabla']['3600']) for c in _GR), 0,
      'fracción de la carga útil que se suelta en 1 h, las tres grandes [%]')


def _medio_h(c):
    import fig_presupuesto_masa as _fpm
    from scipy.optimize import brentq as _bq
    n_ = _dd['naves'][c]
    return _bq(lambda ts: _fpm.suelta_unica_max(n_['carga_t'], _dd['b_sin_comp_t_por_t'], _a * _Qd / ts,
                                                _dd['compresor']) - 0.5 * n_['limite_t'], 1.0, 1e9) / 3600.0


rango('DisenoMedioH', min(_medio_h(c) for c in _GR), max(_medio_h(c) for c in _GR), 1,
      'tiempo en que se llega a la mitad del techo, las tres grandes [h]')
val('DisenoTechoNT', _dd['naves']['ZNT']['limite_t'], 2, 'lo más que podría soltar el Zeppelin NT, con t largo [t]')
val('DisenoTechoFW', _dd['naves']['FW']['limite_t'], 0, 'lo mismo en el LCA60T [t]')
# la misma carta a 1 atm, sin el compresor del escape (comparación del panel b; la base sigue a 14 atm, autor 2026-10-11)
_u = _dd['una_atm']
val('DisenoAUnaAtm', _u['a_t_por_MW'], 1, 'a 1 atm, sin compresor del escape: lo que escala con el rechazo [t por MW]')
val('DisenoQUnaAtm', _u['Q_MJ_por_t'] / 1e3, 2, 'a 1 atm: calor por tonelada liberada, phi_max [GJ por t]')
val('DisenoUnaAtmPct', 100.0 * (1.0 - _u['razon_t']), 0, 'a 1 atm: cuánto más corta es la recuperación de la misma '
    'carga, igual en toda nave [%]')
val('DisenoUnaAtmFWMin', _dd['naves']['FW']['t_1t_s_una_atm'] / 60.0, 1, 'a 1 atm: recuperación de 1 t en el LCA60T con '
    'el rechazo más grande que cabe [min]')

# ── EXP-20 (autor, 2026-10-12): repartir la carga en sueltas con el freno (sección 5.3 y suplementario) ──
seccion('EXP-20: repartir la carga en sueltas con el freno (exp20_reparto_sueltas.py)')
_r20 = rutas.leer_datos('exp20_reparto_sueltas')
_s20 = lambda c, s: next(f for f in _r20['naves'][c]['sueltas'] if abs(f['s_t'] - s) < 1e-9)
val('RepartoTauH', _r20['aQ_s'] / 3600, 1, 'a Q: escala de tiempo del reparto (y de la carta de diseño) [h]')
val('RepartoBsinComp', _r20['b_sin_comp_t_por_t'], 2, "b': cajas fijas sin el compresor de H2 (tanques y PEM) [t por t]")
# D (2026-10-13): el tanque de agua (MIL-T-53029C), por tonelada liberada; en el reparto, por tonelada entregada
_ta = rutas.leer_datos('balance_masa')['tanque_agua']
mismo(_r20['w_t_por_t'], _ta['kg_por_t'][0] / 1e3, 'tanque de agua: el b_w de EXP-20 = el de balance_masa (clase liviana)')
rango('TanqueAguaKgT', *_ta['kg_por_t'], 1, 'tanque de agua: kg por tonelada liberada (tela de 30-62 oz/yd2)')
rango('TanqueAguaPct', *[100 * x for x in _ta['kg_por_kg_agua']], 2, 'tanque de agua: % de la masa del agua que guarda')
rango('TanqueAguaTela', *[x * p.OZ / p.YD2 for x in p.TANQUE_AGUA['tela_ozyd2']], 1, 'tela de MIL-T-53029C [kg/m2]')
val('TanqueAguaW', 1e3 * _r20['w_t_por_t'], 1, 'b_w: tanque de agua por tonelada (la clase liviana) [kg por t]')
val('TanqueAguaPorT', _ta['agua_por_t_kg'], 0, 'agua retenida que guarda por tonelada liberada (phi objetivo) [kg]')
# B2 (2026-10-13): margen de diseño del 10 % sobre la entrega de referencia (sueltas de 1 t de 10 min, LCA60T)
_mg = _r20['margen']
mismo(_mg['N_t'], next(f for f in _r20['naves']['FW']['sueltas'] if abs(f['s_t'] - 1.0) < 1e-9)['N_t']['600'],
      'B2: la carga sin margen = la de RepartoFWDiez')
val('MargenMasaDiezT', _mg['N_t'] - _mg['N_masa_t'], 1, 'B2: carga que cuesta un 10 % más de masa en TODO el equipo, '
    'LCA60T, sueltas de 1 t de 10 min [t]')
val('MargenCalorDiezT', _mg['N_t'] - _mg['N_calor_t'], 1, 'B2: carga que cuesta un 10 % más de a Q (calor por tonelada o '
    'rechazo de sobra), el mismo caso [t]')
val('RepartoTauHw', _r20['aQ_s'] / (1.0 + _r20['w_t_por_t']) / 3600, 1,
    'a Q / (1 + b_w): la mitad de C\' / (1 + b_w) se entrega en este tiempo [h]')
# (2026-10-13) el compresor de H2 en unidades de catálogo: una por cada tanto de suelta
val('RepartoCompUnidadT', _r20['compresor']['cap_kgph'] / _r20['compresor']['kgph_por_t'], 1,
    'suelta que una unidad del compresor de H2 alcanza a comprimir al ritmo del electrolizador (1 t en 24 h) [t]')
val('RepartoBLineal', _r20['b_lineal_t_por_t'], 2, "b' + el compresor proporcional a su capacidad (familia universal) [t por t]")
_fw1 = _s20('FW', 1.0)
val('RepartoFWDiez', _fw1['N_t']['600'], 0, 'LCA60T, sueltas de 1 t de 10 min: carga de negocio [t]')
val('RepartoFWTreinta', _fw1['N_t']['1800'], 0, 'lo mismo con 30 min [t]')
val('RepartoFWMax', _fw1['N_max_t'], 1, 'lo mismo con recuperación muy larga (solo cajas fijas) [t]')
val('RepartoEquipoDiez', _r20['naves']['FW']['carga_t'] - _fw1['N_t']['600'], 0,
    'equipo para sueltas de 1 t de 10 min (igual en toda nave) [t]')
# (2026-10-13) los límites de las hélices salen de EXP-22 a 14 atm (la fuente de los números de energía; las hélices
# sujetan la suelta menos lo que la quema absorbe sola), ya no de EXP-20
_ls22 = rutas.leer_datos('exp22_energia_potencia')['limites_sueltas_s']
val('RepartoHelUnoPes', _ls22['1']['potencia_pes'] / 60, 0, 'sueltas de 1 t, 14 atm: las hélices alcanzan hasta esta '
    'recuperación, pesimista (EXP-22) [min]')
val('RepartoHelUnoOpt', _ls22['1']['potencia_opt'] / 3600, 1, 'lo mismo, optimista [h]')
val('RepartoHelDosPes', _ls22['2']['potencia_pes'] / 60, 0, 'sueltas de 2 t, pesimista (EXP-22) [min]')
val('RepartoDiezMin', _s20('FW', 10.0)['t_min_s'] / 60, 0, 'LCA60T, sueltas de 10 t: recuperación mínima [min]')
val('RepartoDiezHelOpt', _ls22['10']['potencia_opt'] / 3600, 1, 'sueltas de 10 t: hélices hasta, optimista (EXP-22) [h]')
val('RepartoVeinteMin', _s20('FW', 20.0)['t_min_s'] / 3600, 1, 'LCA60T, sueltas de 20 t: recuperación mínima [h]')
val('RepartoVeinteHelOpt', _ls22['20']['potencia_opt'] / 3600, 1, 'sueltas de 20 t: hélices hasta, optimista (EXP-22) '
    '[h]')
val('RepartoHindTreinta', _s20('Hind', 1.0)['N_t']['1800'], 0, 'LZ 129, sueltas de 1 t de 30 min: carga de negocio [t]')
val('RepartoZRUnaHora', _s20('ZR', 0.25)['N_t']['3600'], 1, 'ZR-3, sueltas de 250 kg de 1 h: carga de negocio [t]')
_nt = _s20('ZNT', 0.025)
val('RepartoNTUnaHora', _nt['N_t']['3600'], 2, 'Zeppelin NT, sueltas de 25 kg de 1 h: carga de negocio [t]')
val('RepartoNTMax', max(f['N_max_t'] for f in _r20['naves']['ZNT']['sueltas']), 2,
    'Zeppelin NT: lo más que reparte (una unidad del compresor de H2) [t]')
val('RepartoNTUnicaMax', _r20['naves']['ZNT']['unica_max_t'], 2, 'Zeppelin NT: la suelta única más grande con el '
    'compresor de H2 en unidades de catálogo (t muy largo) [t]')
val('RepartoNTUnaHoraProp', (_r20['naves']['ZNT']['carga_t'] - 0.025 * (_r20['b_t_por_t'] - _r20['w_t_por_t']
                                                                      + _r20['aQ_s'] / 3600)) / (1.0 + _r20['w_t_por_t']), 2,
    'lo mismo con el compresor proporcional a la suelta, 25 kg de 1 h [t]')
val('RepartoUnivUnoTreinta', 100 * (1 - 0.01 * (_r20['b_lineal_t_por_t'] + _r20['aQ_s'] / 1800)) / (1.0 + _r20['w_t_por_t']), 0,
    'universal: sueltas del 1 % de la carga útil de 30 min, carga de negocio [% de la carga útil]')

# ── EXP-22 (autor, 2026-10-12): energía y potencia por etapa; la entrega de referencia (final de la sección 5.3) ──
seccion('EXP-22: la entrega de referencia y la energía por etapa (exp22_energia_potencia.py)')
_e22 = rutas.leer_datos('exp22_energia_potencia')
_En, _Ep, _Et = _e22['E'], _e22['E']['paso'], _e22['E']['total']
mismo(_En['n'], int(_fw1['N_t']['600'] + 1e-9), 'EXP-22: sueltas de la entrega = la carga de negocio de EXP-20 (entera)')
mismo(_En['equipo_t'], _r20['naves']['FW']['carga_t'] - _fw1['N_t']['600'], 'EXP-22: equipo = el de RepartoEquipoDiez')
_MWh = lambda x: [v / 1e3 for v in x]
val('EntregaPasos', _En['n'], 0, 'entrega de referencia, LCA60T: sueltas de 1 t de 10 min (= t de carga de negocio)')
val('EntregaHoras', _En['n'] * _En['t_paso_s'] / 3600, 0, 'entrega de referencia: duración, sueltas una tras otra [h]')
# (2026-10-13) la misma entrega a 1 atm, sin el compresor del escape (ec. 12 con a y Q de 1 atm; auditoría, 5.6)
_u20, _w20 = _r20['una_atm'], _r20['w_t_por_t']
_Cp1 = _r20['naves']['FW']['carga_t'] - _fw1['comp_t'] - _r20['b_sin_comp_t_por_t'] * 1.0      # C' de sueltas de 1 t
_T6 = _En['n'] * _En['t_paso_s']
val('RepartoTauHwUnaAtm', _u20['aQ_s'] / (1.0 + _w20) / 3600, 1, 'a Q / (1 + b_w) a 1 atm [h]')
val('RepartoFWSeisUnaAtm', _Cp1 * _T6 / ((1.0 + _w20) * _T6 + _u20['aQ_s']), 0, 'LCA60T, sueltas de 1 t, a 1 atm: '
    'carga de negocio en la duración de la entrega de referencia (ec. 12) [t]')
val('RepartoEquipoDiezUnaAtm', _r20['naves']['FW']['carga_t'] - (_Cp1 - _u20['aQ_s'] / 600.0) / (1.0 + _w20), 0,
    'equipo para sueltas de 1 t de 10 min, 1 atm [t]')
# (2026-10-14, autor) la cota pesimista de la entrega: el sistema ACC entero en vez de sus haces (sección 5.3, conclusión 5)
_acc20 = _r20['acc_entero']
mismo(_acc20['naves']['FW']['t_min_s'] / 60.0, _apn['FW']['max_acc']['t_rec_min'],
      'ACC entero: la suelta de 1 t más rápida del reparto = la recuperación de PilaMaxACCFWMin (otro camino)')
val('RepartoTauHwACC', _acc20['aQ_s'] / (1.0 + _w20) / 3600, 1, 'a Q / (1 + b_w) con el sistema ACC entero [h]')
val('RepartoFWSeisACC', _Cp1 * _T6 / ((1.0 + _w20) * _T6 + _acc20['aQ_s']), 0, 'LCA60T, sueltas de 1 t, ACC entero: '
    'carga de negocio en la duración de la entrega de referencia (ec. 12) [t]')
val('RepartoFWTreintaACC', _acc20['naves']['FW']['N_t']['1800'], 0, 'LCA60T, sueltas de 1 t de 30 min, ACC entero: '
    'carga de negocio [t]')
val('RepartoHindUnaHoraACC', _acc20['naves']['Hind']['N_t']['3600'], 1, 'LZ 129, sueltas de 1 t de 1 h, ACC entero: '
    'carga de negocio [t]')
val('EntregaPrechMW', _Ep['P_rech_W'] / 1e6, 1, 'entrega de referencia: rechazo instalado [MW]')
val('EntregaHdosPaso', _Ep['H2_kg'], 0, 'entrega de referencia: H2 quemado por suelta [kg]')
val('EntregaHdosTotal', _Et['H2_kg'] / 1e3, 2, 'entrega de referencia: H2 quemado en toda la entrega [t]')
val('EntregaAguaPaso', _Ep['agua_kg'], 0, 'entrega de referencia: agua retenida por suelta [kg]')
val('EntregaAguaTotal', _Et['agua_kg'] / 1e3, 1, 'entrega de referencia: agua retenida en toda la entrega [t]')
rango('EntregaHelPasoMW', _Ep['P_hel_W'][0] / 1e6, _Ep['P_hel_W'][1] / 1e6, 2,
      'entrega de referencia: hélices al soltar, discos mayores-menores [MW]')
rango('EntregaHelPaso', *_Ep['hel_kWh'], 0, 'entrega de referencia: energía de las hélices por suelta [kWh]')
rango('EntregaHelTotal', *_MWh(_Et['hel_kWh']), 1, 'entrega de referencia: energía de las hélices, toda la entrega [MWh]')
rango('EntregaNetoPaso', *_Ep['neto_kWh'], 0, 'entrega de referencia: neto a la planta por suelta, ya descontados '
      'ventiladores, compresor del escape y hélices [kWh]')
rango('EntregaNetoTotal', *_MWh(_Et['neto_kWh']), 1, 'entrega de referencia: neto a la planta, toda la entrega [MWh]')
rango('EntregaRegenPaso', *_MWh(_Ep['regen_kWh']), 2, 'entrega de referencia: regeneración (electrólisis + compresión) '
      'por suelta [MWh]')
rango('EntregaRegenTotal', *_MWh(_Et['regen_kWh']), 0, 'entrega de referencia: regeneración, toda la entrega [MWh]')
rango('EntregaExtPaso', *_MWh(_Ep['externa_kWh']), 2, 'entrega de referencia: electricidad de tierra neta por suelta [MWh]')
rango('EntregaExtTotal', *_MWh(_Et['externa_kWh']), 0, 'entrega de referencia: electricidad de tierra neta, toda la '
      'entrega [MWh]')
val('EntregaReposKg', _Et['H2_repos_kg'], 0, 'entrega de referencia: H2 de reposición desde tierra (agua no retenida) [kg]')
val('EntregaPemKW', _En['P_pem_kW'], 0, 'PEM del equipo: regenera una suelta de 1 t en t_regeneracion_h [kW]')
val('EntregaRegenDias', _En['t_regen_pem_h'] / 24, 0, 'entrega de referencia: días de regeneración con ese PEM')
rango('EntregaBaseMW', *_En['P_base_MW'], 1, 'entrega de referencia: electrólisis en tierra para regenerar toda la '
      'entrega en t_regeneracion_h [MW]')
val('EntregaGasIni', _En['gas_kg'][0] / 1e3, 1, 'entrega de referencia: H2 de los gases al salir [t]')
val('EntregaGasPct', 100 * (_En['gas_kg'][0] - _En['gas_kg'][-1]) / _En['gas_kg'][0], 0,
    'entrega de referencia: H2 de los gases quemado en la entrega [%]')
for _c, _nv in _e22['F']['naves'].items():
    val('SinFrenoCabe' + _c, 1000 * _nv['cruce_t'], 0, f"{_nv['nombre']}: la suelta más grande que se compensa sin freno "
        'con un equipo que cabe en la carga útil [kg]')
_lim = _e22['limites_s']
val('EtapasLimPotPes', _lim['potencia_pes'] / 60, 0, 'LCA60T 1 t: hélices al soltar <= potencia neta del motor hasta, '
    'pesimista (fuerza tras lo que absorbe la quema) [min]')
val('EtapasLimPotOpt', _lim['potencia_opt'] / 3600, 1, 'lo mismo, optimista [h]')
val('EtapasLimEnePes', _lim['energia_pes'] / 3600, 1, 'LCA60T 1 t: energía de las hélices <= neto del motor hasta, '
    'pesimista [h]')
val('EtapasLimEneOpt', _lim['energia_opt'] / 3600, 1, 'lo mismo, optimista [h]')

# ── EXP-21/22 para el paper (autor, 2026-10-13): sección 6 "Energy and regeneration" y la ley de entrega de la 5.3 ──
seccion('EXP-21/22 en el paper: energía de la suelta, regeneración y ley de entrega')
_f21 = rutas.leer_datos('exp21_balance_energia')['filas']
_k14, _k1 = _f21[1]['kWh'], _f21[2]['kWh']                      # HD 14 atm 10 MW; HD 1 atm 10 MW
# (2026-10-13, objeción 3.18) regenerar en vuelo ventea el O2 del agua y aligera la nave; compensarlo comprimiendo gas
# de las celdas en los tanques pide m_O2 / L_comp kg de H2, y a los tanques (dimensionados para todo el H2 quemado por
# tonelada) solo les queda lo que no ocupa el regenerado
_tqT = next(c for c in rutas.leer_datos('balance_masa')['cajas'] if c['caja'] == 'tanques_por_t')
_tqKg = rutas.leer_datos('balance_masa')['tanque_por_kg_H2']
_Lc = rutas.leer_datos('exp12_helio')['h2_comprimido']['L']
for _suf, _f in (('Cat', _f21[1]), ('Una', _f21[2])):
    _mr = _f['masa_kg']['H2_regenerado']
    _o2 = _mr * p.ratio_O2_H2
    # el regenerado y el gas que compensa su O2 (r / L_comp por kg) van a los mismos tanques: se llenan cuando se
    # regeneró capacidad / (1 + r / L_comp)
    _xm = [c / _tqKg / (1.0 + p.ratio_O2_H2 / _Lc) for c in (_tqT['min_kg'], _tqT['max_kg'])]   # [kg por t]
    _p = 'a 14 atm' if _suf == 'Cat' else 'a 1 atm'
    val('RegenOdosVent' + _suf, _o2, 0, f'O2 venteado al regenerar el H2 de una tonelada {_p} (sustentación que gana) [kg]')
    val('RegenHdosComp' + _suf, _o2 / _Lc, 1, f'H2 a comprimir en los tanques para compensarlo todo {_p} [kg por t]')
    rango('RegenVueloPct' + _suf, 100 * min(1.0, _xm[0] / _mr), 100 * min(1.0, _xm[1] / _mr), 0,
          f'parte del H2 de una tonelada que se puede regenerar en vuelo compensando solo con los tanques {_p} [%]')
    # tanques que guardan también el gas que compensa: el piso de la ec. 13 sube en (1 + r / L_comp)
    val('RegenPisoComp' + _suf, _tqKg * _mr * (1.0 + p.ratio_O2_H2 / _Lc) / 1e3, 2,
        f'piso del costo por tonelada recuperada con tanques también para el gas que compensa el O2 {_p} [t por t]')
mismo(_k1['electrolisis'], rutas.leer_datos('balance_masa')['E_el_por_t_kWh'],
      'EXP-21: electrólisis por tonelada a 1 atm = balance_masa (ElecPorT)', tol=1e-3)
for _suf, _k in (('Cat', _k14), ('Una', _k1)):
    _reg = [_k['electrolisis'] + min(_k['compresion_H2']), _k['electrolisis'] + max(_k['compresion_H2'])]
    _ext = [_reg[0] - _k['neto_planta'][1], _reg[1] - _k['neto_planta'][0]]
    if _suf == 'Una':
        _ext_una = list(_ext)
    else:
        _ext_cat = list(_ext)
    _p = 'a 14 atm' if _suf == 'Cat' else 'a 1 atm'
    val('EnergiaQuim' + _suf, _k['quimica_PCI'], 0, f'energía química quemada por tonelada {_p} [kWh]')
    rango('EnergiaMotor' + _suf, *_k['trabajo_motor'], 0, f'trabajo del motor por tonelada {_p} [kWh]')
    rango('EnergiaNeto' + _suf, *_k['neto_planta'], 0, f'neto a la planta por tonelada {_p}, HD 10 MW, ya sin '
          'ventiladores, compresor del escape (14 atm) y hélices [kWh]')
    val('EnergiaElec' + _suf, _k['electrolisis'], 0, f'electrólisis por tonelada {_p} [kWh]')
    rango('EnergiaComp' + _suf, *_k['compresion_H2'], 0, f'compresión del H2 regenerado por tonelada {_p} [kWh]')
    rango('EnergiaRegen' + _suf, *_reg, 0, f'regeneración por tonelada {_p} (electrólisis + compresión) [kWh]')
    rango('EnergiaPaga' + _suf, 100 * _k['neto_planta'][0] / _reg[1], 100 * _k['neto_planta'][1] / _reg[0], 0,
          f'parte de la regeneración que paga el neto del motor {_p} [%]')
    rango('EnergiaExt' + _suf, *_ext, 0, f'electricidad de tierra neta por tonelada {_p} [kWh]')
val('EnergiaCompEsc', _k14['compresor_escape'], 0, 'compresor del escape por tonelada a 14 atm [kWh]')
# M7 frente a ventear, con la electricidad de tierra de EXP-22 a 1 atm
rango('CicloExtMWh', _ext_una[0] / 1e3, _ext_una[1] / 1e3, 1,
      'electricidad externa neta del ciclo por tonelada, 1 atm (EXP-22) [MWh]')
rango('VentAhorroPct', 100 * (1 - _ext_una[1] / 1e3 / _Event), 100 * (1 - _ext_una[0] / 1e3 / _Event), 0,
      'energía que el ciclo ahorra frente a reponer lo venteado, 1 atm (EXP-22) [%]')
rango('VentAhorroPctCat', 100 * (1 - _ext_cat[1] / 1e3 / _Event), 100 * (1 - _ext_cat[0] / 1e3 / _Event), 0,
      'lo mismo a 14 atm (EXP-22) [%]')
# (2026-10-13) frente a ventear, con los datos de la tabla 6 (auditoría, objeción 1.1: el 90 % salía de EXP-11, escalones
# sin freno a 1 atm con retardo de 5 s, que el paper ya no muestra)
val('VentEvitaUnaAtm', 100 * (1 - _c1['perdido_kg_por_t'] / _vent), 0, 'H2 que el ciclo deja de perder frente a '
    'ventear, 1 atm (tabla 6) [%]')
val('VentEvitaCatorce', 100 * (1 - _c99['perdido_kg_por_t'] / _vent), 0, 'lo mismo a 14 atm, phi = 0.99 (tabla 6) [%]')
val('DoceInstantePctUnaAtm', _c1['al_instante_kg'] / 10.0, 0, 'ciclo a 1 atm: sustentación devuelta al instante '
    '[% de 1 t]')
_heC = next(x for x in r12['rutas'] if 'He comprimido' in x['ruta'])
mismo(_c99['calor_GJ_por_t'] / _heC['calor_GJ_por_t'], r12['factor_calor_20MPa'], 'factor de calor H2/He: la ruta de '
      'He comprimido = la de factor_calor_20MPa')
val('DoceFactorCalorUnaAtm', _c1['calor_GJ_por_t'] / _heC['calor_GJ_por_t'], 0, 'calor por tonelada a 1 atm / He a '
    '20 MPa')
_l1 = _e22['limites_1atm_s']
val('EnergiaHelPotUnaPes', _l1['potencia_pes'] / 3600, 1, '1 atm, 1 t: hélices al soltar <= potencia neta del motor '
    'hasta, pesimista [h]')
val('EnergiaHelPotUnaOpt', _l1['potencia_opt'] / 3600, 1, 'lo mismo, optimista [h]')
val('EnergiaHelEneUnaPes', _l1['energia_pes'] / 3600, 1, '1 atm, 1 t: energía de las hélices <= neto del motor hasta, '
    'pesimista [h]')
val('EnergiaHelEneUnaOpt', _l1['energia_opt'] / 3600, 0, 'lo mismo, optimista [h]')
# el módulo de regeneración (por kg de H2 por hora, a tope) y el costo por tonelada recuperada
_I = _e22['I']
rango('RegenKWh', *_I['kW_por_kgph'], 1, 'módulo de regeneración: kWh por kg de H2 (electrólisis + compresión)')
rango('RegenPemKg', *_I['pem_kg_por_kgph'], 0, 'PEM por kg de H2 por hora [kg]')
val('RegenCompCap', _I['cap_comp_kgph'], 2, 'una unidad de compresor de catálogo [kg de H2 por hora]')
val('RegenCompKg', _I['comp_kg_unidad'], 0, 'masa de esa unidad [kg]')
val('RegenTanqueKg', _I['tanques']['350']['kg_por_kg'], 1, 'tanques a 350 bar por kg de H2 guardado [kg]')
val('RegenHdosT', _f21[1]['masa_kg']['H2_regenerado'], 1, 'H2 regenerado por tonelada liberada a 14 atm (el del '
    'nomograma y el costo por tonelada) [kg]')
val('RegenTanqueSiete', _I['tanques']['700']['kg_por_kg'], 1, 'tanques a 700 bar por kg de H2 guardado [kg]')
_G = _e22['G']['G_kWh_m2_dia']
val('RegenSolAnual', _G['anual'], 1, 'horas de sol pleno por día sobre los paneles del casco, Lima, promedio anual')
val('RegenSolInv', _G['invierno'], 1, 'lo mismo, invierno')
val('RegenSolVer', _G['verano'], 1, 'lo mismo, verano')
_C4 = _e22['I4']
val('CostoPiso', _C4['piso_t_por_t'], 2, 'carga útil por tonelada recuperada a bordo: piso de los tanques [t por t]')
for _f, _suf in (('a-Si', 'ASi'), ('CIGS', 'Cigs'), ('red', 'Red')):
    _x = _C4['fuentes'][_f]
    rango('CostoSemana' + _suf, *_x['costo']['7'], 2, f'carga útil por tonelada recuperada en una semana, {_f} [t por t]')
    rango('CostoEquil' + _suf, *_x['T_equilibrio_d'], 0 if _x['T_equilibrio_d'][0] > 5 else 1,
          f'cuesta menos de 1 t por t recuperada desde estos días, {_f}')
# la ley de entrega de la 5.3: el techo de la suelta única con el compresor de catálogo fijo
val('RepartoUnicaMaxFW', _r20['naves']['FW']['unica_max_t'], 0, 'LCA60T: la suelta única más grande con el compresor '
    'de H2 fijo (recuperación muy larga) [t]')

# ── Escribir numeros.tex ─────────────────────────────────────────────────────
lineas = [
    '% ' + '=' * 77,
    '% numeros.tex — ARCHIVO GENERADO por modelo/exportar_numeros.py. NO EDITAR A MANO.',
    '% Regenerar: python correr_todo.py   (desde la carpeta modelo/)',
    '% Uso: \\input{.../modelo/resultados/numeros.tex} en el preámbulo; en el texto',
    '%      se escribe \\valAlphaStar en vez de 0.497 (en modo texto: \\valAlphaStar{}).',
    '% ' + '=' * 77,
]
n = 0
for nombre, texto, comentario in MACROS:
    if nombre is None:
        lineas += ['', f'% --- {comentario} ---']
    else:
        lineas.append(f'\\newcommand{{\\{nombre}}}{{{texto}}}  % {comentario}')
        n += 1
with open(rutas.os.path.join(rutas.RESULTADOS, 'numeros.tex'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lineas) + '\n')
print(f'numeros.tex: {n} macros')
