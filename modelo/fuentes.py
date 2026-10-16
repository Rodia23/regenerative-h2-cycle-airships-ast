# -*- coding: utf-8 -*-
# [EN] Register and audit of the published sources used by the model (Supplementary S4.2).
# =============================================================================
# fuentes.py — REGISTRO Y AUDITORÍA DE LAS FUENTES DE VALIDACIÓN (VAL-00).
#
# El pipeline verifica NUESTRO código; este registro verifica las FUENTES con las que
# lo validamos. Una fuente con un error que no se detecta entra "validada" y contamina
# todo lo que se compara con ella (caso real: Zakhvatkin 2022, tabla 1, ver ZAK22).
#
# PROTOCOLO (cada fuente pasa los 6 pasos ANTES de usarse; queda escrito aquí):
#   1 Identidad   cita, DOI, año, si es arbitrada (revisada por pares) o preprint,
#                 archivo local + SHA-256 (sabemos que usamos EXACTAMENTE ese PDF),
#                 formato (texto nativo o escaneo con OCR), erratas/retractaciones
#                 (Crossref, con fecha de consulta).
#   2 Datos       cada dato con su ubicación (tabla, figura, página, ecuación) y el
#                 método (transcrito; leído de figura con su incertidumbre). Si el
#                 texto extraído del PDF sale desalineado o es OCR, el dato se
#                 confirma contra la IMAGEN de la página.
#   3 Consistencia interna: recalcular sus números con sus propias ecuaciones y
#                 unidades; tabla vs apéndice vs figura vs resumen/conclusiones.
#   4 Contraste   para lo crítico, una segunda fuente o un límite físico.
#   5 Hallazgos   todo lo que no cuadra, escrito (aunque no cambie la decisión).
#   6 Decisión    se usa / se usa con corrección / solo cualitativo / se descarta.
#
# PEDIGRÍ (calificación, al estilo del factor "input pedigree" de NASA-STD-7009A; no se
# caza cada error de cada paper, se califica cuánto se puede confiar en lo que se usa):
#   A  dato de referencia o medición primaria, revisada, reproducida por nosotros
#      (tablas NIST/JANAF/IAPWS; un punto publicado que el modelo reproduce sin ajustar)
#   B  medición o modelo validado publicado y revisado, con alguna limitación
#      declarada (escala, condiciones, discrepancias menores de redacción)
#   C  cualitativo, preliminar, o sin los datos que respaldan sus conclusiones
# Un experimento a escala reducida que valida un MODELO (caso LXS14) no es un defecto:
# es la práctica normal; lo que se transfiere a escala real es el modelo validado, no
# los números del experimento.
#
# La verificación (grupo P de verificacion.py) comprueba que cada PDF existe y no
# cambió (SHA-256), y que los datos que usa el modelo siguen siendo IDÉNTICOS a los
# registrados (las transcripciones largas, por su huella). `python fuentes.py`
# escribe el reporte legible resultados/fuentes_registro.txt.
# Los PDF de terceros no se distribuyen (derechos de autor): si faltan, su identidad queda "sin comprobar".
# =============================================================================

import hashlib
import os
import sys

RAIZ = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
REF = 'referencias_fuentes/'

FUENTES = {
    # ── Ya usadas por el modelo ─────────────────────────────────────────────
    'YG03': {
        'ambito': 'paper 1',
        'cita': 'Yamasaki, T., Goto, N., Identification of blimp dynamics via flight tests, Trans. Japan Soc. '
                'Aero. Space Sci. 46 (153) (2003) 195-205',
        'doi': '10.2322/tjsass.46.195', 'arbitrada': True, 'acceso': 'abierto (J-STAGE)',
        'archivo': REF + 'Yamasaki_Goto_2003_46_195.pdf',
        'sha256': '7d5d67eb4bea3369b05b68cfd6e0b95a930777abb646d13a352d52f660b17689',
        'formato': 'texto nativo', 'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['VAL-01 (masa añadida)', 'paper 1, discusión'],
        'datos': {  # tabla 4, pág. 201; CONFIRMADA contra la imagen
            'analitica': (38.44, 'tabla 4, MLG(2,2) Ana. [kg]'),
            'restringido': (29.65, 'tabla 4, MLG(2,2) Con. [kg]'),
            'vuelo_libre': (46.19, 'tabla 4, MLG(2,2) Fr. [kg]'),
        },
        'modelo': ('val01_componentes', 'YG03_MASA_VERTICAL'),
        'datos_no_usados': {'ALG(2,2) amortiguamiento vertical [N s/m]': (-2.52, -6.66, -0.66),
                            'dirigible': 'Sky Probe-J, 18.0 m3, 5.48 m x 2.50 m (tabla 1)'},
        'consistencia': ['los dos ensayos (restringido y vuelo libre) difieren 1.6x en la masa vertical y 10x '
                         'en el amortiguamiento: se usa el rango, no un valor'],
        'hallazgos': ['el texto extraído corre la columna "Ana." de la tabla 4; confirmado contra la imagen',
                      'el amortiguamiento es de vuelo hacia adelante a 0.5 m/s (sustentación del casco), no de '
                      'un dirigible quieto: no valida el arrastre vertical del modelo'],
        'decision': 'se usa (masa añadida, como rango; propagado en VAL-03)', 'pedigri': 'B',
    },
    'ZAK22': {
        'ambito': 'paper 1',
        'cita': 'Zakhvatkin, L., Schechter, A., Avrahami, I., Water recuperation from hydrogen fuel cell during '
                'aerial mission, Energies 15 (2022) 6848',
        'doi': '10.3390/en15186848', 'arbitrada': True, 'acceso': 'abierto (MDPI)',
        'archivo': REF + 'energies-15-06848-v2.pdf',
        'sha256': '573aec2af35771239c9ae6aed3180c46718106ce04adcaea6ef4d1178a639172',
        'formato': 'texto nativo', 'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['VAL-02 (Ec. 5, condensación)', 'paper 1, validación'],
        # su propio cociente medido / esperado (contra SU modelo), para el texto de la validación (hallazgo M2)
        'datos': {'recoleccion_sobre_su_modelo': ((0.70, 0.85), 'p. 10 ("which is 70-85% of the theoretical amounts '
                                                                 'of available water") y p. 13 (conclusiones)')},
        'transcripciones': {'tabla 1 (6 experimentos)': ('val02_condensacion', 'EXPERIMENTOS',
                                                         '9d7529f332838c2ff9b8c2878b870087de1dc7fa58f2602ae64c7a1143dc8bf1')},
        'consistencia': ['su figura 4 reproducida con sus ecuaciones (0.12 puntos; Tetens vs IF97). Su modelo (p. 6) es '
                         'el mismo balance de Dalton que la ec. 5 con la presión de saturación de Tetens: la '
                         'comparación con sus curvas es de CONSISTENCIA, no validación (hallazgo M2); lo medido es '
                         'el agua recogida',
                         'agua esperada recalculada con sus flujos medidos: 1.5-2.6 % por debajo de su apéndice A'],
        'hallazgos': ['experimento 5: el agua esperada de la tabla no coincide con la de su apéndice',
                      'experimento 6: la humedad de saturación de la tabla no corresponde a su temperatura',
                      'los lambda publicados no salen de sus propios flujos: no se usan; se usan los flujos'],
        'decision': 'se usa con corrección (flujos medidos, no sus lambda; eficiencia de recolección 0.68-1.0 '
                    'propagada en VAL-03)', 'pedigri': 'B',
    },
    # ── Nuevas (2026-09-25): candidatas para validar dinámica, gas isotérmico y tau ──
    'IZ25': {
        'ambito': 'paper 1',
        'cita': 'Izraelevitz, J. S., et al., Flight demonstration and model validation of a prototype '
                'variable-altitude Venus aerobot, J. Aircraft (2025)',
        'doi': '10.2514/1.C038314', 'arbitrada': True,
        'acceso': 'PDF local = preprint arXiv 2411.06643v2 (8 oct 2025); la versión de la revista no se tiene',
        'archivo': REF + 'Izraelevitz_2024_arXiv_2411.06643v2.pdf',
        'sha256': '0a98c58e61c967194ef356e98c257bc0c521b39076f12c203517b1ced5773f63',
        'formato': 'texto nativo', 'erratas': 'Crossref: el DOI de arXiv no está en Crossref; revisar el de la revista',
        'usos': ['candidata: contraste de parámetros de la dinámica vertical'],
        'datos': {'Cd_vertical': (0.8, 'sec. III.B.3 (C_D,top, ajustado a vuelos en hangar, ref. 13)'),
                  'Cd_lateral': (1.0, 'sec. III.B.3'),
                  'coef_masa_virtual': (0.2, 'sec. III.B.3 (datos de elipsoides de Hoerner)'),
                  'error_termico_C': (7.0, 'sec. V (modelo vs medido, +-7 C)')},
        'consistencia': ['pendiente: confirmar los números con la versión de J. Aircraft'],
        'hallazgos': ['preprint posterior a la publicación; los coeficientes son AJUSTADOS a sus vuelos, no '
                      'medidos directamente; forma de globo (no casco esbelto)'],
        'decision': 'contraste de parámetros (no validación directa): su Cd 0.8 entra en la sensibilidad de '
                    'VAL-03', 'pedigri': 'B',
    },
    'RM93': {
        'ambito': 'paper 1',
        'cita': 'Robbins, E., Martone, M., Determination of balloon gas mass and revised estimates of drag and '
                'virtual mass coefficients, Adv. Space Res. 13 (2) (1993) 55-57',
        'doi': '10.1016/0273-1177(93)90275-G', 'arbitrada': True, 'acceso': 'pago (descargado por el autor)',
        'archivo': REF + 'Robbins_Martone_1993_ASR_13_55.pdf',
        'sha256': 'fdac577c9443633513c61964c89ac7a1a8c4afadfc2163fb94321e67cabf940c',
        'formato': 'ESCANEO con OCR (TextBridge): números solo contra la imagen',
        'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['respaldo cualitativo del Cd'],
        'datos': {'Cd_globo': (0.5, 'pág. (2)57, Results: "consistently around 0.5", Re 3.5e5-9e5')},
        'consistencia': ['el texto cita las figuras 1-3 (aceleración, velocidad y altitud medidas) y NO aparecen '
                         'en las 3 páginas publicadas (imagen revisada)'],
        'hallazgos': ['resultados declarados preliminares por los autores; la masa virtual solo "mayor que 0.5", '
                      'sin número; no hay datos del transitorio al soltar lastre'],
        'decision': 'solo cualitativo: NO sirve para la validación de punta a punta', 'pedigri': 'C',
    },
    'LXS14': {
        'ambito': 'paper 1',
        'cita': 'Li, D.-F., Xia, X.-L., Sun, C., Experimental investigation of transient thermal behavior of an '
                'airship under different solar radiation and airflow conditions, Adv. Space Res. 53 (2014) 862-869',
        'doi': '10.1016/j.asr.2013.12.032', 'arbitrada': True, 'acceso': 'pago (descargado por el autor)',
        'archivo': REF + 'Li_Xia_Sun_2014_ASR_53_862.pdf',
        'sha256': '3f4c68b6d63510d6c26fb744d6689173f94de109d44e7072feb343f503f9a094',
        'formato': 'texto nativo', 'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['candidata: orden de magnitud del calentamiento solar del gas (supuesto isotérmico)'],
        'datos': {'modelo_D_m': (0.47, 'sec. 2'), 'modelo_L_m': (1.42, 'sec. 2'),
                  'gas_cerca_sol_C': ((13.7, 38.6), 'sec. 3.1, fig. 3c (T18, 972 W/m2, conv. natural)'),
                  'gas_cerca_sombra_dT_C': (6.7, 'sec. 3.1'),
                  'casco_estable_s': ((600.0, 800.0), 'sec. 3.1 (972 y 842.4 W/m2)')},
        'consistencia': ['conclusiones vs cuerpo: sol/sombra 76.4 C en resultados vs "40-70 C" en conclusiones',
                         'casco-gas 72.4 C en resultados vs "cerca de 60 C" en conclusiones',
                         'casco estable a los 800 s (13 min) con 842.4 W/m2 vs "5-10 min" en resumen y conclusiones',
                         'velocidad del aire: 5.5 m/s en el resumen, 5.3 m/s en el pie de la figura 6'],
        'hallazgos': ['experimento a escala reducida (1.42 m, 0.25 m3, simulador solar) para VALIDAR un modelo '
                      'térmico (sus figs. 9-10: 0.2-1.1 C): práctica normal, no un defecto; sus tiempos no son de '
                      'escala real (el tiempo térmico crece con el tamaño)',
                      'las mediciones del cuerpo son coherentes; lo flojo son el resumen y las conclusiones (cuatro '
                      'discrepancias arriba, confirmadas contra las páginas 1, 4, 5 y 8): se usan los del cuerpo'],
        'decision': 'respaldo del tipo de modelo térmico (convección + radiación) que VAL-03 aplica a escala real '
                    'con correlaciones estándar; no se escalan sus temperaturas', 'pedigri': 'B',
    },
    'GH23': {
        'ambito': 'paper 1',
        'cita': 'Ghatos, S., Venegas, M., de Vega, M., Taha-Janan, M., García-Hernando, N., Experimental evaluation '
                'of the transient and steady state performance of an air-cooled square minichannel heat exchanger, '
                'Int. J. Therm. Sci. 188 (2023) 108255',
        'doi': '10.1016/j.ijthermalsci.2023.108255', 'arbitrada': True, 'acceso': 'pago (descargado por el autor)',
        'archivo': REF + 'Ghatos_2023_IJTS_188_108255.pdf',
        'sha256': '9e043f7aee3eadb0dbc2c4fc7c5a1bdd41b9ccd9db879daafdfe2973f4d58e0f',
        'formato': 'texto nativo', 'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['candidata: respaldo del retardo tau del condensador'],
        'datos': {'tau_forzada_min': ((0.26, 0.99), 'tabla 8 (aire a 0.75 m/s)'),
                  'tau_natural_min': ((1.16, 1.78), 'tabla 8'),
                  'transitorio_forzada_min': ((5.0, 16.0), 'sec. 4.2'),
                  'capacidad_J_K': (217.0, 'tabla 8 (216.9-217.6)')},
        'consistencia': ['tau = C/(UA): con C = 217 J/K y tau 16-60 s, UA = 3.6-12 W/K, plausible para aire a '
                         '0.75 m/s (pendiente: contrastar con sus resistencias de la fig. 5)'],
        'hallazgos': ['intercambiador chico y con aire lento: su tau es cota alta para un radiador compacto '
                      'aeronáutico', 'la tabla 8 tiene celdas vacías (casos solo con convección forzada)'],
        'decision': 'cota superior de tau (orden de magnitud), usada en VAL-03', 'pedigri': 'B',
    },
    # ── Referencias de las constantes del modelo (VAL-01); sin PDF local: son normas y bases
    #    de datos públicas. 'modelo': (módulo, None) = los datos son atributos del módulo.
    'ISA': {
        'ambito': 'paper 1',
        'cita': 'ISO 2533:1975 Standard Atmosphere; U.S. Standard Atmosphere 1976 (NOAA-S/T 76-1562), tabla 3 '
                '(composición del aire seco)',
        'doi': 'no tiene (norma ISO y reporte técnico NOAA/NASA/USAF)', 'arbitrada': True,
        'acceso': 'público (NTRS 19770009539)', 'archivo': None, 'sha256': None, 'formato': 'norma',
        'erratas': 'no aplica (norma vigente)', 'usos': ['VAL-01', 'parametros (atmósfera, aire)'],
        'datos': {'g': (9.80665, 'gravedad estándar'), 'T_ref_atm': (288.15, 'T a nivel del mar'),
                  'P_atmos': (101325, 'P a nivel del mar'), 'R_aire': (287.05287, 'R del aire seco'),
                  'gradiente_ISA': (0.0065, 'gradiente de la troposfera'), 'rho_aire_std': (1.225, 'rho a nivel del mar'),
                  'X_O2_aire': (0.209476, 'USSA 1976, tabla 3'), 'M_aire': (28.9644e-3, 'USSA 1976')},
        'modelo': ('parametros', None),
        'consistencia': ['rho = P/(R T) = 1.2250 kg/m3 (grupo H)'], 'hallazgos': [],
        'decision': 'se usa', 'pedigri': 'A',
    },
    'CODATA': {
        'ambito': 'paper 1',
        'cita': 'Tiesinga et al., CODATA recommended values 2018 (R); Cox, Wagman, Medvedev, CODATA Key Values for '
                'Thermodynamics, Hemisphere 1989 (entalpía de formación del H2O gas)',
        'doi': '10.1103/RevModPhys.93.025010', 'arbitrada': True, 'acceso': 'público', 'archivo': None,
        'sha256': None, 'formato': 'valores recomendados', 'erratas': 'no aplica (valores recomendados)',
        'usos': ['VAL-01', 'parametros (R_univ, PCI del H2)'],
        'datos': {'R_univ': (8.314462618, 'CODATA 2018, exacta')},
        'modelo': ('parametros', None),
        'consistencia': ['PCI = 241.826 kJ/mol / M_H2 = entalpía de formación de JANAF (grupo K)',
                         'PCS - PCI = (M_H2O/M_H2) h_fg a 25 C de NIST (grupo H)'],
        'hallazgos': [], 'decision': 'se usa', 'pedigri': 'A',
    },
    'IUPAC': {
        'ambito': 'paper 1',
        'cita': 'Wieser, M. E., Atomic weights of the elements 2005 (IUPAC Technical Report), Pure Appl. Chem. 78 '
                '(2006) 2051-2066: H 1.00794, O 15.9994',
        'doi': '10.1351/pac200678112051', 'arbitrada': True, 'acceso': 'público', 'archivo': None, 'sha256': None,
        'formato': 'reporte técnico IUPAC', 'erratas': 'reemplazado por intervalos en 2009+; se usan los valores '
                                                        'convencionales de 2005, dentro de esos intervalos',
        'usos': ['VAL-01', 'parametros (masas molares)'],
        'datos': {'M_H2': (2.01588e-3, '2 x 1.00794'), 'M_O2': (31.9988e-3, '2 x 15.9994'),
                  'M_H2O': (18.01528e-3, '2 x 1.00794 + 15.9994')},
        'modelo': ('parametros', None),
        'consistencia': ['M_H2O = M_H2 + M_O2/2 exacto; ratio O2/H2 = 7.937'], 'hallazgos': [],
        'decision': 'se usa', 'pedigri': 'A',
    },
    'IF97': {
        'ambito': 'paper 1',
        'cita': 'IAPWS R7-97(2012), Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic '
                'Properties of Water and Steam, ec. 30 (presión de saturación) y tabla 35 (verificación)',
        'doi': 'no tiene (release de IAPWS, iapws.org)', 'arbitrada': True, 'acceso': 'público', 'archivo': None,
        'sha256': None, 'formato': 'norma', 'erratas': 'no aplica (release vigente)',
        'usos': ['VAL-01', 'Ec. 5 (p_sat)', 'VAL-02'], 'datos': {},
        'transcripciones': {'coeficientes n1-n10 de la ec. 30': ('exp05_condensacion_tiempo_real', '_IF97',
                                                                '56c886dc74a70f3b36a0ef1f595c8852b49917d0cd9a27e55f692d9625baf2c0')},
        'consistencia': ['la implementación reproduce los 3 valores de la tabla 35 a 1e-8 (grupo H)'],
        'hallazgos': [], 'decision': 'se usa', 'pedigri': 'A',
    },
    'NIST_agua': {
        'ambito': 'paper 1',
        'cita': 'NIST Chemistry WebBook, SRD 69, agua (IAPWS-95, Wagner & Pruss 2002) y N2, H2: calor latente y '
                'calores específicos a 25-40 C',
        'doi': '10.18434/T4D303', 'arbitrada': True, 'acceso': 'abierto', 'archivo': None, 'sha256': None,
        'formato': 'base de datos de referencia', 'erratas': 'no aplica (base de datos de referencia)',
        'usos': ['VAL-01', 'q_rech (calor a rechazar)'],
        'datos': {'L_v_agua': (2441.7e3, 'h_fg a 25 C'), 'cp_agua_liq': (4180, 'líquido 25-40 C'),
                  'cp_vapor': (1920, 'vapor saturado 25-40 C'), 'cp_N2': (1041, 'N2 a 27 C'),
                  'cp_H2': (14300, 'H2 a 15-27 C')},
        'modelo': ('parametros', None),
        'consistencia': ['PCI + (M_H2O/M_H2) L_v = PCS de CODATA (grupo H)'], 'hallazgos': [],
        'decision': 'se usa', 'pedigri': 'A',
    },
    'IEA': {
        'ambito': 'paper 1',
        'cita': 'IEA, Global Hydrogen Review 2025, Assumptions annex, p. 6 (electrólisis 63 % del PCI en 2024)',
        'doi': 'no tiene (reporte de la IEA)', 'arbitrada': False,
        'acceso': 'público (reporte institucional, no revisado por pares)', 'archivo': None, 'sha256': None,
        'formato': 'reporte', 'erratas': 'no aplica', 'usos': ['VAL-01', 'regeneración (EXP-11)'],
        'datos': {'eta_electrolisis': (0.63, 'p. 6 del anexo')},
        'modelo': ('parametros', None),
        'consistencia': ['63 % del PCI = 52.9 kWh/kg, dentro del rango publicado de electrolizadores PEM/alcalinos'],
        'hallazgos': ['reporte institucional, no artículo revisado: se califica B'],
        'decision': 'se usa', 'pedigri': 'B',
    },
    # ── Datos de referencia (pedigrí A), guardados en modelo/datos_fuente (en git) ──
    'NIST_H2': {
        'ambito': 'paper 1',
        'cita': 'NIST Chemistry WebBook, SRD 69, propiedades de fluidos del H2 normal (EOS de Leachman et al. '
                '2009; viscosidad de Muzny et al. 2013; conductividad de Assael et al. 2011), isobara 0.101325 MPa',
        'doi': '10.18434/T4D303', 'arbitrada': True, 'acceso': 'abierto (consulta 2026-09-25)',
        'archivo': 'modelo/datos_fuente/nist/H2_isobara_0.101325MPa_280-320K.txt',
        'sha256': '146048e6c4c8e1470a8d360e08efc437185d60bcb9ad43a32f57235da0570e0e',
        'formato': 'tabla de texto del servicio de NIST', 'erratas': 'no aplica (base de datos de referencia)',
        'usos': ['VAL-03 (convección natural interna del gas)'],
        'datos': {'k_300K_W_mK': (0.18670, 'fila 300 K'), 'mu_300K_Pa_s': (8.9385e-06, 'fila 300 K'),
                  'cp_300K_J_gK': (14.313, 'fila 300 K'), 'rho_300K_kg_m3': (0.081841, 'fila 300 K')},
        'consistencia': ['rho de NIST vs gas ideal a 300 K y 1 atm: diferencia < 0.1 % (H2 casi ideal)',
                         'cp de NIST 14.313 J/(g K) vs cp_H2 del modelo 14.30 J/(g K)'],
        'hallazgos': [],
        'decision': 'se usa', 'pedigri': 'A',
    },
    'JANAF_Al': {
        'ambito': 'paper 1',
        'cita': 'Chase, M. W., NIST-JANAF Thermochemical Tables, 4th ed., J. Phys. Chem. Ref. Data Monograph 9 '
                '(1998), tabla Al-001 (aluminio)',
        'doi': '10.18434/T42S31', 'arbitrada': True, 'acceso': 'abierto (janaf.nist.gov, consulta 2026-09-25)',
        'archivo': 'modelo/datos_fuente/janaf/Al-001.txt',
        'sha256': 'ce4bb740c061a75998b3b4262bbdabf76c010360a7eeb73188476005878a1e58',
        'formato': 'tabla de texto', 'erratas': 'no aplica (tabla de referencia)',
        'usos': ['VAL-03 (tau = c/k del radiador)'],
        'datos': {'Cp_298K_J_molK': (24.209, 'fila 298.15 K')},
        'consistencia': ['24.209 J/(mol K) / 26.98 g/mol = 897 J/(kg K), el valor de manual para Al a 25 C'],
        'hallazgos': [],
        'decision': 'se usa', 'pedigri': 'A',
    },
    # ── Tecnología, ambiente y geometría del paper 1 (2026-09-25): los datos de clase de tecnología
    #    del balance de masa, la regeneración y la inercia. Cada dato confirmado contra la IMAGEN de
    #    su página cuando el texto extraído salía corrido.
    'FZ22': {
        'ambito': 'paper 1',
        'cita': 'Aerospace Technology Institute, FlyZero: Thermal Management Roadmap Report, FZO-PPN-COM-0019, '
                'ATI, Cranfield, UK, marzo 2022',
        'doi': 'no tiene (reporte del ATI)', 'arbitrada': False,
        'acceso': 'abierto (ati.org.uk); hoja de ruta de la industria, no revisada por pares',
        'archivo': REF + 'FlyZero_2022_FZO-PPN-COM-0019_Thermal_Management_Roadmap.pdf',
        'sha256': 'f2bbcd5ba7e81441201b61f17bc3a2254f5c738b398a0cdd4202b4857ecb7e51',
        'formato': 'texto nativo; las tablas son gráficos: números confirmados contra la imagen (pp. 4, 5, 7)',
        'erratas': 'no aplica (reporte sin DOI)',
        'usos': ['P9: masa del rechazo (kWkg_TMS_sistema, desde 2026-09-27; perdida_rad)',
                 'EXP-14 y tabla S8 (kWkg_radiador, temperaturas de rechazo)',
                 'VAL-03 (tau = c_Al/k)', 'paper 1, secs. 3.3 y 5'],
        'datos': {'kWkg_radiador': ((10.0, 20.0), 'p. 5, "Air radiator", specific heat rejection: 10 (2026, LT 80 C), '
                                                  '20 (2035, HT 160 C)'),
                  'perdida_rad': ((0.20, 0.15), 'p. 5, "Air radiator", power loss: 20 % (2026), 15 % (2035)'),
                  'kWkg_TMS_sistema': ((5.0, 10.0), 'p. 4, "FC thermal management system specific heat rejection" '
                                       '(heat rejection rate / thermal management total system mass): 5 (2026, '
                                       'rechazo a 100 C, compresion de vapor), 10 (2035, LT-PEM, 100 C, liquido)')},
        'modelo': ('parametros', None),
        'transcripciones': {'temperaturas de rechazo 100 / 160 C (p. 4)': (
            'exp14_calor_por_temperatura', 'T_RECHAZO_FLYZERO',
            'a990abe785383cd67b807efa24aa9699a9e259828a3f993e9d453e353f78d041'),
                            'banda caliente 20 / 20 kW/kg (p. 7)': (
            'exp14_calor_por_temperatura', 'KWKG_CALIENTE',
            'a7e6c71035b7481b5923b5f5e83e0a12665bc40241fb33eeabe4745833761cc7')},
        'datos_no_usados': {'radiador de aire 2050 [kW/kg], pérdida': ((30, 0.10), 'p. 5'),
                            'intercambiador aceite/combustible del ciclo de turbina de H2 [kW/kg]':
                                ('20 (2026), >20 (2035), >20* (2050)', 'p. 7'),
                            'sistema térmico HT-PEM 2035 [kW/kg]': (20, 'p. 4, rechazo a 160 C')},
        'consistencia': ['p. 4 y p. 5 coinciden: rechazo a 100 C en 2026 (LT-PEM 80 C + ~20 C con compresión de '
                         'vapor, notas de la p. 4) y a 160 C en 2035 (HT-PEM)',
                         'pérdida de 20/15/10 % decreciente con el kW/kg creciente: coherente (radiador más '
                         'caliente, menos aire por kW)'],
        'hallazgos': ['son PROYECCIONES de tecnología (2026, 2035), no equipos medidos: el paper debe decirlo '
                      '(PIRT P9)',
                      'EXP-14 usaba 20 y 30 kW/kg para la banda caliente; FlyZero (p. 7) da 20 en 2026 y solo ">20" '
                      'en 2035: el 30 no estaba impreso (el 30 de la p. 5 es el radiador de aire de 2050). '
                      'Corregido el 2026-09-25 (decisión del autor): 20 kW/kg en los dos años',
                      'AUDITORIA 2026-09-27: I5, los 20 kW/kg de la p. 7 son intercambiadores aceite/combustible del '
                      'ciclo de turbina de H2 (el sumidero es el combustible), no aplican a rechazar al aire; S1, los '
                      'indicadores son de aviones en vuelo (aire de impacto) y la calibracion a 30 C no esta en la '
                      'fuente. Paper 1: masa del rechazo con el indicador de SISTEMA de la p. 4 (caja negra)'],
        'decision': 'se usa (proyección declarada): paper 1 con el indicador de sistema de la p. 4 y la pérdida de la '
                    'p. 5',
        'pedigri': 'B',
    },
    'CHIN20': {
        'ambito': 'paper 1',
        'cita': 'Chin, J. C., Aretskin-Hariton, E. D., Ingraham, D. J., Hall, D. L., Schnulo, S. L., Gray, J. S., '
                'Hendricks, E. S., Battery evaluation profiles for X-57 and future urban electric aircraft, AIAA '
                'Propulsion and Energy 2020 Forum (AIAA/IEEE EATS), paper AIAA 2020-3567',
        'doi': '10.2514/6.2020-3567', 'arbitrada': False,
        'acceso': 'abierto (NTRS 20205004497, versión del 15 jul 2020); artículo de congreso, no de revista',
        'archivo': REF + 'Chin_2020_NTRS_20205004497_X57_Battery_Evaluation.pdf',
        'sha256': '93197c20ecfd69791c6d63597cfb5189c6f5a6fab9bca1ee01fe22235bc7aa57', 'formato': 'texto nativo',
        'erratas': 'Crossref 2026-09-25: CORRECCIÓN publicada el 16 sep 2020 (10.2514/6.2020-3567.c1), posterior a '
                   'nuestro PDF; su contenido no se consultó (AIAA, de pago)',
        'usos': ['P11: batería a nivel paquete (f_celda_paquete: buffer de EXP-11, balance_masa)'],
        'datos': {'f_celda_paquete': (149.0 / 225.0, 'p. 5: celda "nominal" 225 Wh/kg, paquete 149 Wh/kg '
                                                     '(55.3 kWh, 370.86 kg)')},
        'modelo': ('parametros', None),
        'consistencia': ['55.3 kWh / 370.86 kg = 149.1 Wh/kg (su propio cociente)',
                         '5120 celdas Samsung 30Q (HL18, p. 9) x 3.0 Ah x 3.60 V = 55.3 kWh nominales: la misma '
                         'energía; a 225 Wh/kg son 246 kg de celdas y 125 kg de empaque, el "approximately 125 kg" '
                         'que ellos dan; 246 / 370.86 = 0.663',
                         'contraste: su presentación de NASA del 28 jul 2020 (NTRS 20205005267, diap. 6; local: '
                         'Chin_2020_NTRS_20205005267_presentacion.pdf) repite '
                         '"Pack *0.663 (structure) -> 149 Wh/kg"',
                         'Pattanayak 2025 (PAT25, p. 7) da los mismos 225 -> 149 Wh/kg (66 %)'],
        'hallazgos': ['parametros.py citaba como primaria NTRS 20180005737 (HL18): ese documento NO trae 225 ni 149 '
                      'Wh/kg; la primaria es esta. Cita corregida',
                      'Pattanayak cita estos números a Chin et al. 2019 (AIAA SciTech 2019-0784), que no se tiene; '
                      'están aquí, en la p. 5',
                      'es energía NOMINAL de celda; la usable es menor (47 de 55.3 kWh en HL18). El modelo usa la '
                      'fracción de MASA celda/paquete, que no depende de esa base',
                      'un solo paquete (un diseño de vuelo): se toma como representativo de la clase'],
        'decision': 'se usa (fracción de masa celda/paquete 0.662); pendiente: ver la corrección de AIAA',
        'pedigri': 'B',
    },
    'HL18': {
        'ambito': 'paper 1',
        'cita': 'Hernandez Lugo, D., Clarke, S., Miller, T., Redifer, M., Foster, T., X-57 Maxwell battery from cell '
                'level to system level design and testing, presentación, NASA NTRS 20180005737, 2018',
        'doi': 'no tiene (NTRS 20180005737)', 'arbitrada': False, 'acceso': 'abierto (NTRS)',
        'archivo': REF + 'NASA_X57_NTRS_20180005737.pdf',
        'sha256': '6b4c7f193712796283f79545f998f34dcb3a78f9b5e2191c9065b5588407a4e9',
        'formato': 'texto nativo (diapositivas)', 'erratas': 'no aplica (presentación)',
        'usos': ['contraste de CHIN20 (el modelo no usa sus números)'],
        'datos': {'n_celdas': (5120, 'p. 9: "(2) 20P128S 5120 cells", Samsung INR18650-30Q, 3.0 Ah, 3.60 V'),
                  'energia_kWh': (47.0, 'p. 9: "Energy = 47.0 kWh" (efectiva)'),
                  'masa_asignada_kg': (386.0, 'p. 9: "Battery Overall Mass Allocation = 386 kg"')},
        'consistencia': ['5120 x 3.0 Ah x 3.60 V = 55.3 kWh nominales = la energía nominal de CHIN20'],
        'hallazgos': ['antes figuraba en parametros.py como la primaria del 225 -> 149 Wh/kg y no trae esos '
                      'números; su 47 kWh / 386 kg = 122 Wh/kg es energía EFECTIVA sobre masa ASIGNADA de diseño '
                      '(otra base, no comparable con 149)'],
        'decision': 'solo contraste', 'pedigri': 'C',
    },
    'PAT25': {
        'ambito': 'paper 1',
        'cita': 'Pattanayak, T., Mavris, D., Battery technology for sustainable aviation: a review of current trends '
                'and future prospects, Appl. Energy 397 (2025) 126356',
        'doi': '10.1016/j.apenergy.2025.126356', 'arbitrada': True, 'acceso': 'abierto (CC BY-NC)',
        'archivo': REF + 'Pattanayak_Mavris_2025_ApplEnergy_397_126356.pdf',
        'sha256': 'd029d0af69cb67aa8256eb98f07127f053fb9f40453027924b415fb6155655d9', 'formato': 'texto nativo',
        'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['P11: potencia específica de celda (buffer de EXP-11)', 'paper 1, sec. 5 (cita de la batería)'],
        'datos': {'kWkg_bateria_celda': ((1.0, 2.0), 'p. 9: "the cell-level specific power is about 1-2 kW/kg"')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'energía específica de celda [Wh/kg]': ((90, 260), 'p. 3'),
                            'X-57 celda -> paquete [Wh/kg]': ((225, 149), 'p. 7, su ref. [49]')},
        'consistencia': ['el 225 -> 149 Wh/kg de la p. 7 coincide con la primaria (CHIN20, p. 5)',
                         'la X-57 (225 Wh/kg de celda) cae dentro de su rango de celdas 90-260 Wh/kg'],
        'hallazgos': ['es una REVISIÓN (fuente secundaria): el 1-2 kW/kg no trae la medición detrás',
                      'Whkg_bateria_celda = 200 Wh/kg del modelo NO es un número de la fuente: es un valor elegido '
                      'dentro de su rango 90-260 (debajo de los 225 de la X-57); se usa en el chequeo de energía del '
                      'buffer de EXP-11: declararlo como supuesto',
                      'el comentario de parametros.py cita además "1.8 kW/kg, IEEE 10598954", sin ficha'],
        'decision': 'se usa (rango de potencia de celda, secundario)', 'pedigri': 'B',
    },
    'CT7': {
        'ambito': 'paper 1',
        'cita': 'EASA, Type-Certificate Data Sheet No. IM.E.010, General Electric CT7 series engines, Issue 10, '
                '19 dic 2025',
        'doi': 'no tiene (hoja de datos del certificado de tipo de EASA)', 'arbitrada': False,
        'acceso': 'público (easa.europa.eu); dato certificado, no revisión por pares',
        'archivo': REF + 'EASA_TCDS_IM.E.010_CT7.pdf',
        'sha256': 'e779c4a7f054a1586832c286bb0b7b986fe6991302140e97b452c98deb35b9ad',
        'formato': 'texto nativo; tablas desalineadas en el texto extraído: confirmadas contra la imagen (pp. 10-12)',
        'erratas': 'no aplica (se usa la edición vigente, Issue 10)',
        'usos': ['balance de masa por caja negra: unidad de conversión (balance_masa.py)'],
        'datos': {'TURBOEJE_TCDS': (((1879e3, 243.6), (1212e3, 194.6)),
                                    'potencia de despegue 5 min (SL): CT7-8 1879 kW (p. 11), CT7-2A 1212 kW (p. 12); '
                                    'peso seco, sec. 5 (p. 10): CT7-8 243.6 kg, CT7-2A 194.6 kg')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'potencia máxima continua (SL) [kW]': ({'CT7-8': 1523, 'CT7-2A': 1189}, 'pp. 11-12')},
        'consistencia': ['kW = shp x 0.7457: 2520 shp = 1879 kW y 1625 shp = 1212 kW (las dos columnas de la tabla)',
                         'con la potencia continua (1523 y 1189 kW) el kW/kg baja de 7.71 y 6.23 a 6.25 y 6.11',
                         'las quemas por escalón de EXP-10 duran 18-157 s, menos que los 5 min del régimen de '
                         'despegue: el régimen aplica'],
        'hallazgos': ['el texto extraído corre UNA fila las etiquetas de la sec. 5: leído del texto, CT7-8 saldría '
                      '239.5 kg y CT7-2A 355.2 kg; la imagen de la p. 10 confirma 243.6 y 194.6',
                      'es el régimen de DESPEGUE (5 min), no el continuo: válido para quemas cortas; para operación '
                      'continua usar 1523 y 1189 kW'],
        'decision': 'se usa (régimen de despegue, por la duración de las quemas)', 'pedigri': 'B',
    },
    'PV18': {
        'ambito': 'paper 1',
        'cita': 'Pande, D., Verstraete, D., Impact of solar cell characteristics and operating conditions on the '
                'sizing of a solar powered nonrigid airship, Aerosp. Sci. Technol. 72 (2018) 353-363',
        'doi': '10.1016/j.ast.2017.11.013', 'arbitrada': True, 'acceso': 'pago (entregado por el autor)',
        'archivo': REF + 'Pande_Verstraete_2018_AST_72_353.pdf',
        'sha256': '4a454e44bad46682940814004200cb699b43bea6675649c220b62aa1d16ee14f', 'formato': 'texto nativo',
        'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['regeneración solar: paneles a-Si (balance_masa.py)', 'paper 1, sec. 6'],
        'datos': {'a-Si': ((0.08, 0.120), 'p. 5 (texto) y p. 9 (pie de la fig. 11): silicio amorfo de película '
                                          'delgada, 8 %, 120 g/m2, su ref. [8]')},
        'modelo': ('parametros', 'PV_TEC'),
        'datos_no_usados': {'monounión / triple unión (eficiencia, g/m2)': (((0.24, 253), (0.32, 345)), 'p. 5')},
        'consistencia': ['el mismo 8 % y 120 g/m2 en el texto (p. 5) y en el pie de la fig. 11 (p. 9)'],
        'hallazgos': ['dato de SEGUNDA mano: lo toman de Colozza & Dolce 2003 (reporte NASA, su ref. [8]); no se '
                      'tiene la primaria',
                      '120 g/m2 es la masa del ARREGLO de película delgada, sin encapsulado ni adhesivo; el 2.4 kg/m2 '
                      'de CIGS (MIA14) es un MÓDULO comercial con adhesivo: el paper compara las dos tecnologías a '
                      'distinto nivel de integración (declararlo)'],
        'decision': 'se usa (a-Si 8 %, 0.12 kg/m2, como arreglo sin encapsulado)', 'pedigri': 'B',
    },
    'MIA14': {
        'ambito': 'paper 1',
        'cita': 'MiaSolé Hi-Tech Corp., MiaSolé FLEX module benefits (hoja comercial FLEX-02W / FLEX-02N), Santa '
                'Clara, CA (URL de 2014; el PDF lleva (c) 2016)',
        'doi': 'no tiene (hoja comercial)', 'arbitrada': False, 'acceso': 'público (miasole.com)',
        'archivo': REF + 'MiaSole_2014_FLEX_Benefits.pdf',
        'sha256': 'c8f60f5dca3b0c3becc1b8eca7c9584f7fdae1f29bb279927edad0328b1704a2', 'formato': 'texto nativo',
        'erratas': 'no aplica (hoja comercial)',
        'usos': ['regeneración solar: módulos CIGS flexibles (balance_masa.py)', 'paper 1, sec. 6'],
        'datos': {'CIGS': ((0.146, 2.4), 'p. 2: FLEX-02W, 380 W (el mejor) sobre 2.598 m x 1.000 m = 14.6 % a '
                                         '1000 W/m2 sobre el área del módulo, y 2.4 kg/m2 con adhesivo (p. 1: '
                                         '"less than 2.4 kg/m2"; apertura "up to 17%")')},
        'modelo': ('parametros', 'PV_TEC'),
        'datos_no_usados': {'FLEX-02W: largo x ancho [m], masa [kg], potencia [W]': ((2.598, 1.000, 6.2, (340, 380)),
                                                                                       'p. 2'),
                            'FLEX-02N: largo x ancho [m], masa [kg], kg/m2': ((2.598, 0.370, 2.7, 2.9), 'p. 2')},
        'consistencia': ['FLEX-02W: 6.2 kg / 2.60 m2 = 2.39 kg/m2 (cuadra con 2.4)',
                         'FLEX-02W: 340-380 W / 2.60 m2 = 131-146 W/m2 a 1000 W/m2: eficiencia sobre el área TOTAL '
                         'del módulo 13.1-14.6 %, no 17 %',
                         'FLEX-02N: 2.7 kg / 0.96 m2 = 2.8 kg/m2 (la hoja dice 2.9): mayor que el "less than 2.4" '
                         'de la p. 1'],
        'hallazgos': ['el 17 % es eficiencia de APERTURA ("up to", solo las celdas); el modelo la aplicaba al área '
                      'total cubierta y sobrestimaba la energía 16-30 %. Corregido el 2026-09-25 (decisión del autor): '
                      '14.6 %, el mejor módulo sobre su área total',
                      'hoja comercial (techos, camiones), no aeronáutica; 2.4 kg/m2 incluye el adhesivo: es cota alta '
                      'de masa para CIGS; el paper dice "up to 2.4" y el FLEX-02N pesa 2.9',
                      'año: la URL es de 2014/09 y el PDF dice (c) 2016'],
        'decision': 'se usa con corrección: eficiencia de módulo 0.146 (el mejor FLEX-02W), no el 0.17 de apertura',
        'pedigri': 'B',
    },
    'TK26': {
        'ambito': 'paper 1',
        'cita': 'Tuckerman, L. B., Inertia factors of ellipsoids for use in airship design, NACA Report No. 210, '
                'Washington, DC, 1926',
        'doi': 'no tiene (reporte NACA; NTRS 19930091277)', 'arbitrada': True,
        'acceso': 'abierto (NTRS); reporte técnico con revisión interna del NACA',
        'archivo': REF + 'Tuckerman_1926_NACA_TR210.pdf',
        'sha256': 'b688633061461511dc95a6cdf7b55e26a0ecbbf9a8619f0d21fc571fe2910be3',
        'formato': 'ESCANEO con OCR (NTRS): el OCR destruye las fórmulas; se leyeron de la imagen',
        'erratas': 'no aplica (reporte de 1926 sin DOI)',
        'usos': ['P4: masa añadida vertical k (inercia.py, EXP-13, todas las simulaciones)', 'paper 1, sec. 3.1'],
        'datos': {'k_masa_anadida': (0.86, 'p. 9, sec. 4 (esferoide prolato): k2 = beta0 / (2 - beta0) con sus '
                                           'integrales alfa0, beta0; evaluado por nosotros en esbeltez 4 = 0.860')},
        # desde 2026-09-28 el modelo evalúa la FÓRMULA con la esbeltez de cada nave (parametros.k_lamb): el grupo W
        # de verificacion.py la compara con 0.804, 0.860 y 0.917 en esbeltez 3, 4 y 6 (ya no hay un k fijo que comparar)
        'consistencia': ['alfa0 + 2 beta0 = 2 exacto en esbeltez 3, 4 y 6 (identidad de las integrales)',
                         'k2 = 0.804, 0.860, 0.917 en esbeltez 3, 4, 6: el rango 0.80-0.92 del paper',
                         'k1 (longitudinal) = 0.122, 0.082, 0.045: los valores clásicos (Lamb 1932, art. 373)',
                         'contraste con vuelo: Yamasaki & Goto 2003 (YG03) medido/analítico 0.77-1.20'],
        'hallazgos': ['la fuente da FÓRMULAS, no tablas: el 0.86 es nuestro cálculo con sus fórmulas, no un número '
                      'impreso (el paper debe decir "evaluado con")'],
        'decision': 'se usa', 'pedigri': 'A',
    },
    'SEN03': {
        'ambito': 'paper 1',
        'cita': 'SENAMHI, MINEM, Atlas de Energía Solar del Perú, SENAMHI, Lima, 2003',
        'doi': 'no tiene (atlas institucional)', 'arbitrada': False,
        'acceso': 'público (senamhi.gob.pe); mediciones con piranómetro de la red nacional',
        'archivo': REF + 'SENAMHI_2003_Atlas_Energia_Solar_Peru.pdf',
        'sha256': 'ea93f6aded4ab55b38f2a4f421237cbca2618f62f6a56013bd0b21c04c6d758e',
        'formato': 'texto nativo (tablas); los anexos son gráficos, leídos de la imagen',
        'erratas': 'no aplica (atlas sin DOI)',
        'usos': ['regeneración solar: irradiancia medida en Lima (irradiancia_casco.py, balance_masa.py)',
                 'paper 1, sec. 6'],
        'datos': {'SENAMHI_GHI_ANUAL': (4.06, 'p. 16, tabla 2: estación A. von Humboldt (La Molina, Lima; '
                                              'piranómetro 1968-1999)'),
                  'SENAMHI_LAT': (-12.08, 'p. 16, tabla 1 (12.08 S, 76.95 W, 238 m)')},
        'modelo': ('parametros', None),
        'transcripciones': {'H/H0 mensual de A. von Humboldt (anexo 6, p. 29, leído +-0.01)': (
            'parametros', 'SENAMHI_KT_MES', '1832d313d3e724f9e378cd3d5209a49c13373a6813ca4ab80494b667687acae7')},
        'datos_no_usados': {'heliofanía relativa': (0.40, 'p. 16, tabla 2')},
        'consistencia': ['anual rearmado con sus anexos 6 y 7 (12 meses de H/H0 x H0 a 12 S): 4.09 kWh/m2/día, '
                         '+0.8 % sobre la tabla 2 (dentro de la lectura de la figura; grupo T)',
                         'H0 calculado a 12 S: 28.2 MJ/m2 en junio y 40.0 en diciembre, lo que muestra su anexo 7 '
                         '(p. 31; grupo T)',
                         'invierno (jun-ago) 2.42, verano (dic-feb) 5.62 kWh/m2/día'],
        'hallazgos': ['HISTORIA: el modelo usaba rendimiento_efectivo = 2.0 kWh/m2/día (1.2-2.8) citando a SENAMHI, '
                      'que no salía de esta fuente (la estación mide 4.06 de promedio anual). Reemplazado el '
                      '2026-09-25 (decisión del autor) por la irradiancia medida x el factor de forma del casco',
                      'el GHI medido ya incluye la nubosidad; lo que faltaba era el factor de forma de la mitad '
                      'superior del casco (0.81-0.84, irradiancia_casco.py)',
                      'H/H0 mensual leído de un gráfico (+-0.01): se registra por su huella'],
        'decision': 'se usa (irradiancia medida por mes; invierno y verano en el paper)', 'pedigri': 'B',
    },
    'PAGE61': {
        'ambito': 'paper 1',
        'cita': 'Page, J. K., The estimation of monthly mean values of daily total short wave radiation on vertical '
                'and inclined surfaces from sunshine records for latitudes 40N-40S, in: Proceedings of the United '
                'Nations Conference on New Sources of Energy, Rome, 21-31 August 1961, vol. 4 (Solar Energy I), '
                'United Nations, New York, 1964, pp. 378-390 (trabajo E/CONF.35/S/98)',
        'doi': 'no tiene (actas de la ONU, símbolo E/CONF.35/5)', 'arbitrada': False,
        'acceso': 'ORIGINAL consultado 2026-09-27 en la Biblioteca Digital de la ONU (digitallibrary.un.org/record/'
                  '3824017: vol. 4 completo, 676 pp., 48.8 MB, escaneado con OCR). La descarga directa la bloquea el '
                  'filtro anti-bot del sitio: se leyó en el navegador, sin copia local. Antes (2026-09-25) solo se '
                  'tenía de segunda mano: Duffie & Beckman (DB13) y el paquete solaR de CRAN (FdKtPage)',
        'archivo': None, 'sha256': None, 'formato': 'escaneo con OCR; la ec. [6] (p. 383) y la tabla 3 (p. 382) '
                                                      'confirmadas contra la imagen', 'erratas': 'no aplica',
        'usos': ['fracción difusa de las medias mensuales (irradiancia_casco.py)'],
        'datos': {'PAGE_B': (1.13, 'p. 383, ec. [6]: "The mean regression equation for the ten stations is '
                                   'D_H/Q = 1.00 - 1.13 Q/Q_o" (medias mensuales)')},
        'modelo': ('irradiancia_casco', None),
        'datos_no_usados': {'tabla 3: c y d de las 10 estaciones': ('c 0.72-1.20, d -0.67 a -1.43; medias 1.00 y -1.13',
                                                                   'p. 382 (confirmada contra la imagen; medias '
                                                                   'recalculadas)'),
                            'máximo de la difusa [fracción de Q_o]': (0.221, 'p. 383, en Q/Q_o = 0.44')},
        'consistencia': ['medias de la tabla 3 recalculadas: c = 10.00/10 = 1.00, d = -11.27/10 = -1.127: es la ec. [6]',
                         'el máximo de D_H = c Q + d Q^2/Q_o con c = 1.00 y d = -1.13 cae en Q/Q_o = 1/2.26 = 0.44 y vale '
                         '1/4.52 = 0.221 Q_o, lo que dice el trabajo (p. 383)',
                         'a K_T 0.29-0.52 da f_d 0.41-0.67; Liu & Jordan (medias mensuales, mismo código de solaR) '
                         'da 0.36-0.60: el factor de forma cambia < 1 % entre las dos',
                         'con toda la luz directa o toda difusa el factor va de 0.78 a 0.86: la correlación elegida '
                         'casi no importa (grupo T)'],
        'hallazgos': ['el trabajo es de 1961 (conferencia de Roma, documento del 16 de mayo de 1961) pero las actas '
                      'se publicaron en 1964: se cita el volumen de 1964, pp. 378-390 (la cita vieja decía 1961, p. 378)',
                      'Page aclara (p. 383) que la ecuación, como la de Angstrom, no vale para valores diarios sino '
                      'para medias mensuales: es el uso del modelo',
                      'los meses de invierno de Lima (K_T 0.29) están dentro de los datos de Page: la columna "Range of '
                      'monthly mean values %" de la tabla 3 (se lee como Q/Q_o) va de 24-45 (Kew) y 32-52 '
                      '(Leopoldville) a 57-78 (Windhoek); Uccle, ajustada a ojo, llega a 9'],
        'decision': 'se usa', 'pedigri': 'B',
    },
    'DB13': {
        'ambito': 'paper 1',
        'cita': 'Duffie, J. A., Beckman, W. A., Solar Engineering of Thermal Processes, 4th ed., Wiley, 2013: '
                'declinación (Cooper 1969), H0 diaria, días representativos (Klein 1977), cielo isotrópico sobre '
                'superficies inclinadas (Liu & Jordan)',
        'doi': '10.1002/9781118671603', 'arbitrada': True, 'acceso': 'libro de texto (no local)',
        'archivo': None, 'sha256': None, 'formato': 'libro de texto', 'erratas': 'no aplica',
        'usos': ['irradiancia_casco.py: H0, geometría solar y modelo isotrópico'],
        'datos': {'G_SC': (1367.0, 'constante solar de los atlas de esa época (la 4.a ed. usa 1361: -0.4 % en H0)'),
                  'albedo_suelo': (0.2, 'reflectancia del suelo, valor estándar')},
        'modelo': ('parametros', None),
        'consistencia': ['H0 a 12 S = anexo 7 de SENAMHI (grupo T)',
                         'la extensión al casco curvo es nuestra: forma cerrada (|s_p| + s_z)/pi verificada contra '
                         'integración directa y contra dos casos analíticos (grupo T)'],
        'hallazgos': ['el albedo 0.2 es un SUPUESTO (aporta 3.6 % del GHI); sobre el mar sería menor'],
        'decision': 'se usa', 'pedigri': 'A',
    },
    'DOE17': {
        'ambito': 'paper 1',
        'cita': 'U.S. DOE, Fuel Cell Technologies Office, Target Explanation Document: Onboard Hydrogen Storage for '
                'Light-Duty Fuel Cell Vehicles, 2017',
        'doi': 'no tiene (documento del DOE)', 'arbitrada': False, 'acceso': 'público (energy.gov)',
        'archivo': REF + 'DOE_2017_Onboard_H2_Storage_Targets_Explanation.pdf',
        'sha256': '4b39a1eb85b4eb8c8aa69eeecca12c00ff7760ccf547ef6c7ad27efbd41f1bf1', 'formato': 'texto nativo',
        'erratas': 'no aplica (documento de metas)',
        'usos': ['hasta 2026-09-29: tanques (balance_masa.py); desde entonces solo contraste de DOE24'],
        'datos': {'wt_COPV_meta_2025': (0.055, 'p. 6, tabla 1: META 2025 de capacidad gravimétrica del SISTEMA, '
                                               '1.8 kWh/kg (0.055 kg H2 / kg sistema)')},
        'datos_no_usados': {'metas 2020 / final [kg H2 / kg sistema]': ((0.045, 0.065), 'p. 6, tabla 1')},
        'consistencia': ['1.8 kWh/kg / 33.3 kWh/kg (PCI del H2) = 0.054: cuadra con el 0.055 entre paréntesis',
                         'la misma meta en las tablas de las pp. 6, 8 y 12',
                         'la misma meta 2025 (1.8 kWh/kg) en la tabla 1 del Record 24006 (DOE24)'],
        'hallazgos': ['es una META para 2025 (vehículos livianos), no un valor logrado',
                      'cuenta todo el sistema (tanque, válvulas, regulador, soportes), no solo el tanque',
                      '2026-09-29: según el propio DOE (DOE24), en 2023 el sistema daba 1.71 kWh/kg: la meta no se '
                      'alcanzó'],
        'decision': 'reemplazada por DOE24 (autor 2026-09-29: lo que existe, no la meta); queda como contraste',
        'pedigri': 'B',
    },
    'DOE24': {
        'ambito': 'paper 1',
        'cita': 'Simmons, K., Hulvey, Z., Celestine, A.-D., Onboard Type IV compressed hydrogen storage system - '
                'cost and performance status, DOE Hydrogen Program Record 24006, U.S. DOE, 28 abr 2025 (actualiza '
                'el Record 19008)',
        'doi': 'no tiene (Program Record del DOE)', 'arbitrada': False,
        'acceso': 'público (hydrogen.energy.gov/docs/hydrogenprogramlibraries/pdfs/24006-onboard-storage-cost-'
                  'performance-status.pdf); revisado por B. D. James (SA), D. Byerly (Hexagon Agility), R. Ahluwalia '
                  '(ANL), D. Hill (OPmobility), T. Abdel-Baset y M. Gras (Forvia); no es revisión por pares de revista',
        'archivo': REF + 'DOE_2025_Record24006_TypeIV_Storage_Status.pdf',
        'sha256': 'b08885e9d3a967f1b61f3013e42c73b2325fdabeec660b1122e569444d5f082b', 'formato': 'texto nativo',
        'erratas': 'no aplica (Program Record); consultado 2026-09-29',
        'usos': ['balance de masa: tanques para reponer el H2 perdido y guardar el regenerado (balance_masa.py)',
                 'paper 1, sec. 6 y tabla 6'],
        'datos': {'kWhkg_tanque_sistema': (1.71, 'p. 1 ("Gravimetric Energy Density: 1.71 kWh/kg system") y tabla 1 '
                                                 '(p. 2), columna "2023 Status"; sistema de 700 bar, un tanque, 5.6 kg '
                                                 'de H2 útil'),
                  'TANQUE_DOE24': ({'P_bar': 700.0, 'P_min_bar': 15.0, 'T_C': 15.0, 'V_L': 144.2, 'D_cm': 39.6,
                                    'util_kg': 5.6, 'guardado_kg': 5.8, 'compuesto_kg': 72.7},
                                   'tabla 2 (p. 4, confirmada contra la imagen): presión nominal 700 bar, mínima 15, '
                                   '15 C, 144.2 L, diámetro interior 39.6 cm, 5.6 kg útiles y 5.8 guardados; fig. 1 '
                                   '(p. 5): compuesto 72.7 kg, valor impreso sobre la barra "2023 Baseline"')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'estados 2013 / 2015 / 2019 [kWh/kg sistema]': ((1.5, 1.40, 1.48), 'tabla 1'),
                            'metas 2025 / final [kWh/kg sistema]': ((1.8, 2.2), 'tabla 1'),
                            'capacidad volumétrica 2023 [kWh/L sistema]': (0.86, 'tabla 1'),
                            'factor de seguridad mínimo; liner 2023 (PA) [cm]': ((2.25, 0.38), 'tablas 2 y 3 (p. 4)')},
        'consistencia': ['1.71 kWh/kg / 33.3 kWh/kg (PCI) = 5.1 % de H2 en masa del sistema',
                         'el estado 2015 de la tabla (1.40) = el del Record 15013 (1.40 +- 0.04 kWh/kg, Ordaz, Houchins '
                         'y Hua, 30-09-2015; PDF en referencias/DOE_2015_Record15013_TypeIV_Storage_Status.pdf)',
                         '"~16 % increase since 2019": 1.71 / 1.48 = 1.155'],
        'hallazgos': ['es una PROYECCIÓN por análisis (elementos finitos y DFMA de ANL y SA) para 500 000 o 100 000 '
                      'unidades al año, no la masa pesada de un tanque de serie',
                      'es de 700 bar y el modelo guarda a 350 bar: parametros.tanque_a_presion lleva ESTE diseño a 350 '
                      'con sus propias reglas (diámetro fijo, compuesto ~ P V, liner ~ área con domos hemisféricos, '
                      'componentes de DOE19 iguales, densidades NIST a 15 C): 1.80 kWh/kg, razón 1.052 (ANL-10/24, '
                      'HUA10, da 1.058)',
                      'la meta 2025 del DOE (1.8, DOE17) no se alcanzó'],
        'decision': 'se usa (estado y diseño 2023 del sistema completo, llevado a 350 bar)', 'pedigri': 'B',
    },
    'DOE19': {
        'ambito': 'paper 1',
        'cita': 'Adams, J., Houchins, C., Ahluwalia, R., Onboard Type IV compressed hydrogen storage system - cost and '
                'performance status, DOE Hydrogen and Fuel Cells Program Record 19008, U.S. DOE, 25 nov 2019',
        'doi': 'no tiene (Program Record del DOE)', 'arbitrada': False,
        'acceso': 'público (hydrogen.energy.gov/docs/hydrogenprogramlibraries/pdfs/19008_onboard_storage_cost_'
                  'performance_status.pdf); revisado por B. D. James (SA), M. Veenstra (Ford), N. Newhouse (Hexagon '
                  'Lincoln) y M. Weimer (PNNL)',
        'archivo': REF + 'DOE_2019_Record19008_TypeIV_Storage_Status.pdf',
        'sha256': 'a9abb21df6456f6265caf53c03bd02843c283b2e760ff1e92ed84b165c6079e8', 'formato': 'texto nativo; tabla suplementaria 5 confirmada contra la imagen (p. 12), 2026-09-29',
        'erratas': 'no aplica (Program Record)',
        'usos': ['balance de masa: componentes del sistema que no escalan con la presión (parametros.tanque_a_presion)'],
        'datos': {'TANQUE_DOE19_FIJOS': ({'boss': 0.9, 'valvula': 3.0, 'regulador': 3.6, 'otros_BOP': 15.4},
                                         'tabla supl. 5 (p. 12): tank boss 0.9, integrated in-tank valve 3.0, '
                                         'integrated regulator 3.6, other BOP 15.4 kg')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'H2 guardado / útil, compuesto (fibra + resina), liner HDPE, total [kg]':
                                ((5.8, 5.6, 90.3, 67.1, 23.2, 8.0, 126.1), 'tabla supl. 5 (p. 12)'),
                            'volumen del tanque 2019 [L]': (147, 'tabla supl. 5')},
        'consistencia': ['5.6 x 33.3 / 126.1 = 1.479: el 1.48 kWh/kg de la misma tabla',
                         'el diseño 2023 (DOE24, p. 3) mantiene "identical storage system design and BOP components" '
                         'de este registro: por eso sus componentes se toman fijos',
                         'el liner 2023 que sale por diferencia (7.7 kg) es menor que este (8.0 kg): cuadra con el PA '
                         'de 0.38 cm que reemplazó al HDPE de 0.5 cm (DOE24, tabla 3)'],
        'hallazgos': ['las filas de la tabla supl. 5 suman 127.0 kg (con 5.8 de H2) y el total impreso es 126.1; la '
                      'capacidad 1.48 usa el 126.1 (sin el boss cuadraría). Si el boss sobrara, el liner 2023 por '
                      'diferencia sería 0.9 kg mayor y el sistema a 350 bar ~0.5 % más pesado: no se corrige'],
        'decision': 'se usa (componentes fijos del sistema)', 'pedigri': 'B',
    },
    'HUA10': {
        'ambito': 'paper 1',
        'cita': 'Hua, T. Q., Ahluwalia, R. K., Peng, J.-K., Kromer, M., Lasher, S., McKenney, K., Law, K., Sinha, J., '
                'Technical assessment of compressed hydrogen storage tank systems for automotive applications, '
                'ANL-10/24, Argonne National Laboratory, septiembre de 2010',
        'doi': '10.2172/1010895 (informe OSTI); versión de revista: Int. J. Hydrogen Energy 36 (2011) 3037-3049, '
               'doi 10.1016/j.ijhydene.2010.11.090 (no leída: muro de pago)', 'arbitrada': False,
        'acceso': 'abierto (energy.gov/sites/prod/files/2014/03/f9/compressedtank_storage.pdf); informe de laboratorio '
                  'nacional que evalúa los tanques tipo IV de Quantum (p. 7); la versión de revista sí es arbitrada',
        'archivo': REF + 'Hua_2010_ANL-10-24_Compressed_H2_Storage.pdf',
        'sha256': 'd4772b16a705ebc39961760caf1aa7a76bce8622b817fa47ba6a6ba00109c4f9',
        'formato': 'texto nativo; la tabla I sale desalineada en el texto extraído: los dos valores se toman del texto '
                   'de la p. 16, que los da en prosa',
        'erratas': 'no aplica (informe); erratas de la versión de revista no consultadas',
        'usos': ['verificación, grupo L: contraste independiente de la razón 350/700 bar que el modelo calcula con el '
                 'diseño 2023 del DOE (parametros.tanque_a_presion); hasta el 2026-09-29 fue la razón del modelo'],
        'datos': {'wt_ANL_350_700': ((5.5, 5.2), 'p. 16: "the gravimetric capacity of the 350-bar system is 5.5 wt%" y '
                                                 '"The gravimetric capacity of the base case 700-bar system is 5.2 wt%"; '
                                                 'tabla I (p. 6); tipo IV, un tanque, 5.6 kg útiles')},
        'consistencia': ['wt% = H2 útil (5.6 kg) / masa total del sistema (nota 2, p. 14): la misma definición que la '
                         'capacidad gravimétrica de sistema del DOE (kWh/kg = wt% x 33.3)',
                         'la diferencia la explica la fibra: capa de 25.9 mm a 700 bar frente a 14.7 mm a 350 (p. 14)',
                         'CONTRASTE con el diseño 2023 (DOE24, 2026-09-29): compuesto 72.7 kg de un sistema de 109.1 kg '
                         '(5.6 x 33.3 / 1.71), 5.8 kg de H2 y 30.6 kg de resto (liner, bosses, espuma, planta). A 350 bar '
                         'y 15 C (NIST: 23.995 kg/m3; 1.2509 a 15 bar) el tanque pasa de 144.2 a 246.2 L y el compuesto, '
                         'que a igual diámetro escala con P V, a 62.1 kg. Con el resto fijo la razón 350/700 es 1.107; '
                         'con el resto creciendo con el volumen, 0.908. El 1.058 de ANL cae dentro. (Con las mismas '
                         'densidades, la tabla 2 de DOE24 se reproduce: 5.61 kg útiles y 5.79 guardados.)'],
        'hallazgos': ['análisis de 2010 (tanques de Quantum): la razón entre presiones se supone vigente para el estado '
                      '2023; por eso el modelo ya no la usa: convierte con el diseño 2023 (DOE24) y los componentes de '
                      '2019 (DOE19) y obtiene 1.052, a 0.6 % de esta'],
        'decision': 'contraste (autor 2026-09-29: usar los datos 2023; la razón 1.052 del modelo coincide con 1.058)',
        'pedigri': 'B',
    },
    'SUN24': {
        'ambito': 'paper 1',
        'cita': 'Sun, K., Jiang, Y., Liu, H., Configuration optimization of PV array for stratospheric airship under '
                'nonuniform radiation condition, Aerosp. Sci. Technol. 155 (2024) 109550',
        'doi': '10.1016/j.ast.2024.109550', 'arbitrada': True, 'acceso': 'pago (local)',
        'archivo': REF + 'Sun_Jiang_Liu_2024_AST_155_109550.pdf',
        'sha256': 'e660afb4f44a454e1302795faf094bca1020a91b93bc310b5aed5bfd7dbf0dec', 'formato': 'texto nativo',
        'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['regeneración solar: geometría del arreglo sobre el casco (A_paneles_fraccion)', 'paper 1, sec. 6'],
        'datos': {'angulo_central_arreglo_deg': ((-90, 90), 'p. 2, sec. 2 ("covers the curved upper surface") y '
                                                           'fig. 9, p. 5'),
                  'casco_largo_diametro_m': ((220, 50), 'p. 4, tabla 1')},
        'consistencia': ['con su casco (220 x 50 m) como esferoide, la mitad superior entre y = 20 y 200 m (el tramo '
                         'de irradiancia pareja de su fig. 10, p. 7) es 0.449 del área total: el 0.45 del modelo es '
                         'coherente con esa geometría'],
        'hallazgos': ['la fuente NO da la fracción 0.45 que el paper le atribuye ("covering a fraction 0.45 of the '
                      'hull [Sun 2024]"): da la geometría (mitad superior, -90 a +90 grados); 0.45 es nuestro '
                      'supuesto (mitad superior sin las puntas)',
                      'dirigible estratosférico (20 km): solo se toma la geometría del arreglo, no sus irradiancias',
                      'la clave bib del paper es "Gao2024" pero los autores son Sun, Jiang y Liu (solo la clave; el '
                      'texto impreso está bien)'],
        'decision': 'solo cualitativo (geometría); la fracción 0.45 queda como supuesto declarado', 'pedigri': 'C',
    },
    'ZH25': {
        'ambito': 'paper 1',
        'cita': 'Zhang, X., Yue, H., Zheng, H., Tan, L., Zhang, Z., Li, F., Proton exchange membrane fuel cells for '
                'aircraft applications: a comprehensive review of key challenges and development trends, Hydrogen 6 '
                '(2025) 116',
        'doi': '10.3390/hydrogen6040116', 'arbitrada': True, 'acceso': 'abierto (MDPI; entregado por el autor)',
        'archivo': REF + 'Zhang_2025_Hydrogen_6_116.pdf',
        'sha256': 'bac0f0cc5b890960b7058d45a7bf0fd3e73176f55960d1f3097dcb73791103cc', 'formato': 'texto nativo',
        'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['balance de masa: unidad de regeneración (electrolizador PEM, PROXY)', 'paper 1, tabla de masas'],
        'datos': {'kWkg_PEM': ((0.5, 2.2), 'p. 2 y tabla 1 (p. 3): pila PEM enfriada por aire 0.5-1.0 kW/kg, por '
                                           'agua 1.2-2.0 kW/kg; p. 9: sistema turboalimentado de HyPoint con NASA, '
                                           '2.2 kW/kg')},
        'modelo': ('parametros', None),
        'consistencia': ['el mismo 0.5-1.0 / 1.2-2.0 kW/kg en el texto (p. 2) y en la tabla 1 (p. 3)'],
        'hallazgos': ['son PILAS de combustible, no electrolizadores: el modelo las usa como PROXY de la masa por '
                      'kW de una celda PEM con su planta (declarado [S-proxy] en parametros.py)',
                      'el extremo optimista 2.2 kW/kg es un solo demostrador (HyPoint), citado, no medido por ellos',
                      'es una revisión (secundaria)'],
        'decision': 'se usa como proxy declarado', 'pedigri': 'C',
    },
    'VW09': {
        'ambito': 'paper 1',
        'cita': 'Verhelst, S., Wallner, T., Hydrogen-fueled internal combustion engines, Prog. Energy Combust. Sci. '
                '35 (2009) 490-527',
        'doi': '10.1016/j.pecs.2009.08.001', 'arbitrada': True,
        'acceso': 'manuscrito aceptado (133 pp., formato de preimpresión; entregado por el autor), no la versión '
                  'maquetada: las páginas citadas son las del manuscrito',
        'archivo': REF + 'Verhelst_Wallner_2009_PECS_35_490_manuscrito.pdf',
        'sha256': '65296580123c9c3cab069848995a182c38d43646daca8b9127f1385187622888', 'formato': 'texto nativo',
        'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['cierre de energía: eficiencia de motores de H2 reales (eta_motor_H2)', 'paper 1, sec. 5'],
        'datos': {'eta_motor_H2': ((0.30, 0.45), 'pp. 53-54: MAN H2876 UH01 (12.8 L, bus), 31 % al freno, '
                                                 'estequiométrico sin turbo; p. 57: inyección directa en motor Ford de '
                                                 'un cilindro, "más de 45 %" al freno estimado; conclusiones (p. 80): '
                                                 '"45 % is achievable"')},
        'modelo': ('parametros', None),
        'consistencia': ['el 45 % aparece en el cuerpo (p. 57) y en las conclusiones (p. 80); rendimientos '
                         'indicados de 42-44 % (pp. 58-59) son coherentes con un 45 % al freno solo como pico estimado'],
        'hallazgos': ['el 0.30 del modelo redondea hacia abajo el 31 % del motor de bus (conservador, 1 punto)',
                      'el 45 % es ESTIMADO (motor de un cilindro con fricción supuesta de 0.7 bar), no medido en un '
                      'motor completo: es el extremo optimista'],
        'decision': 'se usa (rango 0.30-0.45: motor de serie a pico estimado)', 'pedigri': 'B',
    },
    'CC75': {
        'ambito': 'paper 1',
        'cita': 'Churchill, S. W., Chu, H. H. S., Correlating equations for laminar and turbulent free convection from '
                'a vertical plate, Int. J. Heat Mass Transfer 18 (1975) 1323-1329',
        'doi': '10.1016/0017-9310(75)90243-4', 'arbitrada': True, 'acceso': 'pago (no local); correlación de manual',
        'archivo': None, 'sha256': None, 'formato': 'correlación empírica', 'erratas': 'Crossref 2026-09-25: ninguna',
        'usos': ['VAL-03: tiempo térmico del gas (convección natural interna), paper 1, discusión'],
        'datos': {'Nu': ('(0.825 + 0.387 Ra^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27))^2',
                         'correlación para todo Ra, placa vertical (val03_cotas.py)')},
        'consistencia': ['a Ra alto Nu crece como Ra^(1/3): h no depende de la longitud elegida (el diámetro)',
                         'h = 5.5 W/(m2 K) para el H2 con las propiedades de NIST (NIST_H2)'],
        'hallazgos': ['es de placa vertical; aplicada al interior curvo del casco es una aproximación: se cubre con '
                      'h x0.5 a x2 (tiempo térmico 13-53 min)'],
        'decision': 'se usa (orden de magnitud, con el rango de h propagado)', 'pedigri': 'B',
    },
    'FIS37': {
        'ambito': 'paper 1',
        'cita': 'Fischer, W. O. L. (Goodyear-Zeppelin Corp.), Ballast water recovery for airships, U.S. Patent '
                '2,078,532, 27 abr 1937 (solicitud 27 ene 1934)',
        'doi': 'no tiene (patente de EE. UU.)', 'arbitrada': False,
        'acceso': 'público (Google Patents); examinada por la oficina de patentes, no revisada por pares',
        'archivo': REF + 'Fischer_1937_US2078532.pdf',
        'sha256': 'e3e2e73699066d4f4e32f8512db150b8f85a8ce4526cef314f0629166da43120',
        'formato': 'ESCANEO con OCR (el OCR mezcla las dos columnas): confirmado contra la imagen de la p. 2',
        'erratas': 'no aplica (patente)',
        'usos': ['paper 1, introducción: antecedente del ciclo (quemar H2 con aire y retener su agua como lastre)'],
        'datos': {'practica_previa': ('recuperar el agua del escape de los motores como lastre "has been successfully '
                                      'applied in helium-inflated airships"', 'p. 2, col. 1'),
                  'propuesta': ('quemar hidrógeno bajo una caldera; el producto es vapor de agua mezclado con el N2 del '
                                'aire de combustión; el condensado compensa el peso del combustible líquido gastado',
                                'p. 2, col. 2, l. 15-30'),
                  'agua_por_libra': ('los gases de la caldera traen ~3 veces el agua por libra que el escape de un motor',
                                     'p. 2, col. 2')},
        'consistencia': ['el "~3 veces" cuadra con la estequiometría del modelo: fracción másica de agua en el escape '
                         '0.253 con H2 y 0.088 con octano, razón 2.88'],
        'hallazgos': ['el H2 va en bolsas aparte rodeadas de helio, o es el de la envolvente de H2; el vapor mueve '
                      'turbinas que dan potencia de propulsión',
                      'propósito: compensar lentamente el combustible de crucero; no trata tasas ni sueltas de carga, '
                      'y no trae mediciones'],
        'decision': 'se usa como antecedente (cualitativo)', 'pedigri': 'C',
    },
    'DOE9013': {
        'ambito': 'paper 1',
        'cita': 'Gardiner, M., Energy requirements for hydrogen gas compression and liquefaction as related to '
                'vehicle storage needs, DOE Hydrogen and Fuel Cells Program Record 9013, U.S. DOE, 2009',
        'doi': 'no tiene (registro del programa de hidrógeno del DOE)', 'arbitrada': False,
        'acceso': 'público (hydrogen.energy.gov); aprobado por el programa, no revisado por pares',
        'archivo': REF + 'DOE_2009_Record9013_H2_Compression.pdf',
        'sha256': '7b997a04f9747535a8ca4390c68e266b943f78a1984cf5a29e7e8eced524555c', 'formato': 'texto nativo',
        'erratas': 'no aplica (registro de 2009)',
        'usos': ['EXP-16: electricidad de compresión del H2 regenerado hacia los tanques'],
        'datos': {'W_iso_20_350': (1.05, 'p. 1: trabajo isotérmico mínimo 20 -> 350 bar (exergía con NIST)'),
                  'eta_compresion_H2': (0.52, 'p. 2: eficiencia de compresión HDSAM para 350 bar (compresor '
                                              'alternativo, 56 % isentrópica x 92 % motor)'),
                  'e_comp_medido': (3.1, 'p. 2: promedio medido en las estaciones del programa (2-4 kWh/kg sin '
                                         'atípicos; 1.7-6.4 en total)')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'1 atm -> 20 bar, mínimo [kWh/kg]': (1.02, 'p. 2'),
                            'HDSAM a 440 bar; APCI 350(440) bar [kWh/kg]': ((2.23, 2.05), 'pp. 2-3')},
        'consistencia': ['gas ideal isotérmico 20 -> 350 bar a 300 K = 0.98 kWh/kg; su 1.05 con gas real de NIST es '
                         '7 % mayor, lo esperable a 350 bar (grupo V)',
                         '1.05 / 0.52 = 2.02 kWh/kg, entre su HDSAM a 440 bar (2.23) y APCI (2.05): coherente',
                         'contraste independiente: Giuffrida & Colbertaldo 2025 (GC25) dan 2-3 kWh/kg con datos de '
                         'fabricantes'],
        'hallazgos': ['las cifras son de estaciones de servicio con llenado rápido y sobrepresión (440 bar): para un '
                      'llenado lento a bordo son cota alta; el promedio medido 3.1 se usa como extremo pesimista',
                      'documento de 2009: los compresores actuales pueden ser algo mejores'],
        'decision': 'se usa (eficiencia típica y promedio medido como rango)', 'pedigri': 'B',
    },
    'GC25': {
        'ambito': 'paper 1',
        'cita': 'Giuffrida, A., Colbertaldo, P., Overview of diaphragm compressors for hydrogen service: capacity, '
                'discharge pressure and operational challenges, J. Energy Storage 129 (2025) 117286',
        'doi': '10.1016/j.est.2025.117286', 'arbitrada': True, 'acceso': 'pago (entregado por el autor)',
        'archivo': REF + 'Giuffrida_Colbertaldo_2025_JEST_129_117286.pdf',
        'sha256': '01fe866e18a59211309af802f377e51eca153124c23f2c798ef6fdb5996242cd', 'formato': 'texto nativo',
        'erratas': 'Crossref 2026-09-27: ninguna',
        'usos': ['EXP-16: contraste de la electricidad de compresión (DOE9013)'],
        'datos': {'e_comp_kWhkg': ((2.0, 3.0), 'p. 13: "2-3 kWh/kgH2" con potencia nominal y flujo de fabricantes, '
                                              'compresores de más de 30 kW; "around a 5 % increase" sobre la '
                                              'electrólisis'),
                  'tabla3_caso1': ((5.5, 0.5, 15.1, 2.6), 'p. 13, tabla 3: 5.5 Nm3/h, 0.5 -> 15.1 MPa, 2.6 kW '
                                                          'nominales, 2 etapas (confirmada contra la imagen)')},
        'datos_no_usados': {'eficiencia isentrópica medida (con aire) [%]': ((42, 55), 'p. 13, Zhao et al.'),
                            'capacidad máxima actual [kg/h]': (180, 'p. 13 (2000 Nm3/h)')},
        'consistencia': ['recalculado desde su tabla 3 (flujo x 0.0899 kg/Nm3, potencia nominal): casos 3-6 dan '
                         '2.0-3.0 kWh/kg, lo que declaran',
                         'contraste: coincide con el 2.0-3.1 kWh/kg del DOE (DOE9013) que usa el modelo (grupo V)'],
        'hallazgos': ['NO trae masas: los fabricantes no publican datos detallados ("rationally confidential", '
                      'p. 13); la masa del compresor sigue sin fuente liviana',
                      'el compresor chico de la tabla (caso 1, 0.49 kg/h; el pedido es ~1.7 kg/h) da 5.3 '
                      'kWh/kg con su potencia NOMINAL: los motores chicos van sobredimensionados y el consumo real '
                      'no se publica; un compresor pequeño podría gastar más que 2-3 kWh/kg'],
        'decision': 'contraste de la energía de compresión (confirma el rango del DOE); no resuelve la masa',
        'pedigri': 'B',
    },
    'BAUER': {
        'ambito': 'paper 1',
        'cita': 'BAUER Kompressoren, H-Series 2, compresor de hidrógeno H15.2-20 (hoja de datos en la web del '
                'fabricante, consulta 2026-09-27)',
        'doi': 'no tiene (catálogo comercial)', 'arbitrada': False,
        'acceso': 'público (bauercomp.com/products/hydrogen/h-series-2); sin PDF: datos transcritos de la página',
        'archivo': None, 'sha256': None, 'formato': 'página web del fabricante', 'erratas': 'no aplica (catálogo)',
        'usos': ['EXP-16: contraste (compresor industrial más grande que el pedido)'],
        'datos': {'COMPRESOR_BAUER': ({'P_bar': 350.0, 'flujo_kgph': (1.78, 2.22), 'motor_kW': 16.0,
                                     'masa_kg': (1066.0, 1814.0)},
                                    'tabla técnica del modelo H15.2-20: 350 bar, 1.78-2.22 kg/h, motor 16 kW, '
                                    '1066-1814 kg')},
        'modelo': ('parametros', None),
        'consistencia': ['el flujo pedido (~1.7 kg/h, todo el H2 regenerado en un día) es algo menor que el de este '
                         'equipo (1.78-2.22 kg/h), pero pesa 1.1-1.8 t; el de referencia es el Hydro-Pac C06-15 (HPAC, '
                         '0.55 t)'],
        'hallazgos': ['equipo INDUSTRIAL de tierra (bastidor, refrigeración, controles); el motor de 16 kW es la '
                      'potencia instalada, no el consumo; no hay compresores de H2 aeronáuticos con datos publicados',
                      'la masa escala con el equipo, no con el flujo: no se escala hacia abajo; se usa entera como cota'],
        'decision': 'solo contraste (reemplazado como referencia por HPAC, 2026-09-27)', 'pedigri': 'C',
    },
    'HPAC': {
        'ambito': 'paper 1',
        'cita': 'Hydro-Pac Inc., LX-Series gas compressors for hydrogen, 1500 to 12000 psi (folleto técnico), '
                'Fairview, PA, 01/2011',
        'doi': 'no tiene (catálogo comercial)', 'arbitrada': False,
        'acceso': 'público (hydropac.com/literature/LX-SERIES Hydrogen Brochure 01_2011- reduced.pdf)',
        'archivo': REF + 'HydroPac_2011_LX-Series_Hydrogen_Brochure.pdf',
        'sha256': '534fe6f9b52d2d53c6bdeae56145874224ce256255b326787665e982fb2c0ed6',
        'formato': 'texto nativo con tablas desalineadas: capacidades y pesos confirmados contra la imagen (pp. 11, 13, '
                   '14)', 'erratas': 'no aplica (catálogo; "subject to change", p. 14)',
        'usos': ['EXP-16: compresor eléctrico de referencia (masa y tamaño) para el H2 regenerado'],
        'datos': {'COMPRESOR_REF': ({'P_bar': 414.0, 'flujo_kgph': (2.69, 2.94), 'motor_kW': 11.0,
                                     'masa_kg': (545.0, 545.0)},
                                    'modelo C06-15 (15 HP = 11 kW, 6000 psi máx.): 19.2 scfm a 300 psig de entrada '
                                    '(p. 13, serie 140/300) y 21 scfm a 450 psig (p. 14, serie 300/600), x 0.14 '
                                    'kg/h por scfm (p. 14, nota 8); peso con aceite 1200 lb = 545 kg (p. 11); '
                                    'confirmado contra la imagen 2026-09-27')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'C06-03 (3 HP) a 450 psig [scfm]': (4.4, 'p. 14'),
                            'C06-05 (5 HP), 260 kg, a 300 / 450 psig [scfm]': ((5.8, 7.1), 'pp. 11, 13 y 14: la '
                                                                             'referencia anterior (0.81-0.99 kg/h)'),
                            'C06-10 (10 HP), 410 kg, a 280 / 300 / 450 psig [scfm]': ((10.4, 11.1, 13.4),
                                                                                    'pp. 11, 13 y 14')},
        'consistencia': ['su conversión 0.14 kg/h por scfm = densidad del H2 a 60 F y 1 atm (0.1446 con el R del modelo)',
                         'lb y kg de su tabla coinciden (1200 lb = 544 kg; 570 lb = 259 kg)',
                         'potencia nominal / flujo del C06-15 = 3.7-4.1 kWh/kg: por encima del consumo típico (2-3), '
                         'como los motores sobredimensionados de GC25'],
        'hallazgos': ['equipo INDUSTRIAL (bastidor, reservorio de aceite de 148 l, refrigeración): su masa es cota '
                      'superior para un compresor de vuelo',
                      'folleto de 2011; capacidades a 60 Hz y gas a 27 C (p. 14)',
                      'con ~20 bar de entrada (275 psig, la presión de partida de la compresión en el modelo) el '
                      'caudal pedido de todo el H2 regenerado (~1.7 kg/h) no lo da el C06-10 (10.4-11.1 scfm a '
                      '280-300 psig = 1.46-1.55 kg/h); el C06-15 sí (18.0-19.2 scfm = 2.5-2.7 kg/h)',
                      '2026-09-27: antes se usaba el C06-05 (0.26 t) para la mitad del H2; el autor aclaró que se '
                      'comprime TODO al ritmo del electrolizador'],
        'decision': 'se usa como compresor de referencia el C06-15 (cota superior de masa, 0.55 t)', 'pedigri': 'C',
    },
    'MILT53029': {
        'ambito': 'paper 1',
        'cita': 'U.S. Department of Defense, MIL-T-53029C, Military specification: Tanks, fabric, collapsible: 3,000, '
                '10,000, 20,000, and 50,000 gallon, drinking water, 11 de mayo de 1994 (reemplaza a MIL-T-53029B, 1989)',
        'doi': 'no tiene (especificación militar; EverySpec MIL-T-53029C.051899)', 'arbitrada': False,
        'acceso': 'público ("Distribution Statement A. Approved for public release", p. 1); descargado de everyspec.com '
                  'el 2026-10-13; estado: inactiva (no se usa para diseños nuevos), lo que no cambia sus requisitos',
        'archivo': REF + 'MIL-T-53029C_1994_Tanks_Fabric_Collapsible_Drinking_Water.pdf',
        'sha256': '7c191984f6176557f5e02573b99c442bc93d7d9c8afa2420c81bd65d18a6fbd8',
        'formato': 'ESCANEO sin capa de texto (72 pp.): datos leídos de la imagen de las pp. 7 y 37, 2026-10-13',
        'erratas': 'tabla II, ítem 1: "30 in./62 maximum" se lee 30 mínimo / 62 máximo (errata tipográfica de "min.")',
        'usos': ['paper 1, D de PENDIENTES_PAPER1_2026-10-13: el tanque de agua (lastre) como caja fija por tonelada '
                 'liberada (balance_masa.tanque_agua); en el reparto, por tonelada entregada (EXP-20, EXP-22)'],
        'datos': {'capacidad_gal': (3000.0, 'p. 1, 1.2: el tamaño más chico de la norma'),
                  'lado_vacio_ft': (14.0, 'p. 37, fig. 1, vista superior: 14\'-0" (+0/-6") por lado, tanque vacío'),
                  'caras': (2, 'p. 37, fig. 1: almohada (dos caras de tela); llena, 12\'-6" x 12\'-6" x 4\'-0" de alto'),
                  'tela_ozyd2': ((30.0, 62.0), 'p. 7, tabla II, ítem 1 "Weight (oz/sq.yd)", tanques de 3000 y 10 000 '
                                               'gal (FED-STD-191, método 5041)')},
        'modelo': ('parametros', 'TANQUE_AGUA'),
        'consistencia': ['30-62 oz/yd2 = 1.02-2.10 kg/m2, en el rango de la tela de los tanques de agua de emergencia '
                         '(Oxfam: tela de PVC sobre poliéster de al menos 1100 g/m2)',
                         'la almohada llena (12.5 x 12.5 x 4 ft = 625 ft3 = 17.7 m3 de caja) contiene los 11.4 m3 de su '
                         'capacidad (una almohada no es una caja: 64 %)',
                         'tela de las dos caras: 36.4 m2 x 1.02-2.10 kg/m2 = 37-77 kg para 11.3 t de agua: 0.33-0.68 % de '
                         'la masa del agua (balance_masa.tanque_agua)'],
        'hallazgos': ['la norma no da la masa del tanque armado: solo la tela y la geometría; las válvulas de aluminio, '
                      'las asas, los parches y la lona de suelo no se cuentan (cota optimista, como las demás cajas)',
                      'contraste, equipo completo: la ficha del USMC del tanque de 3000 gal (TAMCN B21307B) da 617 lb de '
                      'envío y 25 653.2 lb operativo; operativo menos 3000 gal de agua = 617 lb: es el tanque con '
                      'accesorios, lona y caja, 2.5 % del agua; el tambor de 500 gal (B05717B), rodable, 7.6-8.8 %',
                      'contraste, dirigibles: el renglón "Ballast and Water System" de Goodyear 1975 (fuentes.GY75, '
                      'p. 20, ítem 8) pesa 950 lb (Los Angeles), 1229 (Graf Zeppelin), 3430 (Macon) y 3126 '
                      '(Hindenburg), 0.47-0.85 % de la sustentación bruta; sin la capacidad de agua de cada nave no se '
                      'pasa a kg por kg de agua (con la del Graf Zeppelin de airships.net, sin fuente primaria, ~5 %)'],
        'decision': 'se usa: la tela de la norma como la clase más liviana encontrada (0.33-0.68 % del agua); el '
                    'tanque del USMC y los sistemas de lastre de Goodyear, solo contraste', 'pedigri': 'B',
    },
    'SGT800': {
        'ambito': 'paper 1',
        'cita': 'Siemens AG, Gas turbine SGT-800 for power generation applications (ficha técnica), Article No. '
                'PGDG-T10040-00-7600, Power and Gas Division, Erlangen, 2017',
        'doi': 'no tiene (catálogo comercial)', 'arbitrada': False,
        'acceso': 'público (folleto de Siemens AG; descargado por el autor el 2026-10-11)',
        'archivo': 'referencias/equipos/Siemens_2017_SGT-800_datasheet_PGDG-T10040-00-7600.pdf',
        'sha256': '2ca4b94870e9bfea7a141919a0dbe02b064d3edc5d3695cf933078effb9d96a0',
        'formato': 'texto nativo con la tabla desalineada: columnas y pesos confirmados contra la imagen de la p. 2, '
                   '2026-10-11', 'erratas': 'no aplica (catálogo; "Subject to changes and errors", p. 2)',
        'usos': ['paper 1, sección 5.2 (crítica "compresión sin masa"): escala de la masa del compresor del escape a '
                 '14 atm; todo el paquete se atribuye a su compresor y se escala por el trabajo de compresión'],
        'datos': {'SGT800': ({'P_MWe': (47.5, 50.5, 54.0, 57.0), 'PR': (20.1, 21.0, 21.4, 21.8),
                              'm_escape_kgps': (132.8, 134.2, 135.5, 136.6),
                              'heat_rate_kJkWh': (9547.0, 9389.0, 9206.0, 8970.0), 'masa_kg': 285000.0},
                             'p. 2, "Simple cycle power generation": potencia, "Pressure ratio", "Exhaust gas flow" y '
                             '"Heat rate" de las cuatro columnas; "Physical dimensions", "Approx. weight" del paquete '
                             'clásico, 285 000 kg; nota: "For power generation, AC generator is included"; '
                             'condiciones ISO, gas natural')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'paquete Single Lift [kg]': (305000.0, 'p. 2: para barcazas, sobre patín; el clásico es '
                                                                  'el de generación en tierra'),
                            'temperatura de escape [C]': ((541, 553, 563, 565), 'p. 2'),
                            'compresor': ('15 etapas, álabes guía variables en las 3 primeras', 'p. 1')},
        'consistencia': ['eficiencia bruta = 3600 / heat rate en las cuatro columnas (37.7, 38.3, 39.1 y 40.1 %, como '
                         'las impresas)',
                         '285 000 kg = 628 300 lb, como en la tabla',
                         'el combustible (heat rate x potencia / PCI del metano, 2.5-2.8 kg/s) es ~2 % del escape: '
                         'aire de 130.3-133.8 kg/s'],
        'hallazgos': ['la masa es del paquete clásico entero (turbina, caja reductora, auxiliares, cerramiento y el '
                      'generador): atribuirla toda al compresor la sobreestima a esa escala',
                      'su compresor es 71-76 veces más potente que el del escape (trabajo isentrópico de una etapa) y '
                      'las máquinas chicas pesan más por kW: el texto lo declara con el doble de kg/kW',
                      'el aire de enfriamiento que se sangra del compresor no se descuenta (no se publica): su trabajo '
                      'real es algo menor y la masa estimada algo baja; lo cubre el factor 2',
                      'el caudal de la ficha es el de escape: el aire se obtiene restando el combustible'],
        'decision': 'se usa como escala de la masa del compresor del escape (cota del paquete entero), con el factor '
                    '2 por tamaño declarado', 'pedigri': 'C',
    },
    'ROTAX23': {
        'ambito': 'paper 1',
        'cita': 'BRP-Rotax GmbH & Co KG, Rotax 916 iS/c: product folder (folleto de producto), Gunskirchen, Austria, 2023',
        'doi': 'no tiene (catálogo comercial)', 'arbitrada': False,
        'acceso': 'público (flyrotax.com/assets/files/products/Rotax-916iSc-Product_Folder.pdf)',
        'archivo': REF + 'Rotax_2023_916iSc_product_folder.pdf',
        'sha256': '2667dbda82845a46cdcce3f12a733a20a567d09a0ea937c57b2946cc98d77db4',
        'formato': 'texto nativo; tabla de pesos y potencias confirmada contra la imagen (p. 2), 2026-09-27',
        'erratas': 'no aplica (catálogo; "for information only - without commitment")',
        'usos': ['paper 1, tabla 6: unidad de conversión, clase motor alternativo de H2 (extremo liviano, proxy de '
                 'aviación)'],
        'datos': {'MOTOR_ROTAX': (((117e3, 85.8), (101e3, 91.2)),
                                  'p. 2: 117 kW (160 hp) a 5800 1/min y 101 kW (137 hp) de potencia máxima continua a '
                                  '5500 1/min; motor con reductora i = 2.54, turbo, intercooler y escape '
                                  '85.8 kg; bastidor 2.0, alternador externo 3.0 y deflector 0.4 kg (91.2 kg en total)')},
        'modelo': ('parametros', None),
        'consistencia': ['160 hp = 117.7 kW y 137 hp = 102.2 kW: coinciden con los kW impresos (redondeo)',
                         '85.8 kg = 189 lb, como en la tabla'],
        'hallazgos': ['motor de GASOLINA: Bosch Aviation Technology lo convirtió a H2 (115 kW, Aviation Week, '
                      '31-03-2025) sin publicar su masa; se usa como PROXY de la clase de aviación',
                      'certificado EASA CS-E (iSc); TBO 2000 h'],
        'decision': 'se usa como extremo liviano (proxy de aviación) de la clase motor alternativo de H2',
        'pedigri': 'C',
    },
    'DEUTZ24': {
        'ambito': 'paper 1',
        'cita': 'DEUTZ AG, DEUTZ TCG 7.8 H2: product datasheet (ficha técnica, versión v1.0), Colonia, Alemania, 2024',
        'doi': 'no tiene (catálogo comercial)', 'arbitrada': False,
        'acceso': 'público (deutz.com.au/wp-content/uploads/2025/05/DEUTZ_ProductDatasheet_TCG_7.8_H2.pdf); el '
                  'consumo de H2 está en deutz.com/en/products/hydrogen-engines/',
        'archivo': REF + 'DEUTZ_2024_TCG_7.8_H2_datasheet.pdf',
        'sha256': '6e82fcc918a32e8d4431261453b3ad25176c393805767446619504310ea5a1bd',
        'formato': 'texto nativo con columnas desalineadas: potencia y peso confirmados contra la imagen (p. 2), '
                   '2026-09-27',
        'erratas': 'no aplica (catálogo)',
        'usos': ['paper 1, tabla 6: unidad de conversión, clase motor alternativo de H2 (extremo pesado)'],
        'datos': {'MOTOR_DEUTZ': (((220e3, 770.0),),
                                  'p. 2: power output acc. ISO 14396 220 kW a 2200 1/min; weight (DIN 70020 part 7A) '
                                  '~770 kg')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'consumo de H2 a plena potencia [kg/h]': (19, 'deutz.com: "~19 kg/h @max. power (~12 kg/h '
                                                                           'GenSet)"'),
                            'par máximo [N m]': (1000, 'p. 2, 1400-1600 1/min')},
        'consistencia': ['220 kW con ~19 kg/h de H2 (PCI 120 MJ/kg): eta = 0.35, dentro de 0.30-0.45 de '
                         'Verhelst & Wallner (parametros.eta_motor_H2)'],
        'hallazgos': ['motor industrial (tren, estacionario, maquinaria), no de aviación: extremo pesado de la clase',
                      'la ficha no dice lambda; el peso es "~" (aproximado)'],
        'decision': 'se usa como extremo pesado de la clase motor alternativo de H2 (único motor de H2 de serie con '
                    'ficha publicada)',
        'pedigri': 'C',
    },
    'HASK': {
        'ambito': 'paper 1',
        'cita': 'Haskel International, Gas boosters (catálogo de compresores neumáticos de gas), Burbank, CA',
        'doi': 'no tiene (catálogo comercial)', 'arbitrada': False, 'acceso': 'público (catálogo PDF del fabricante)',
        'archivo': REF + 'Haskel_gas_boosters_catalog.pdf',
        'sha256': 'c2c5ae2c10fd4de6d2d855eb3435be534db8a78120d2cfb477ef5a0159127ca7', 'formato': 'texto nativo',
        'erratas': 'no aplica (catálogo)',
        'usos': ['referencia: existen compresores de H2 de decenas de kg (no entra al balance)'],
        'datos': {'masa_booster_kg': ((17, 70), 'pp. 7-9: modelos AGT y 8/14AGT, 38 lb (17 kg) a 154 lb (70 kg)'),
                  'P_H2_bar': (345, 'pp. 7-8: salida máxima para hidrógeno de varios modelos AGT, 5000 psi')},
        'consistencia': ['libras y kilogramos del catálogo coinciden (38 lb = 17 kg; 121 lb = 55 kg)'],
        'hallazgos': ['son compresores NEUMÁTICOS: se mueven con aire comprimido, que exige otro compresor a bordo y '
                      'rinde poco; por eso no se usan en el balance, solo muestran que la masa puede ser mucho menor '
                      'que la cota industrial'],
        'decision': 'solo cualitativo (orden de magnitud de masa)', 'pedigri': 'C',
    },
    'PAS15': {
        'ambito': 'paper 1',
        'cita': 'Pasternak, I., Flight system for a constant volume, variable buoyancy air vehicle, U.S. Patent '
                '9,016,622 B1, 28 abr 2015 (solicitud 13/415,182, 8 mar 2012)',
        'doi': 'no tiene (patente de EE. UU.)', 'arbitrada': False,
        'acceso': 'público (Google Patents); examinada por la USPTO, no revisada por pares',
        'archivo': REF + 'Pasternak_2015_US9016622B1.pdf',
        'sha256': '217082d431d83bf184e2c158e1e2170fcf2856c6bf5b4157a91d3a520b8f6e76',
        'formato': 'ESCANEO con OCR (ImageMagick; flujos comprimidos dañados): se leyó la imagen de la p. 1',
        'erratas': 'no aplica (patente)',
        'usos': ['paper 1, sec. 7 (helio): que el gas de sustentación se puede comprimir a bordo'],
        'datos': {'concepto': ('comprimir y guardar el gas de sustentación en tanques a presión para controlar la '
                               'sustentación estática al cargar y descargar ("COSH")',
                               'p. 1, resumen (57); descripción')},
        'consistencia': ['el resumen (p. 1) y la descripción dicen lo mismo: gas de sustentación comprimido a '
                         'tanques, volumen casi constante, sin balonets'],
        'hallazgos': ['patente: describe el concepto, sin mediciones ni masas; el paper solo la usa para decir que '
                      'el helio se puede comprimir a bordo'],
        'decision': 'solo cualitativo (existencia del concepto)', 'pedigri': 'C',
    },

    # ── Datos de las naves de la Tabla 2 (hallazgo S3, autor 2026-09-28) ──
    # parametros.PLATAFORMAS toma estos datos; el grupo W de verificacion.py comprueba que sean los mismos.
    'GY75': {
        'ambito': 'paper 1',
        'cita': 'Goodyear Aerospace Corporation, Feasibility study of modern airships (Phase I), Vol. III: '
                'Historical overview (Task I), NASA CR-137692 (Vol. 3), NTIS N75-32034, Akron, OH, agosto de 1975 '
                '(la carátula del NTIS dice "Phase II"; la portada interior y el NTRS, Phase I)',
        'doi': 'no tiene (NTRS 19750023961)', 'arbitrada': False,
        'acceso': 'abierto (NTRS); informe de contratista de la NASA, sin revisión por pares',
        'archivo': REF + 'Goodyear_1975_NASA_Modern_Airships_PhII_Vol3_Historical.pdf',
        'sha256': '92776cb5d97284b25e1e961a3a107ae297b4d28a4c9d3da5b4a70e26c698681b',
        'formato': 'ESCANEO sin capa de texto: tablas de las pp. 19-20 transcritas por el autor y confirmadas '
                   'contra la imagen de la página',
        'erratas': 'del informe: encabezado "LZ-127 Graf Zeppelin (1298)" por 1928; desglose del casco del '
                   'Hindenburg (a-f) suma 110 166 lb, impreso 110 435; Macon vacío 236 497 (p. 20) vs 236 493 (p. 19)',
        'usos': ['paper 1, Tabla 2: volumen total, peso vacío, largo y diámetro del LZ 129 y del ZR-3 (H2)',
                 'paper 1, ruta nueva: recuperación de agua de los rígidos como candidato de la tabla de brecha '
                 '(req_masa_rechazo.py)'],
        'datos': {'LZ129_V_total_ft3': (7_650_000, 'p. 19, Total Volume'),
                  'LZ129_vacio_lb': (249_000, 'p. 19, Weight Empty (= p. 20, total)'),
                  'LZ129_L_ft': (814, 'p. 19, Length'), 'LZ129_D_ft': (135.2, 'p. 19, Diameter'),
                  'ZR3_V_total_ft3': (2_800_000, 'p. 19, Total Volume'),
                  'ZR3_vacio_lb': (86_400, 'p. 19, Weight Empty, versión H2 (He: 89 239)'),
                  'ZR3_L_ft': (660, 'p. 19, Length'), 'ZR3_D_ft': (90.7, 'p. 19, Diameter'),
                  # equipos de recuperación de agua del escape (p. 20, ítem 18) y sustentación bruta (p. 20)
                  'agua_lb': ({'LA_He': 4_000, 'LZ127': 5_876, 'Macon': 12_922, 'LZ129': 8_336},
                              'p. 20, item 18 "Water Recovery"'),
                  'bruto_lb': ({'LA_He': 153_140, 'LZ127': 260_300, 'Macon': 403_465, 'LZ129': 456_076},
                               'p. 20, "Gross Lift"'),
                  # consumo y velocidad de las mismas cuatro naves (ruta nueva, tabla de brecha; autor 2026-09-28)
                  'consumo_lbh': ({'LA_He': 450, 'LZ127': 490, 'Macon': 780, 'LZ129': 660},
                                  'p. 19, "Fuel Cons *** - lb/hr" (nota ***: a 50 nudos; LA: columna He)'),
                  'vmax_kn': ({'LA_He': 65.2, 'LZ127': 69, 'Macon': 72.2, 'LZ129': 72},
                              'p. 19, "Max. Speed - knots" (LA: columna He)')},
        'consistencia': ['bruto = vacío + útil EXACTO en las 10 naves de la p. 19 (ZR-3 He: 1 lb de redondeo)',
                         'L/D recalculado = la columna del informe en las 10 naves',
                         'gas / total = 0.923-0.939 en las 10 naves',
                         'p. 20: suma de ítems = peso vacío en las 5 columnas; % del bruto coincide'],
        'hallazgos': ['Hindenburg: gas 7 060 000 ft3 = informe de 1937 (7 063 000); largo 814 ft vs 803.8 ft y '
                      'bruto 206.9 t vs 214.5 t del informe de 1937 (otra forma de medir y otras condiciones)',
                      'ZR-3: gas 73 624 m3 vs 70 000 m3 de Dürr 1924 (DUR24); se usa el volumen TOTAL de esta fuente',
                      'ítem 18 de la p. 20 ("water recovery"): equipos de recuperación de agua del escape, '
                      '1.8-5.9 t (antecedente del ciclo)'],
        'decision': 'se usa (volumen total, peso vacío y medidas del LZ 129 y del ZR-3)', 'pedigri': 'B',
    },
    'DOC37': {
        'ambito': 'paper 1',
        'cita': 'R. W. Knight, The Hindenburg accident: a comparative digest of the investigations and findings, with '
                'the American and translated German reports included, Report No. 11, Safety and Planning Division, '
                'Bureau of Air Commerce, U.S. Department of Commerce, Washington, DC, agosto de 1938 (portada; ROSA P '
                'dot:81236 lo cataloga como serie CAA Safety and Planning Reports)',
        'doi': '10.21949/1511799 (comprobado 2026-09-28: resuelve a rosap.ntl.bts.gov/view/dot/81236, el mismo archivo)',
        'arbitrada': False,
        'acceso': 'abierto (repositorio ROSA P del DOT); informe gubernamental oficial',
        'archivo': REF + 'CAA_1938_Hindenburg_accident_dot81236.pdf',
        'sha256': '7eba1d6c1a0e177d937d3329e350bb92021293a56ab7506a6d0e8305878ce8a8',
        'formato': 'PDF del DOT; datos leídos por el autor (páginas por anotar)',
        'erratas': 'el OCR de otra copia daba 7,043,000 ft3; el autor lee 7,063,000',
        'usos': ['paper 1, Tabla 2: carga máxima del LZ 129 ("total payload")'],
        'datos': {'LZ129_carga_lb': (41_990, '"total payload"'),
                  'LZ129_equipado_lb': (430_950, '"weight of ship with necessary equipment and fuel"'),
                  'LZ129_sustentacion_lb': (472_940, '"total lift (under standard conditions)"')},
        'consistencia': ['430 905 + 41 990 = 472 895 lb: 45 lb (0.01 %) menos que la sustentación impresa; '
                         'descuadre del propio informe, sin efecto'],
        'hallazgos': ['el peso del modelo viejo (263.6 t) superaba esta sustentación total (214.5 t)'],
        'decision': 'se usa (carga máxima del LZ 129)', 'pedigri': 'B',
    },
    'DUR24': {
        'ambito': 'paper 1',
        'cita': 'Dürr, L., The American airship ZR-3, NACA Technical Memorandum No. 286, Washington, DC, '
                'noviembre de 1924',
        'doi': 'no tiene (NTRS 19930086540)', 'arbitrada': False,
        'acceso': 'abierto (NTRS); Dürr fue el diseñador jefe de Zeppelin',
        'archivo': REF + 'Durr_1924_NACA_TM286_ZR3.pdf',
        'sha256': '8e99ce4121687a0de29dd4f19f67324a0f8495ca08616fa5c6fb38b77c74e7fb',
        'formato': 'ESCANEO con OCR (NTRS)',
        'erratas': 'un resumen automático leyó "6,000 kg" donde el texto dice 46,000 kg (101 412 lb)',
        'usos': ['paper 1, Tabla 2: carga comercial del ZR-3'],
        'datos': {'ZR3_carga_kg': (5_000, '"pay load of 5 metric tons (11,023 lb.) of passengers, baggage and '
                                          'freight"'),
                  'ZR3_capacidad_kg': (46_000, 'capacidad de carga total con H2 de gravedad específica 0.1, '
                                               '760 mmHg, 0 C y 60 % de humedad')},
        'consistencia': ['2 472 000 ft3 = 70 000 m3; 656.2 ft = 200 m; 90.68 ft = 27.64 m',
                         'sustentación con sus condiciones: 70 000 x (1.2905 - 0.129) = 81.3 t; menos 46 t = 35.3 t de '
                         'peso muerto (Goodyear 1975 da 39.2 t vacío, en otras condiciones)'],
        'hallazgos': ['gas 70 000 m3 vs 73 624 m3 de GY75: se usa el volumen total de GY75'],
        'decision': 'se usa (carga comercial del ZR-3)', 'pedigri': 'B',
    },
    'EASANT': {
        'ambito': 'paper 1',
        'cita': 'EASA, Type-Certificate Data Sheet No. EASA.AS.001, Zeppelin LZ N07-100 / LZ N07-101, Issue 10',
        'doi': 'no tiene (documento de certificación)', 'arbitrada': False,
        'acceso': 'abierto (https://www.easa.europa.eu/en/document-library/type-certificates/airship-cs-30-cs-31ha/'
                  'easaas001-zeppelin-lz-n07-100, consultado 2026-09-28); dato de certificación',
        'archivo': REF + 'EASA_TCDS_AS.001_ZeppelinNT_Issue10_track_changes.pdf',
        'sha256': 'f99db5453f97ff948d4bd921c5cdda6f77222e8a7591e335c9849d4a88cde923',
        'formato': 'PDF con control de cambios (para citar, bajar la versión limpia de la ed. 10)',
        'erratas': 'ninguna conocida',
        'usos': ['paper 1, Tabla 2: volumen de la envolvente y medidas del Zeppelin NT'],
        'datos': {'NT_V_envolvente_m3': (8_425, 'pp. 5 (N07-100) y 10 (N07-101)'),
                  'NT_L_m': (75.1, 'pp. 5 y 10'), 'NT_D_m': (14.2, 'pp. 5 y 10')},
        'consistencia': ['los dos modelos (-100 y -101) dan la misma envolvente y medidas'],
        'hallazgos': ['no publica peso vacío (la "max. gondola mass", 2 690-3 100 kg, es solo la góndola) y la carga '
                      'máxima figura como "reserved" (sec. 20): el peso vacío del modelo se deriva y la carga es la '
                      'liberación de referencia, ambos declarados'],
        'decision': 'se usa (volumen y medidas)', 'pedigri': 'A',
    },
    'FWH26': {
        'ambito': 'paper 1',
        'cita': 'Flying Whales, LCA60T (https://www.flying-whales.com/en/the-lca60t/) y preguntas frecuentes '
                '(https://www.flying-whales.com/en/faq/), consultados 2026-09-28',
        'doi': 'no tiene (sitio del fabricante)', 'arbitrada': False,
        'acceso': 'abierto (web del fabricante); apoyo de medidas y carga en cc0997en.pdf, p. 3',
        'archivo': REF + 'cc0997en_LCA60T.pdf',
        'sha256': '4df467705a1cccf1149783e5c378e0fab939933c192519a5a3a63d362132644b',
        'formato': 'páginas web (texto); el PDF de apoyo es texto nativo',
        'erratas': 'la FAQ dice 200 000 m3 de helio en 14 celdas y la página del LCA60T 180 000 m3',
        'usos': ['paper 1, Tabla 2: caso base (volumen, peso vacío, medidas, carga)'],
        'datos': {'FW_V_helio_m3': (180_000, 'página del LCA60T, "Volume of helium"'),
                  'FW_vacio_t': (100, 'FAQ, peso de la estructura y los sistemas'),
                  'FW_L_m': (200, 'FAQ; cc0997en p. 3'), 'FW_D_m': (50, 'FAQ; cc0997en p. 3'),
                  'FW_carga_t': (60, 'FAQ y página del LCA60T; cc0997en p. 3')},
        'consistencia': ['100 t + 60 t = 160 t <= 180 000 m3 x 1.056 kg/m3 (helio) = 190 t'],
        'hallazgos': ['solo del fabricante; discrepancia 180 000 / 200 000 m3: el autor eligió 180 000 (el lado '
                      'exigente: casco menor, más flujo) y se declara la otra cifra (el flujo cambia ~5 %)'],
        'decision': 'se usa (caso base), con la discrepancia declarada', 'pedigri': 'C',
    },
    # ── Ruta nueva del paper 1 (autor, 2026-09-28): requisito de masa del rechazo + brecha con lo que existe ──
    'ZSS20': {
        'ambito': 'paper 1',
        'cita': 'Zhang, L., Spatari, S., Sun, Y., Life cycle assessment of novel heat exchanger for dry cooling of '
                'power plants based on encapsulated phase change materials, Appl. Energy 271 (2020) 115227',
        'doi': '10.1016/j.apenergy.2020.115227', 'arbitrada': True,
        'acceso': 'suscripción (ScienceDirect); PDF del autor. Manuscrito abierto en OSTI 1799023',
        'archivo': REF + 'Zhang_Spatari_Sun_2020_ApplEnergy_271_115227.pdf',
        'sha256': 'd1b4b74877e6846046357f9e2faa2315ff1b967e3654c63d4dc90d9f0624dc06',
        'formato': 'texto nativo; la tabla 3 sale desalineada en el texto extraído: confirmada contra la imagen '
                   '(p. 5) el 2026-09-29',
        'erratas': 'Crossref 2026-09-29: no consultado [COMPLETAR]',
        'usos': ['paper 1, ruta nueva: condensador enfriado por aire de central (ACC) como candidato de la tabla de '
                 'brecha (req_masa_rechazo.py)'],
        'datos': {'ACC_Q_MWth': (388.9, 'tabla 3, "Cooling duty, MWth", verano'),
                  'ACC_vent_MWe': (4.01, 'tabla 3, "Fan power, MWe", verano'),
                  'ACC_T_aire_C': (33.33, 'p. 4, tabla del clima del sitio (Chicago), "Dry bulb temperature", verano'),
                  'ACC_ITD_C': (23.6, 'tabla 3, "ITD", verano')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'ACC promedio / invierno: calor [MWth], ventiladores [MWe], ITD [C]':
                                ((375.8, 4.04, 27.5), (378, 3, 40.1), 'tabla 3'),
                            'planta: potencia neta [MWe], eficiencia (verano)': ((646, 0.5085), 'tabla 3')},
        'consistencia': ['COP de la tabla 3 = calor / ventiladores: 388.9 / 4.01 = 97.0 (impreso 97)',
                         'la potencia neta y el ITD de la tabla 3 = los de la tabla 2 (verano, 646 MWe, 23.6 C)'],
        'hallazgos': ['la planta es un ciclo combinado de gas natural de ~688 MWe; el ACC solo condensa el vapor de la '
                      'turbina de vapor (p. 12)',
                      'diseño para tierra: acero al carbono y ventiladores al ~1 % del calor (la electricidad es cara, '
                      'la masa no importa); para una nave la compensación sería la opuesta (más aire, menos área)'],
        'decision': 'se usa como referencia de "la física del paper, hecha como se hace hoy en tierra" (no es un '
                    'diseño para la nave)', 'pedigri': 'B',
    },
    'IBX22': {
        'ambito': 'paper 1',
        'cita': 'Ionblox Inc., Launch (High Energy AV Cell), product sheet AV-Product-Sheet_IB-2022, Fremont, CA, 2022',
        'doi': 'no tiene (ficha del fabricante)', 'arbitrada': False,
        'acceso': 'abierto (ionblox.com/air, "Download datasheet")',
        'archivo': REF + 'Ionblox_2022_Launch_AV_cell_product_sheet.pdf',
        'sha256': '807aa867827db3be941454022efa976e6e31bf50622fd3300aada08bfa22032a',
        'formato': 'texto nativo (1 página)', 'erratas': 'no aplica (ficha comercial)',
        'usos': ['contraste de la búsqueda del ítem 8 (la sensibilidad usa AMP26, mayor potencia a 30 s)'],
        'datos': {'kWkg_celda_30s': (3.8, '"Specific discharge power at 50% SOC from a 12C, 30s pulse": 3800 W/kg'),
                  'Whkg_celda': (340.0, '"Specific energy at C/3 rate": 340 Wh/kg')},
        'datos_no_usados': {'capacidad [Ah], masa [g], voltaje nominal [V]': ((12, 121, 3.45), 'misma ficha'),
                            'ciclos (100 % DOD, 1C-1C)': ('> 700', 'misma ficha')},
        'consistencia': ['12C x 12 Ah = 144 A x ~3.2 V en descarga / 0.121 kg = ~3.8 kW/kg: la potencia sale de su '
                         'capacidad, su masa y la tasa',
                         'búsqueda del ítem 8 (2026-09-29): paquete comercial medido con celdas Molicel P45B '
                         '(Power Battery, ficha 2024) 1.1-1.6 kW/kg, dentro del rango del modelo; Molicel P45B 2.4 '
                         'kW/kg de celda (10 s); Toshiba SCiB LTO ~3.5 kW/kg de celda (10 s)'],
        'hallazgos': ['dato del fabricante, sin ensayo independiente; ciclo de vida corto (> 700 ciclos)',
                      'es de CELDA: a paquete con el 0.66 del X-57 (CHIN20) da ~2.5 kW/kg'],
        'decision': 'contraste (2026-09-29): superada por AMP26 en la sensibilidad; no se cita en el paper',
        'pedigri': 'C',
    },
    'AMP26': {
        'ambito': 'paper 1',
        'cita': 'Amprius Technologies, SiCore SA88 10 Ah ultra-high power cell, technical data sheet AMP-TDS-SA88-R3, '
                'en: Product Portfolio Q3 2026, Fremont, CA, 2026',
        'doi': 'no tiene (catálogo del fabricante)', 'arbitrada': False,
        'acceso': 'abierto (amprius.com/documents/Amprius_Product_Catalog.pdf)',
        'archivo': REF + 'Amprius_2026_Product_Portfolio_Q3.pdf',
        'sha256': 'e5a1dc8208d5a099f0d5b1ab1c884f079b94a3a1c55597de914b4fa418a129c2',
        'formato': 'texto nativo; la ficha del SA88 (p. 26) confirmada contra la imagen el 2026-09-29',
        'erratas': 'no aplica (catálogo; "subject to change without notice")',
        'usos': ['paper 1, ítem 8: sensibilidad del requisito del rechazo a la mejor batería de hoy (req_masa_rechazo.py)'],
        'datos': {'AMP_SA88': ({'I_pulso_30s_A': 153.0, 'V_nom_V': 3.45, 'masa_g': 96.5, 'Whkg': 365.0},
                               'p. 26: "Max Pulse (<= 30 seconds) 153 A (15C)"; "Nominal 3.45 V"; "Weight 96.5 +- 2g"; '
                               '"Gravimetric (@ 0.2C Discharge) 365 Wh/kg"')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'corriente continua máxima [A]; capacidad [Ah]; ciclos a 80 % (+1C/-3C)':
                                ((102, 10.2, 200), 'p. 26')},
        'consistencia': ['energía de la ficha 35.19 Wh / 96.5 g = 365 Wh/kg, la impresa',
                         '153 A x 3.45 V / 96.5 g = 5.5 kW/kg de celda con el voltaje nominal (cota superior: a 15C el '
                         'voltaje cae); 30 s a esa potencia piden 46 Wh/kg, menos que sus 365'],
        'hallazgos': ['dato del fabricante, sin ensayo independiente; 200 ciclos; pide apriete externo de 30 psi',
                      'es de CELDA: a paquete con el 0.66 del X-57 (CHIN20) da ~3.6 kW/kg; con el mejor empaquetado de '
                      'LOB20 (0.742), ~4.1'],
        'decision': 'solo sensibilidad (la clase representativa sigue siendo la de Pattanayak 2025); autor, 2026-09-29',
        'pedigri': 'C',
    },
    'LOB20': {
        'ambito': 'paper 1',
        'cita': 'Löbberding, H., Wessel, S., Offermanns, C., Kehrer, M., Rother, J., Heimes, H., Kampker, A., From cell '
                'to battery system in BEVs: analysis of system packing efficiency and cell types, World Electr. Veh. J. '
                '11 (4) (2020) 77',
        'doi': '10.3390/wevj11040077', 'arbitrada': True, 'acceso': 'abierto (MDPI)',
        'archivo': REF + 'Lobberding_2020_WEVJ_11_77.pdf',
        'sha256': '3deff29895904174efba68beb2a5619abb9033b6eca25348cd9dfd46686088ec',
        'formato': 'texto nativo; la tabla 1 sale desalineada en el texto extraído: confirmada contra la imagen (p. 8)',
        'erratas': 'Crossref 2026-09-29: no consultado [COMPLETAR]',
        'usos': ['contraste del factor celda -> paquete (0.66 del X-57, CHIN20) en la sensibilidad de la batería'],
        'datos': {'f_celda_paquete_bev': ((0.484, 0.742), 'tabla 1, "Efficiency rate of cell level to system level", '
                                          'gravimétrico: el menor (cilíndricas, 0.484) y el mayor (prismáticas y bolsa, '
                                          '0.742) de 25 autos eléctricos 2010-2019')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'promedios gravimétricos (prismáticas y bolsa, cilíndricas)': ((0.575, 0.552), 'tabla 1')},
        'consistencia': ['el 0.66 del X-57 (CHIN20) cae dentro del rango y sobre los promedios (0.55-0.58)'],
        'hallazgos': ['paquetes de AUTOS (energía), no de alta potencia ni de aviación; datos de fuentes externas '
                      '(p. 12, los autores lo dicen)'],
        'decision': 'contraste: el 0.66 se queda; con el mejor empaquetado (0.742) el requisito sigue sobre el '
                    'sistema de FlyZero 2035 (verificación, grupo X)', 'pedigri': 'B',
    },
    'ZSS20S': {
        'ambito': 'paper 1',
        'cita': 'Zhang, L., Spatari, S., Sun, Y., Appl. Energy 271 (2020) 115227, material suplementario (Appendix A, '
                'mmc1): tablas S1-S18',
        'doi': '10.1016/j.apenergy.2020.115227 (suplementario)', 'arbitrada': True,
        'acceso': 'libre desde la página del artículo (ScienceDirect, Appendix A)',
        'archivo': REF + 'Zhang_Spatari_Sun_2020_ApplEnergy_271_115227_suplementario.pdf',
        'sha256': 'aa544e0d6e73a5513723151ce0b272e395212fd79d023ec1dc33cf374801f488',
        'formato': 'texto nativo; la tabla S7 sale desalineada en el texto extraído: confirmada contra la imagen '
                   '(p. S8) el 2026-09-29', 'erratas': 'ver ZSS20',
        'usos': ['paper 1, ruta nueva: masa del ACC (req_masa_rechazo.py)'],
        'datos': {'ACC_masas_t': ({'haces_acero': 550, 'haces_aluminio': 46, 'estructura': 1300, 'ductos': 650,
                                   'tuberias_acero': 150, 'tuberias_inox': 3, 'bomba_aire': 5, 'mecanico_acero': 45,
                                   'mecanico_prfv': 65, 'tanque': 25, 'electrico_cobre': 5, 'electrico_acero': 30,
                                   'electrico_plastico': 5},
                                  'tabla S7, "Mass (ton)": haces (tubo de acero al carbono, aleta de aluminio), '
                                  'estructura de soporte, ductos y carcasa, tuberías y válvulas, bomba de aire, '
                                  'equipo mecánico (acero; ventilador y campana de PRFV), tanque de condensado, '
                                  'eléctrico y control')},
        'modelo': ('parametros', None),
        'datos_no_usados': {'costos de capital por equipo [$]': ('4 544 320 haces ... 1 468 000 eléctrico',
                                                                 'tabla S7')},
        'consistencia': ['la suma de la tabla S7 es 2879 ton; núcleo + estructura = (596 + 1300) / 2879 = 66 %, el '
                         'porcentaje que da el paper (p. 9)'],
        'hallazgos': ['masas ESTIMADAS a partir de los parámetros de diseño (p. 7, "estimated based on the cooling '
                      'system design parameters"), no pesadas',
                      '"ton" sin aclarar (el paper es de EE. UU.: puede ser tonelada corta, 0.907 t; cambia 10 %)'],
        'decision': 'se usa (masa del ACC de referencia), con las dos salvedades declaradas', 'pedigri': 'C',
    },
    'HOW99': {
        'ambito': 'paper 1',
        'cita': 'Howe, D., Ground handling and mooring, cap. 11 de Khoury, G. A., Gillett, J. D. (Eds.), Airship '
                'Technology, 1.a ed., Cambridge University Press, 1999, pp. 297-320',
        'doi': 'no tiene (libro)', 'arbitrada': False,
        'acceso': 'libro editado (CUP); copia digitalizada por Internet Archive (2022)',
        'archivo': REF + 'Khoury_Gillett_1999_Airship_Technology_1ed.pdf',
        'sha256': 'd6dee14bb2a5fa7c7d4e916ca97c112fa9810db6d4b3cfa59149667a09665ec2',
        'formato': 'ESCANEO con OCR (Internet Archive); el pasaje se leyó en el texto extraído (pdftotext) y se '
                   'ubicó por los encabezados de las páginas vecinas (296 y 298)',
        'erratas': 'del índice: da el cap. 11 en la p. 296 y el 12 en la 319; las páginas impresas son 297 y 321 '
                   '(se cita 297-320)',
        'usos': ['paper 1, discusión, "Safety": ventear H2 al descender era peligroso cerca del suelo'],
        'datos': {'venteo_al_descender': ('"Lifting Gas Management": los dirigibles de H2 eran propensos al incendio '
                                          'catastrófico junto al suelo; el aterrizaje era un momento especialmente '
                                          'peligroso por el venteo frecuente de H2 para bajar la sustentación al '
                                          'descender, junto con la descarga estática al tocar tierra; el Graf Zeppelin '
                                          'reducía el venteo con un combustible de la densidad del aire',
                                          'p. 297, primer apartado del cap. 11')},
        'consistencia': ['el cap. 11 empieza en la p. 297 (pdf 313) y el 12 en la 321 (pdf 337); la p. 320 es una '
                         'foto'],
        'hallazgos': ['es una afirmación cualitativa de un capítulo de libro, sin datos de accidentes'],
        'decision': 'se usa como antecedente cualitativo (seguridad)', 'pedigri': 'C',
    },
    'BA19': {
        'ambito': 'paper 1',
        'cita': 'Ben Abdallah, F., Azouz, N., Beji, L., Abichou, A., Modeling of a heavy-lift airship carrying a payload '
                'by a cable-driven parallel manipulator, Int. J. Adv. Robot. Syst. 16 (2019) 1729881419861769',
        'doi': '10.1177/1729881419861769', 'arbitrada': True,
        'acceso': 'abierto (CC BY, SAGE); PDF del editor descargado por el autor el 2026-10-15',
        'archivo': REF + 'BenAbdallah_2019_IJARS_16_1729881419861769.pdf',
        'sha256': 'c9a27414d95f6e98d8113f6792cca2e8733bf609a9393f9d5c985c387753c54f',
        'formato': 'texto nativo',
        'erratas': 'no se buscaron (cita cualitativa)',
        'usos': ['paper 1, introducción: cargar y descargar en vuelo estacionario es una fase peligrosa'],
        'datos': {'carga_en_vuelo_peligrosa': ('"One hazardous task in an airship mission is the loading and unloading '
                                               'of the heavy load by means of motorized winches while the airship '
                                               'hovers above the ground"',
                                               'p. 12 del PDF, inicio de "Simulation results"')},
        'consistencia': ['el resumen no lo dice; lo dice el cuerpo (p. 12), y la introducción habla de la seguridad del '
                         'conjunto por la oscilación de la carga suspendida'],
        'hallazgos': ['no dice que el exceso de sustentación desplace el dirigible "by meters within seconds" ni '
                      'menciona personal en tierra; la oración del paper que se lo atribuía se corrigió el 2026-10-15'],
        'decision': 'se usa como antecedente cualitativo (operación de carga)', 'pedigri': 'C',
    },
    'BAR13': {
        'ambito': 'paper 1',
        'cita': 'Barton, S. A., Airship buoyancy control using inflatable vacuum chambers, J. Aircr. 50 (2013) 1001-1010',
        'doi': '10.2514/1.C031654', 'arbitrada': True,
        'acceso': 'de pago (AIAA); copia no oficial (ilide.info) descargada por el autor el 2026-10-15, maquetada como '
                  'manuscrito aceptado ("Vol. , No. ,"), recibido 6 sep 2011, revisión 5 ene 2013',
        'archivo': REF + 'Barton_2013_JAircraft_50_1001.pdf',
        'sha256': '4e96bf4f36110a26434cd4f7dfe87323eda597b4af60e2c1f5fc933ffe83d3ac',
        'formato': 'texto nativo',
        'erratas': 'no se buscaron (cita cualitativa)',
        'usos': ['paper 1, introducción y sección 7: las cámaras de vacío convencionales son demasiado pesadas'],
        'datos': {'vacio_convencional_pesado': ('"conventional vacuum chambers are impractical for this purpose because '
                                                'of their large mass (high weight)"', 'p. 1, introducción')},
        'consistencia': ['propone la cámara de vacío inflable (IVC) como alternativa liviana: diseño conceptual '
                         '(tubo de ~1 m de radio y ~600 m, ~7.5 t de rango), estable en su análisis elástico; solo se '
                         'fabricaron muestras de tela; en sus conclusiones trabajan en una IVC de ~50 m3. Concluye '
                         'que la IVC pesa lo mismo que un sistema que solo comprime gas; su ventaja es agregar peso '
                         'rápido'],
        'hallazgos': ['no dice que la masa "has prevented flight" ni que toda cámara de vacío sea "prohibitive": el '
                      'paper se corrigió el 2026-10-15 para citar solo las convencionales; las inflables no se '
                      'mencionan (sin prototipo funcional conocido después de 2013, 8 citas en OpenAlex)'],
        'decision': 'se usa como antecedente cualitativo (cámaras de vacío)', 'pedigri': 'C',
    },
    'ZR22': {
        'ambito': 'paper 1',
        'cita': 'Zhou, S. Y., Renaud, M. (Flying Whales), Airship in wood harvesting: the case of the FLYING WHALES\' '
                'LCA60T, FAO, Roma, 2022 (documento cc0997en)',
        'doi': 'no tiene (repositorio de FAO, https://www.fao.org/3/cc0997en/cc0997en.pdf)', 'arbitrada': False,
        'acceso': 'abierto (FAO); la copia local es idéntica a la de FAO (SHA-256 comprobado el 2026-10-15)',
        'archivo': REF + 'cc0997en_LCA60T.pdf',
        'sha256': '4df467705a1cccf1149783e5c378e0fab939933c192519a5a3a63d362132644b',
        'formato': 'texto nativo',
        'erratas': 'no se buscaron (datos del fabricante)',
        'usos': ['paper 1, introducción, sección 3.1 (tabla 2: carga y medidas del LCA60T) y sección 5.2 (carga en '
                 'vuelo estacionario)'],
        'datos': {'LCA60T': ('60 t de carga, rígido de 200 m x 50 m, descarga en vuelo estacionario en zonas sin '
                             'infraestructura', 'p. 1 (resumen) y p. 3')},
        'consistencia': ['las medidas y la carga coinciden con la web del fabricante (FWH26); no da el volumen ni la '
                         'masa en vacío, que siguen citándose de la web'],
        'hallazgos': ['documento de ingenieros del fabricante publicado por FAO, no revisado por pares'],
        'decision': 'se usa para la carga, las medidas y la operación en vuelo estacionario', 'pedigri': 'C',
    },
    'FWEP24': {
        'ambito': 'paper 1',
        'cita': 'Bougon, S., Le François, B. (Flying Whales SAS), Device and method for ballast management, patente '
                'europea EP 3980329 B1, concedida el 17 ene 2024 (solicitud WO 2020/245527 A1)',
        'doi': 'no tiene (patente)', 'arbitrada': False,
        'acceso': 'abierto (Google Patents); PDF de la patente concedida, en francés',
        'archivo': REF + 'FlyingWhales_EP3980329B1_ballast_management.pdf',
        'sha256': '1c0e2ddbca0cf782563d65a9df854efb31b482516d432d6cd9e970c99714686d',
        'formato': 'texto nativo (francés)',
        'erratas': 'no aplica (patente)',
        'usos': ['paper 1, sección 5.2 y discusión: intercambio de la carga por lastre de agua desde un tanque en tierra'],
        'datos': {'lastre_desde_tierra': ('"une pompe 12 raccordée à une citerne (non représentée) présente au sol sur le '
                                          'site de dépose" (el agua del suelo lastra el bastidor de carga al descargar)',
                                          'párrafo [0044]')},
        'consistencia': ['la patente no nombra al LCA60T; la web del fabricante (FAQ) describe para el LCA60T el '
                         'mismo intercambio con unas 60 t de agua de tanques en tierra'],
        'hallazgos': ['método patentado, no una operación demostrada'],
        'decision': 'se usa como antecedente del intercambio de carga por lastre de agua', 'pedigri': 'C',
    },
    'FWUS19': {
        'ambito': 'paper 1',
        'cita': 'Kuhlmann, H. F. (Flying Whales SAS), Lighter-than-air aircraft and deballasting method implemented in '
                'said lighter-than-air aircraft, patente de EE. UU. 10,494,076 B2, 3 dic 2019 (también CA 2995961 C, '
                'EP 3183169 B1)',
        'doi': 'no tiene (patente)', 'arbitrada': False,
        'acceso': 'abierto (Google Patents); PDF de la patente concedida',
        'archivo': REF + 'FlyingWhales_US10494076B2_deballasting.pdf',
        'sha256': 'a2380530ca2a1884aefd3fa94c8bc6ab51af64962d567ce3c8b64ee3b1e40220',
        'formato': 'texto nativo',
        'erratas': 'no se buscaron (ejemplo de realización)',
        'usos': ['paper 1, sección 3.5 (Airship data): el volumen de helio publicado del LCA60T, 180 000 m3, es también '
                 'el volumen desplazado del dirigible de 60 t de esta patente'],
        'datos': {'volumen_desplazado_m3': (180_000, 'descripción, ejemplo: "a dirigible having a load capacity of sixty '
                                                     'tonnes ... a length of 150 metres and a diameter of 45 metres, for '
                                                     'a displacement volume of 180,000 m3"; helio "approximately '
                                                     '110,000 m3"')},
        'consistencia': ['es una configuración anterior y más compacta (150 x 45 m) que la del LCA60T actual (200 x 50 m)'
                         '; por eso 180 000 m3 se toma como cota inferior del casco actual'],
        'hallazgos': ['ejemplo no limitante de una patente, no una ficha del LCA60T'],
        'decision': 'se usa para justificar el volumen del casco del LCA60T como cota inferior', 'pedigri': 'C',
    },
}


def sha256(ruta):
    h = hashlib.sha256()
    with open(ruta, 'rb') as fh:
        for bloque in iter(lambda: fh.read(1 << 20), b''):
            h.update(bloque)
    return h.hexdigest()


def huella_objeto(x):
    """Huella de una transcripción (repr estable de tuplas/listas/dicts de números)."""
    return hashlib.sha256(repr(x).encode('utf-8')).hexdigest()


def reporte():
    L = ['REGISTRO DE FUENTES DE VALIDACION (generado por modelo/fuentes.py; no editar a mano)',
         'Protocolo: 1 identidad (DOI, arbitrada, SHA-256, erratas) - 2 datos con ubicacion - 3 consistencia',
         'interna - 4 contraste - 5 hallazgos - 6 decision. Detalle en la cabecera de fuentes.py.', '']
    for k, f in FUENTES.items():
        arch = f['archivo']
        estado_pdf = 'sin archivo local' if not arch else (
            'OK' if os.path.exists(os.path.join(RAIZ, arch)) and sha256(os.path.join(RAIZ, arch)) == f['sha256']
            else 'NO COINCIDE / FALTA')
        L += [f'[{k}] {f["cita"]}',
              f'    DOI {f["doi"]} | arbitrada: {"si" if f["arbitrada"] else "NO"} | {f["ambito"]} | {f["acceso"]}',
              f'    formato: {f["formato"]} | erratas: {f["erratas"]} | PDF: {estado_pdf}',
              f'    usos: {"; ".join(f["usos"])}']
        for nom, (v, ubic) in f['datos'].items():
            L.append(f'    dato {nom} = {v}  ({ubic})')
        for nom, (mod, attr, h) in f.get('transcripciones', {}).items():
            L.append(f'    transcripcion {nom}: {mod}.{attr} (huella {h[:12]}...)')
        L += [f'    consistencia: {c}' for c in f['consistencia']]
        L += [f'    HALLAZGO: {h}' for h in f['hallazgos']]
        L += [f'    DECISION: {f["decision"]} | PEDIGRI {f["pedigri"]}', '']
    return '\n'.join(L) + '\n'


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    texto = reporte()
    with open(os.path.join(RAIZ, 'modelo', 'resultados', 'fuentes_registro.txt'), 'w', encoding='utf-8') as fh:
        fh.write(texto)
    print(texto)


if __name__ == '__main__':
    main()
