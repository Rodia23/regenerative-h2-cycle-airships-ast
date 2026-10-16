# [EN] Single source of all the constants: airships, physics and equipment classes.
# =============================================================================
# parametros.py — ÚNICA fuente de constantes del modelo (pipeline, Paso 0).
# Todos los scripts importan de aquí; ningún script define valores propios.
# Para simular otro dirigible, edita la sección "DIRIGIBLES" y "CONTROL".
# =============================================================================

import math


# =============================================================================
# CONSTANTES FÍSICAS
# =============================================================================

# Valores PUBLICADOS (VAL-01, val01_componentes.py los compara con su fuente). Adoptados el
# 2026-09-24 en lugar de los redondeados (9.81, 8.314, 120 MJ/kg, 1:8:9, L_v a 100 °C, H_s 8500,
# 45 kWh/kg); los valores viejos están en el historial de git.
g             = 9.80665    # [m/s²]      Gravedad estándar (ISA, ISO 2533)
R_univ        = 8.314462618  # [J/mol·K] Constante universal de los gases (CODATA 2018)
M_H2          = 2.01588e-3 # [kg/mol]    Masas molares (IUPAC 2005: H 1.00794, O 15.9994)
M_O2          = 31.9988e-3 # [kg/mol]
M_H2O         = 18.01528e-3  # [kg/mol]
M_He          = 4.002602e-3  # [kg/mol]   Helio (IUPAC), solo para la comparación de EXP-12
R_H2          = R_univ / M_H2   # [J/kg·K]  ≈ 4124.5
R_aire        = 287.05287  # [J/kg·K]    Constante específica del aire seco (ISA)
cp_H2         = 14300      # [J/kg·K]    Calor específico del H2 (NIST: 14.27 a 15 °C, 14.31 a 27 °C)
Q_comb_H2     = 241.826e3 / M_H2   # [J/kg]  PCI del H2 = 119.96 MJ/kg (CODATA: agua como vapor a 25 °C)
L_v_agua      = 2441.7e3   # [J/kg]      Calor latente del agua a 25 °C (NIST, IAPWS-95): la temperatura
                           #             de referencia del PCI; así PCI + (M_H2O/M_H2)·L_v = PCS de CODATA
ratio_O2_H2   = 0.5 * M_O2 / M_H2   # [—] Estequiometría: 1 kg H2 + 7.937 kg O2 -> 8.937 kg H2O
                           #             (el agua por kg de H2 es 1 + ratio_O2_H2 = M_H2O/M_H2)
cp_N2         = 1041       # [J/kg·K]    Calor específico del N2 (NIST, 27 °C)
cp_agua_liq   = 4180       # [J/kg·K]    Calor específico del agua líquida (NIST, 25-40 °C)
L_fusion_hielo = 333.55e3  # [J/kg]      Calor de fusión del hielo a 0 °C y 1 atm (6.009 kJ/mol): valor de manual. Solo
                           #             para la cota del hielo de la auditoría 2 (M6, sección 5.2); no cambia el modelo
cp_vapor      = 1920       # [J/kg·K]    Calor específico del vapor saturado (NIST: 1912 a 25 °C, 1931 a 40 °C)
T_ref_comb    = 25.0       # [°C]        Temperatura de referencia del PCI y de L_v
X_O2_aire     = 0.209476   # [—]         Fracción MOLAR de O2 en el aire seco (US Standard Atmosphere 1976)
M_aire        = 28.9644e-3 # [kg/mol]    Masa molar del aire seco (US Standard Atmosphere 1976)
w_O2_aire     = X_O2_aire * M_O2 / M_aire   # [—] Fracción másica de O2 ≈ 0.2314 (aire estequiométrico
                           #             = ratio_O2_H2 / w_O2_aire ≈ 34.3 kg de aire por kg de H2;
                           #             el N2 y el resto entran y salen, solo el O2 queda como agua)
M_inerte_aire = (M_aire - X_O2_aire * M_O2) / (1.0 - X_O2_aire)   # [kg/mol] Aire sin O2 (N2 + Ar),
                           #             ≈ 28.16 g/mol, para la dilución del vapor en el escape (EXP-05)
def k_lamb(esbeltez):
    """Masa añadida TRANSVERSAL / masa de aire desplazado de un esferoide prolato de esbeltez L/D (Lamb,
    Hydrodynamics §373; Tuckerman, NACA Report 210, 1926, p. 9): k2 = beta0/(2 - beta0). 0.804, 0.860 y
    0.917 en esbeltez 3, 4 y 6 (fuentes.TUC26)."""
    e = math.sqrt(1.0 - 1.0 / esbeltez ** 2)
    beta0 = 1.0 / e ** 2 - (1.0 - e ** 2) / (2.0 * e ** 3) * math.log((1.0 + e) / (1.0 - e))
    return beta0 / (2.0 - beta0)


k_masa_anadida = k_lamb(4.0)   # [—]     k de REFERENCIA (esbeltez 4, 0.860) para las fórmulas genéricas
                           #             (EXP-02/05/06, VAL-01). Cada nave usa el SUYO, con su esbeltez
                           #             (PLATAFORMAS[...]['k'], desde 2026-09-28, pendiente S2-2). Entra en
                           #             la INERCIA (inercia.py, EXP-13), no en el peso.
phi_agua_retenida = 1.0    # [—]         Fracción del agua de la combustión que queda a bordo AL INSTANTE
                           #             (1 = supuesto del artículo; 0 = el escape sale con el vapor).
                           #             L(phi) = rho_aire/rho_H2 - 1 + 9·phi  (22.3 con phi = 1; 13.4 con 0)



# =============================================================================
# ATMÓSFERA
# =============================================================================

rho_aire_std  = 1.225      # [kg/m³]   Densidad del aire en la altitud de referencia h0
                          #           (estado de referencia = ISA nivel del mar; operación costera,
                          #           h0 sobre el suelo. La ISA a 50 m da 1.2191, -0.5 %: solo mueve
                          #           el arrastre; L depende de rho_aire/rho_H2 a igual T y P)
                          #           rho(h) = rho_aire_std · exp(-(h - h0)/H_escala)
P_atmos       = 101325     # [Pa]      Presión atmosférica (ISA nivel del mar)
T_ref_atm     = 288.15     # [K]       Temperatura ISA a nivel del mar
                          #           Coherente: rho_aire_std = P_atmos/(R_aire·T_ref_atm) = 1.2250
gradiente_ISA = 0.0065     # [K/m]     Gradiente térmico de la troposfera ISA
H_escala      = 1.0 / (g / (R_aire * T_ref_atm) - gradiente_ISA / T_ref_atm)   # [m] ≈ 10 416
                          #           Altura de escala LOCAL de la densidad ISA: la exponencial
                          #           reproduce la ISA a primer orden en la banda de la maniobra.
                          #           (Antes 8500 m, que es la escala de la PRESIÓN isoterma RT/g;
                          #           efecto en la simulación del artículo: 0.2 mm, VAL-01.)


# =============================================================================
# DIRIGIBLES — las 4 plataformas del paper 1 (Tabla 2)
# =============================================================================
# Decisión del autor 2026-09-28 (hallazgo S3 de la auditoría; datos y fuentes de cada nave en fuentes.py).
# Antes cada nave se definía por una masa de estructura sin fuente y el volumen salía de ella (el "LZ 129"
# tenía 232 000 m3 y más masa que su sustentación real); el Pathfinder 3 salió (sin datos citables).
#
# La dinámica depende SOLO del volumen del casco: en equilibrio M_I = (1 + k)·rho_aire·V_casco, pese lo que
# pese la nave (lo que no es gas es aire de balonets, que se acelera con ella; inercia.py). Todas van llenas
# de H2 en el modelo; las de helio (LCA60T, NT) son variantes NOCIONALES con H2 del mismo casco.
#   V_casco   volumen TOTAL del casco [m3] (celdas de gas + aire entre ellas y en balonets)
#   m_vacio   peso vacío publicado [kg]; None = no publicado: se deriva abajo (el máximo que flota con la
#             liberación de referencia y el casco lleno de H2)
#   carga_max carga útil comercial [kg]
#   L, D      largo y diámetro máximo [m]: de ahí salen la esbeltez, k (Lamb), el área de arrastre y la
#             superficie del casco (bloque "campos derivados", más abajo)
#   color     paleta Okabe-Ito (CVD-safe), la misma en todas las figuras
#   corto     sufijo de las macros de numeros.tex (solo letras)
FT  = 0.3048               # [m]   pie (las fuentes históricas están en pies y libras)
FT3 = FT ** 3              # [m3]  pie cúbico
LB  = 0.45359237           # [kg]  libra
GAL = 3.785411784e-3        # [m3]  galón de EE. UU. (la norma de los tanques de agua, MIL-T-53029C)
OZ  = 0.028349523125       # [kg]  onza avoirdupois
YD2 = (3.0 * FT) ** 2      # [m2]  yarda cuadrada
PLATAFORMAS = {
    # Goodyear Aerospace, Feasibility Study of Modern Airships, Phase II, Vol. III: Historical Overview, 1975
    # (NASA, NTRS 19750023961), pp. 19-20: volumen total 7 650 000 ft3, peso vacío 249 000 lb, 814 x 135.2 ft.
    # Carga: "total payload" 41 990 lb, U.S. Civil Aeronautics Authority, The Hindenburg accident, 1938
    # (DOI 10.21949/1511799; fuentes.DOC37)
    'LZ 129 Hindenburg': dict(V_casco=7_650_000 * FT3, m_vacio=249_000 * LB, carga_max=41_990 * LB,
                              L=814 * FT, D=135.2 * FT, color='#009E73', corto='Hind'),
    # flying-whales.com: página del LCA60T (180 000 m3 de helio; la FAQ dice 200 000 m3, se declara) y FAQ
    # (peso vacío 100 t, "estructura y sistemas"; 200 x 50 m; 60 t de carga)
    'Flying Whales LCA60T': dict(V_casco=180_000.0, m_vacio=100_000.0, carga_max=60_000.0,
                                 L=200.0, D=50.0, color='#0072B2', corto='FW'),
    # Goodyear 1975, pp. 19-20, versión con H2: volumen total 2 800 000 ft3, peso vacío 86 400 lb, 660 x 90.7 ft.
    # Carga comercial 5 t: L. Dürr, The American airship ZR-3, NACA TM-286, 1924
    'ZR-3 Los Angeles': dict(V_casco=2_800_000 * FT3, m_vacio=86_400 * LB, carga_max=5_000.0,
                             L=660 * FT, D=90.7 * FT, color='#D55E00', corto='ZR'),
    # EASA TCDS AS.001 (LZ N07-100/-101), ed. 10, pp. 5 y 10: envolvente 8 425 m3, 75.1 x 14.2 m. Peso vacío no
    # publicado (la "max. gondola mass" es solo la góndola); carga "reserved": 1 t = la liberación de referencia
    'Zeppelin NT': dict(V_casco=8_425.0, m_vacio=None, carga_max=1_000.0, carga_publicada=False,
                        L=75.1, D=14.2, color='#E69F00', corto='ZNT'),
}

nombre_dirigible = 'Flying Whales LCA60T'   # caso base del artículo (autor, 2026-09-28)
V_casco_FW_alternativo = 200_000.0   # [m3] la otra cifra del fabricante para el LCA60T (FAQ, fuentes.FWH26):
                                     # sensibilidad del caso base (antes escrita en exportar_numeros.py)

m_carga       = 1000       # [kg]   Masa de la carga útil a liberar (liberación de referencia)
h0            = 50         # [m]    Altura inicial y referencia del controlador


# =============================================================================
# AERODINÁMICA
# =============================================================================

Cd            = 0.5        # [—]    Coeficiente de arrastre (cuerpo elipsoidal)


# =============================================================================
# HIDRÓGENO — almacenamiento y ballonets
# =============================================================================
# La presión, la masa y la compresión de los tanques están en el bloque de la regeneración
# (P_tanque_bar, kWhkg_tanque_sistema, wt_COPV, W_iso_20_350). Auditoría M6 (2026-09-29): salieron
# P_tanque, rho_H2_tanque (27.8 kg/m3 estaba mal: el NIST da 23.3 a 350 bar y 25 C), m_descomp y
# V_delta_max, que ningún script usaba.


# =============================================================================
# TÉRMICA — gas isotérmico
# =============================================================================

T1            = T_ref_atm  # [K]   Temperatura del H2 (isotérmico; = T_ambiente ISA nivel del mar, 15 °C)


# =============================================================================
# DIRIGIBLES — campos DERIVADOS (no se escriben a mano)
# =============================================================================

def superficie_esferoide(L, D):
    """Superficie de un esferoide prolato de largo L y diámetro D [m2]."""
    a, b = L / 2.0, D / 2.0
    e = math.sqrt(1.0 - (b / a) ** 2)
    return 2.0 * math.pi * b ** 2 * (1.0 + a / (b * e) * math.asin(e))


rho_H2_ref = P_atmos / (R_H2 * T1)   # [kg/m3] H2 en el estado de referencia (ISA nivel del mar)
for _d in PLATAFORMAS.values():
    _d['esbeltez'] = _d['L'] / _d['D']
    _d['k']        = k_lamb(_d['esbeltez'])                     # masa añadida de ESTA nave
    _d['Area']     = math.pi / 4.0 * _d['L'] * _d['D']          # [m2] área de arrastre en vertical (planta)
    _d['A_total']  = superficie_esferoide(_d['L'], _d['D'])     # [m2] superficie del casco (paneles)
    _d['m_vacio_publicado'] = _d['m_vacio'] is not None
    _d.setdefault('carga_publicada', True)          # NT: "reserved" en el TCDS (M-3, auditoría 3)
    if _d['m_vacio'] is None:                                   # el máximo que flota con la carga de referencia
        _d['m_vacio'] = (rho_aire_std - rho_H2_ref) * _d['V_casco'] - m_carga
    _d['sustentacion'] = (rho_aire_std - rho_H2_ref) * _d['V_casco']   # [kg] con H2 y el casco lleno
    if _d['m_vacio'] + _d['carga_max'] > _d['sustentacion'] * (1.0 + 1e-9) and _d['m_vacio_publicado']:
        raise ValueError(f"{_d['corto']}: peso vacío + carga máxima no flota con H2 en {_d['V_casco']:.0f} m3")
del _d

NAVE_BASE = PLATAFORMAS[nombre_dirigible]   # caso base (figuras, tabla de masas, cifras del abstract)
Area      = NAVE_BASE['Area']               # [m2]  área de arrastre del caso base
m_vacio   = NAVE_BASE['m_vacio']            # [kg]  peso vacío del caso base (las simulaciones lo cambian por nave)


# =============================================================================
# ENERGÍA — generador
# =============================================================================

rend_gen      = 0.30       # [—]   Rendimiento del generador eléctrico de H2 (30%)
                            # Contexto: E_generada = flujo_H2 · Q_comb_H2 · rend_gen
                            # NOTA: unificado a 30% (valor del código).
                            #   El artículo declarará este valor explícitamente.



# =============================================================================
# CONDENSACIÓN
# =============================================================================
# NOTA: NO hay Q_cond fijo. El artículo reporta Q_cond,min ≈ 1514 kW, que es
# el resultado de invertir  t_cond = m_H2O · L_v_agua / Q_cond  con t_cond = 10 min
# y m_H2O = 9 · (m_carga / L_eff)  (ver condensacion.py, que lo calcula aparte).
# Ref. de escala: compacto 10–50 kW, auto 50–200 kW, industrial 200–2000 kW.
t_cond_objetivo = 10.0     # [min]  Ventana objetivo para condensar el agua del descenso

# Condensación EN TIEMPO REAL (EXP-05, solo referencia de orden de magnitud, no diseño):
T_salida_cond  = 40.0      # [°C]    Temperatura a la que sale la mezcla N2 + vapor no condensado
                           #         (aire exterior "no más de 30-40 °C", supuesto del autor)
T_aire_exterior = 30.0     # [°C]    Aire exterior del caso de diseño del condensador (día caluroso):
                           #         enfría el escape y es el aire de combustión. El H2 entra a esa
                           #         misma temperatura (isotérmico con el ambiente); L no cambia con
                           #         la temperatura ambiente (depende de M_aire/M_H2 a igual T y P).
h_conv_casco   = (5.0, 25.0)  # [W/m²K] Convección exterior sobre el casco: calma (convección
                           #         natural, en vuelo estacionario) a forzada (viento/soplantes)
dT_casco       = 25.0      # [K]     Diferencia media pared-aire de un condensador que recorre el
                           #         casco (pared entre el rocío ~73 °C y la salida, aire <= 30 °C)

# Condensar TODO (EXP-06, criterio del autor): números mínimos, sin elegir equipo
phi_objetivo   = 0.99      # [—]     Fracción del agua a condensar ("todo", en la práctica)
cp_aire        = 1005      # [J/kg·K] Calor específico del aire de enfriamiento
dT_aire_enfr   = 30.0      # [K]     Calentamiento admisible del aire de enfriamiento (30 -> 60 °C)
dT_media_cond  = 40.0      # [K]     Diferencia media de temperatura del condensador (para UA)
dT_buffer_agua = 50.0      # [K]     Calentamiento del agua si se guardara el calor en vez de rechazarlo (EXP-06)

# Dinámica del condensador en las simulaciones del artículo (EXP-09, condensador.py)
tau_condensador = 5.0      # [s]     Retardo térmico de primer orden (tiempo hasta régimen; SUPUESTO,
                           #         rango a explorar 0-30 s)
T_pulmon       = 70.0      # [°C]    Temperatura del vapor retenido como gas en el modo "pulmón"
R_vapor        = 461.5     # [J/kg·K] Constante del vapor de agua


# =============================================================================
# CLASES DE REFERENCIA DE LAS CAJAS NEGRAS: rechazo de calor y su potencia de ventiladores
# ([P] = fuente primaria, [S] = secundaria)
# =============================================================================
gamma_aire     = 1.4       # [—]
kWkg_radiador  = (10.0, 20.0)  # [kW/kg] Radiador de aire, FlyZero 2026 / 2035 [P]
kWkg_TMS_sistema = (5.0, 10.0) # [kW/kg] SISTEMA de gestión térmica completo = calor rechazado / masa TOTAL del
                               #         sistema, FlyZero p. 4 [P]: 5 (2026, rechazo a 100 °C con compresión de
                               #         vapor) y 10 (2035, LT-PEM, rechazo a 100 °C con líquido; el 20 de 2035
                               #         es HT-PEM a 160 °C). Masa REFERENCIAL de la caja negra de rechazo del
                               #         paper 1 (tabla 6, EXP-11): hallazgos I5/S1 de la auditoría (2026-09-27)
                               #         reemplazan las dos bandas (EXP-14)
perdida_rad    = (0.20, 0.15)  # [—]     Potencia de bombeo/arrastre por kW rechazado, 2026 / 2035 [P]
                               #         (FlyZero FZO-PPN-COM-0019, radiador de aire: 20 % en 2026, 15 % en
                               #         2035, 10 % en 2050; antes decía 0.10 = 2050, mezclado con 20 kW/kg
                               #         de 2035. Orden: conservador, optimista, igual que kWkg_radiador)
# Revisión del 2026-10-01 (registro de revisión R-03): f de FlyZero incluye el arrastre de una instalación con aire
# de impacto en vuelo; una unidad con ventiladores en vuelo estacionario gasta solo la caída de presión de su
# intercambiador: f = dp / (eta_vent rho_aire cp_aire dT_aire). Los de FlyZero quedan como cota alta.
dp_rechazo_Pa  = (200.0, 1000.0)  # [Pa] caída de presión del lado del aire: SUPUESTO (de aletas holgadas a compacto)
eta_ventilador = 0.7       # [—]     rendimiento del ventilador: SUPUESTO
v_aire_entrada = 10.0      # [m/s]   velocidad del aire de enfriamiento en la entrada, solo para la escala del área:
                           #         SUPUESTO
# Revisión del 2026-10-01 (R-07): humedad del aire de entrada (el modelo usa aire seco; la sensibilidad la da
# val03_cotas.py, parte E)
HR_aire        = (0.5, 0.9)   # [—]  humedad relativa del día de diseño (30 °C, costa): SUPUESTO
# Revisión del 2026-10-01 (discusión, "Braking"): las hélices vectoriales que la nave ya lleva para maniobrar, fuera de
# la premisa del paper 1, sosteniendo la nave mientras se quema al ritmo del rechazo instalado (helices_vectoriales.py).
# Sin datos de propulsión verificados por nave: hélices genéricas
freno_A_discos = (20.0, 100.0)    # [m2] área total de discos de hélice: SUPUESTO
freno_FM       = 0.7       # [—]     figura de mérito de esas hélices: SUPUESTO
freno_P_bajo   = 2.0e6     # [W]     rechazo menor, para el compromiso potencia-tiempo: SUPUESTO
freno_escalon_grande = 10_000.0   # [kg] escalón grande, para la escala F^1.5 de las hélices: SUPUESTO
# Las hélices como dispositivo de sujeción en la simulación (drop_flujo_optimo_perdrop.correr(helices=...), freno.py):
# retardo de su fuerza y lazo de altitud que la manda
hel_tau  = 1.0     # [s]     retardo de primer orden de la fuerza de las hélices (giro y régimen): SUPUESTO
hel_wn   = 0.3     # [rad/s] frecuencia natural del lazo de altitud de las hélices: SUPUESTO
hel_zeta = 1.0     # [—]     amortiguamiento de ese lazo: SUPUESTO
HP_KW    = 0.745700  # [kW/hp] caballo de fuerza mecánico
hel_margen_autor = (0.30, 0.40)   # [—] márgenes de seguridad sobre la fuerza de las hélices (freno.py): SUPUESTO del autor
# El dispositivo de sujeción (freno.py, sección 5.2 del paper). El freno es la 4.ª caja negra: con un rechazo
# instalado P, el ciclo compensa en tiempo real solo el escalón m_e(P) y el freno toma el resto. Usa freno_A_discos,
# freno_FM (hélices ideales), hel_tau, hel_wn, hel_zeta (su lazo en la simulación) y hel_margen_autor[0] (su fuerza).
freno_P_grilla   = (2e6, 10e6)               # [W] rechazo instalado de la tabla del freno: el menor y el de referencia
                                             #     (P_ref_masa): SUPUESTO (rejilla). freno.py agrega el límite: la
                                             #     potencia con la que los haces de ACC pesan la carga útil del LCA60T
                                             #     (auditoría 2026-10-06, M2: antes 50 MW, que no cabía)
freno_escalones  = {'Flying Whales LCA60T': (1000.0, 10000.0, 30000.0),   # [kg] 30 t ~ un contenedor de 20 pies
                    'Zeppelin NT': (1000.0,)}                             #      (ISO 668: 30.48 t); 1 t = carga del NT
freno_t_fig_h    = 24.0                      # [h] tiempo de recuperación máximo de la figura
HELICES_NAVES = {
    'Zeppelin NT': {
        'n_vert': 3, 'D_m': 2.7, 'P_unidad_kW': 147.0, 'P_planta_kW': 3 * 147.0, 'vectorial': True,
        'liviandad_cert_kg': 200.0, 'pesadez_cert_kg': (400.0, 500.0),
        'nota': 'CONFIRMADO en FAA TCDS AS1CE rev. 3 (8 mar 2021), pp. 1-3, y EASA TCDS EASA.AS.001 ed. 10 (20 oct '
                '2025), secc. 1.5 y 1.10: 3 Lycoming IO-360-C1G6; 3 hélices vectoriales Hoffmann HO-V373 de 2.7 m, '
                'tripalas de paso variable, máx. 1250 rpm (flancos y popa), más una lateral HO-V123F de 2.2 m que no '
                'cuenta. Vectorización: flancos de 0° (adelante) a 120° hacia arriba; popa de 0° (atrás) a 90° hacia '
                'abajo (AS1CE p. 3). LIMITACIONES: liviandad estática máx. 200 kg; pesadez máx. 400 kg en despegue y '
                'aterrizaje (= la carga máx. del tren) y 500 kg en vuelo. Lastre de agua: 700 kg. La potencia, 147 kW '
                'por motor, NO está en estas hojas (está en el TCDS del motor, FAA 1E10 / EASA.IM.E.032): ficha de ZLT, '
                'por confirmar'},
    'Flying Whales LCA60T': {
        'n_vert': 32, 'D_m': 4.0, 'P_unidad_kW': 200.0, 'P_planta_kW': 4000.0, 'vectorial': True,
        'nota': 'CONFIRMADO en la FAQ de Flying Whales (texto guardado por el autor, 2026-10-01): 200 m x 50 m, '
                '~200 000 m3 de helio en 14 celdas, 60 t de carga útil, estructura y sistemas 100 t; propulsión híbrida '
                '(turbogenerador + baterías) con 32 motores eléctricos alrededor del casco; carga y descarga en vuelo '
                'estacionario con viento de hasta 40 km/h; la CARGA se intercambia con LASTRE DE AGUA (~60 t, tanques '
                'en tierra) y los 32 propulsores mantienen la posición. Visto en búsqueda, POR CONFIRMAR (página del '
                'propulsor): 7 puntos de propulsión, más de 200 kW continuos por motor (Evolito D250, 240 kW pico), '
                'hélice de 4 m tripala de paso variable, 4000 kW combinados, generadores de 1 MW (Honeywell). Rango de '
                'vectorización DESCONOCIDO: "todas verticales" es la cota superior'},
    'ZR-3 Los Angeles': {
        'n_vert': 0, 'D_m': None, 'P_unidad_kW': 400 * HP_KW, 'P_planta_kW': 5 * 400 * HP_KW, 'vectorial': False,
        'nota': '5 Maybach VL-I de 400 hp, hélices de paso fijo, empujadoras, con reversa; sin vectorización '
                '(NavSource). SIN VERIFICAR'},
    'LZ 129 Hindenburg': {
        'n_vert': 0, 'D_m': 6.0, 'P_unidad_kW': None, 'P_planta_kW': None, 'vectorial': False,
        'nota': '4 Daimler-Benz DB 602 diésel (potencia nominal en discusión: 1050-1320 hp por motor); hélices de 4 '
                'palas, paso fijo, ~6.0 m (fuente secundaria); sin vectorización. La recuperación de agua del escape '
                'NO se instaló en el LZ 129 (usaba hidrógeno); sí en el LZ 130 tras el accidente. SIN VERIFICAR'},
}
# Condensador enfriado por aire de una central (ACC), candidato de la tabla de brecha (req_masa_rechazo.py; ruta
# nueva del autor, 2026-09-28: requisito + brecha en vez de elegir equipo). NGCC de 646 MWe, verano [P]:
ACC_Q_MWth     = 388.9     # [MW]    calor rechazado (Zhang, Spatari & Sun 2020, tabla 3, verano; fuentes.ZSS20)
ACC_vent_MWe   = 4.01      # [MW]    potencia de ventiladores (misma tabla)
ACC_T_aire_C   = 33.33     # [°C]    aire de verano del sitio (Chicago, tabla del clima, p. 4)
ACC_ITD_C      = 23.6      # [K]     diferencia inicial de temperatura (tabla 3): condensa a ~57 °C
ACC_masas_t    = {'haces_acero': 550, 'haces_aluminio': 46, 'estructura': 1300, 'ductos': 650,
                  'tuberias_acero': 150, 'tuberias_inox': 3, 'bomba_aire': 5, 'mecanico_acero': 45,
                  'mecanico_prfv': 65, 'tanque': 25, 'electrico_cobre': 5, 'electrico_acero': 30,
                  'electrico_plastico': 5}   # [ton] suplementario, tabla S7 (fuentes.ZSS20S); "ton" sin aclarar
# Baterías a NIVEL PAQUETE (decisión del autor, 2026-09-25): valor de celda x fracción celda/paquete.
f_celda_paquete = 149.0 / 225.0  # [—] NASA X-57: celda 225 Wh/kg -> paquete 149 Wh/kg (66 %): BMS, térmica,
                           #         estructura e interconexiones (primaria: Chin et al. 2020, AIAA 2020-3567, p. 5:
                           #         55.3 kWh / 370.86 kg; también Pattanayak 2025, Appl. Energy 397, 126356; ver
                           #         fuentes.CHIN20). Se aplica igual a Wh/kg y a kW/kg (es masa).
Whkg_bateria_celda = 200.0 # [Wh/kg] Celda Li-ion típica: SUPUESTO dentro del rango de celda 90-260 Wh/kg
                           #         (Pattanayak 2025); la celda de la X-57 tiene 225
kWkg_bateria_celda = (1.0, 2.0)  # [kW/kg] Celda Li-ion, conservador / optimista (Pattanayak 2025: 1-2 kW/kg;
                           #         celda de alta potencia 1.8 kW/kg, IEEE 10598954)
Whkg_bateria   = Whkg_bateria_celda * f_celda_paquete                  # [Wh/kg] paquete
kWkg_bateria   = tuple(x * f_celda_paquete for x in kWkg_bateria_celda)  # [kW/kg] paquete (buffer, EXP-11)
# SENSIBILIDAD, no clase representativa (autor 2026-09-29, ítem 8): la mayor potencia específica encontrada para
# un pulso de 30 s, celda Amprius SiCore SA88 de ultra alta potencia (catálogo Q3 2026, ficha AMP-TDS-SA88-R3,
# fuentes.AMP26) [P]; supera a la Ionblox (3.8 kW/kg, fuentes.IBX22). Solo dice cuánto baja el requisito del
# rechazo con la mejor celda de hoy (req_masa_rechazo.py)
AMP_SA88 = {'I_pulso_30s_A': 153.0, 'V_nom_V': 3.45, 'masa_g': 96.5, 'Whkg': 365.0}
kWkg_celda_30s_mejor = AMP_SA88['I_pulso_30s_A'] * AMP_SA88['V_nom_V'] / AMP_SA88['masa_g']   # [kW/kg] = W/g; con
                               # el voltaje NOMINAL: cota superior (bajo 15C el voltaje cae)
Whkg_celda_30s_mejor = AMP_SA88['Whkg']   # [Wh/kg] celda, a C/5
# Contraste del factor celda -> paquete (0.66, X-57): 25 paquetes de autos eléctricos 2010-2019, gravimétrico
# celda -> sistema 0.484-0.742 (Löbberding et al. 2020, tabla 1; fuentes.LOB20) [P]
f_celda_paquete_bev = (0.484, 0.742)

# Balance de masa por caja negra del paper 1 (balance_masa.py): una clase de tecnología
# representativa por caja, sin elegir equipo.
TURBOEJE_TCDS  = ((1879e3, 243.6), (1212e3, 194.6))  # [W, kg] CT7-8 y CT7-2A, EASA TCDS IM.E.010 [P]:
                           #         unidad de conversión (turbina de gas; sin el generador eléctrico)
#                          Paper 1 desde 2026-09-27 (hallazgo I2 de la auditoría, decisión del autor): la unidad
#                          de conversión es de la clase MOTOR ALTERNATIVO de H2 a lambda = 1 (Verhelst & Wallner),
#                          coherente con el escape estequiométrico; la turbina (exceso de aire) queda como nota.
MOTOR_ROTAX    = ((117e3, 85.8), (101e3, 91.2))  # [W, kg] Rotax 916 iS/c (gasolina; Bosch lo pasó a H2, 115 kW, sin
                           #         masa publicada: PROXY de aviación). 117 kW de despegue con el motor solo (reductora,
                           #         turbo, intercooler, escape); 101 kW continuos con bastidor, alternador y deflector
                           #         (+5.4 kg). Folleto de producto 2023, p. 2 [S]
MOTOR_DEUTZ    = ((220e3, 770.0),)  # [W, kg] DEUTZ TCG 7.8 H2, primer motor de H2 de serie (tren, estacionario,
                           #         maquinaria): 220 kW (ISO 14396), ~770 kg (DIN 70020 parte 7A). Ficha 2024, p. 2 [S];
                           #         ~19 kg/h de H2 a plena potencia (deutz.com, eta ~0.35)
MOTORES_ALTERNATIVOS = MOTOR_ROTAX + MOTOR_DEUTZ   # rango de la clase: aviación (proxy) a industrial de H2
kWkg_PEM       = (0.5, 2.2)    # [kW/kg] Electrolizador PEM, pesimista / optimista. PROXY: pilas PEM de
                           #         aviación con su balance de planta (Zhang 2025, Hydrogen 6:116) [S-proxy]
kWhkg_tanque_sistema = 1.71  # [kWh/kg] Capacidad gravimétrica del SISTEMA Tipo IV a 700 bar (5.6 kg útiles),
                           #         ESTADO 2023 del DOE: Program Record 24006 (2025), tabla 1 [S]. La meta 2025
                           #         (1.8, DOE17) no se alcanzó (autor 2026-09-29: lo que existe, no la meta).
TANQUE_DOE24   = {'P_bar': 700.0, 'P_min_bar': 15.0, 'T_C': 15.0, 'V_L': 144.2, 'D_cm': 39.6,
                  'util_kg': 5.6, 'guardado_kg': 5.8, 'compuesto_kg': 72.7}
                           #         Diseño 2023 de ese mismo sistema: Record 24006, tabla 2 (p. 4) y fig. 1 (p. 5,
                           #         valor impreso sobre la barra) [S]. Diámetro fijo, largo ajustado (p. 3).
TANQUE_DOE19_FIJOS = {'boss': 0.9, 'valvula': 3.0, 'regulador': 3.6, 'otros_BOP': 15.4}
                           # [kg]    Componentes que no escalan con la presión: Record 19008 (2019), tabla supl. 5
                           #         (p. 12) [S]; el diseño 2023 conserva los mismos componentes de planta (Record
                           #         24006, p. 3). El liner 2023 sale por diferencia (el PA de 0.38 cm reemplazó al
                           #         HDPE de 0.5 cm y 8.0 kg de 2019).
NIST_RHO_H2_15C = {1.5: 1.2509, 35.0: 23.995, 70.0: 40.172}
                           # [kg/m3] H2 a 15 °C por presión [MPa]: NIST WebBook, SRD 69 (Leachman et al. 2009),
                           #         consulta 2026-09-29. Reproducen la tabla 2 del Record 24006 (5.61 kg útiles).


def tanque_a_presion(P_bar):
    """Sistema Tipo IV del diseño 2023 del DOE (700 bar) llevado a la presión P con las reglas del mismo diseño:
    diámetro fijo y largo ajustado, compuesto proporcional a P·V (pared delgada, misma fibra y factor de seguridad),
    liner proporcional a su área (domos hemisféricos: SUPUESTO), el resto de los componentes igual y el H2 no usable
    (a la presión mínima) proporcional al volumen. A 700 bar devuelve exactamente el sistema del Record 24006.
    Devuelve masas [kg], volumen [L] y w = H2 útil / masa del sistema (la definición del DOE)."""
    d, rho = TANQUE_DOE24, NIST_RHO_H2_15C
    r_min = rho[d['P_min_bar'] / 10.0]
    M_ref = d['util_kg'] * Q_comb_H2 / 3.6e6 / kWhkg_tanque_sistema         # 109.1 kg (1.71 kWh/kg)
    fijos = sum(TANQUE_DOE19_FIJOS.values())
    liner_ref = M_ref - d['compuesto_kg'] - d['guardado_kg'] - fijos       # por diferencia
    V = d['V_L'] * (rho[d['P_bar'] / 10.0] - r_min) / (rho[P_bar / 10.0] - r_min)
    r = d['D_cm'] / 200.0                                                   # [m]

    def area(V_L):
        V_dom = 4.0 / 3.0 * math.pi * r ** 3
        return 4.0 * math.pi * r ** 2 + 2.0 * math.pi * r * (V_L / 1e3 - V_dom) / (math.pi * r ** 2)

    compuesto = d['compuesto_kg'] * P_bar * V / (d['P_bar'] * d['V_L'])
    liner = liner_ref * area(V) / area(d['V_L'])
    h2 = d['guardado_kg'] + (V - d['V_L']) / 1e3 * r_min
    M = compuesto + liner + fijos + h2
    return {'V_L': V, 'compuesto_kg': compuesto, 'liner_kg': liner, 'liner_ref_kg': liner_ref, 'fijos_kg': fijos,
            'H2_kg': h2, 'M_kg': M, 'w': d['util_kg'] / M}
PV_TEC         = {'a-Si': (0.08, 0.120),      # (eficiencia, kg/m2): amorfo, Pande & Verstraete 2018 [P] (arreglo)
                  'CIGS': (0.146, 2.4)}       # módulo flexible comercial MiaSolé FLEX-02W [P]: 380 W sobre 2.598 m2
                           #         = 14.6 % sobre el ÁREA DEL MÓDULO (el 17 % de la hoja es de apertura, solo las
                           #         celdas; decisión del autor 2026-09-25, fuentes.MIA14); 2.4 kg/m2 con adhesivo
# Compresión del H2 regenerado hacia los tanques (EXP-16; arquitectura (b) del autor, 2026-09-27:
# electrolizador -> compresor -> tanques = acumulador -> válvula -> celdas de gas).
W_iso_20_350   = 1.05      # [kWh/kg] Trabajo isotérmico mínimo de 20 a 350 bar (exergía, NIST), DOE Hydrogen
                           #         Program Record 9013 (2009), p. 1 [P]; se parte de 20 bar como el DOE (PEM
                           #         entrega hasta ~30 bar)
eta_compresion_H2 = 0.52   # [—]     Eficiencia de compresión de un compresor alternativo típico para 350 bar
                           #         (56 % isentrópica x 92 % del motor; modelo HDSAM), Record 9013, p. 2 [P]
e_comp_medido  = 3.1       # [kWh/kg] Promedio MEDIDO en las estaciones del programa del DOE (2-4 sin atípicos),
                           #         Record 9013, p. 2 [P]: extremo pesimista
P_tanque_bar   = 350.0     # [bar]   Tanques Tipo IV
                           #         Arquitectura (b) del autor (aclarada 2026-09-27): TODO el H2 regenerado se
                           #         comprime a los tanques, al ritmo del electrolizador, para que no cambie la
                           #         flotabilidad; nada va directo a las celdas. El agua se guarda como agua y se
                           #         regenera cuando hace falta. Los tanques salen llenos de tierra.
wt_COPV        = tanque_a_presion(P_tanque_bar)['w']   # [—] H2 útil / masa del sistema a 350 bar = 0.054
                           #         (autor 2026-09-29: con los datos 2023 del DOE, no con la razón 350/700 de 2010;
                           #         esa razón, ANL-10/24 = HUA10, queda como chequeo independiente en el grupo L)
t_regeneracion_h = 24.0    # [h]     Regenerar el H2 de una tonelada en un día (el mismo que balance_masa); el
                           #         compresor trabaja las mismas horas (su caudal = el del electrolizador)
COMPRESOR_REF  = {'P_bar': 414.0, 'flujo_kgph': (2.69, 2.94), 'motor_kW': 11.0, 'masa_kg': (545.0, 545.0)}
                           #         Compresor ELÉCTRICO de catálogo del tamaño pedido: Hydro-Pac LX C06-15 (15 HP =
                           #         11 kW, 6000 psi máx.): 19.2 scfm con entrada a 300 psig (~22 bar abs.) y 21 a
                           #         450 psig, x 0.14 kg/h por scfm; 1200 lb = 545 kg con aceite (fuentes.HPAC) [S].
                           #         El más chico de 6000 psi que alcanza el caudal con ~20 bar de entrada.
                           #         Equipo industrial (bastidor, reservorio): cota superior de la masa en vuelo
COMPRESOR_BAUER = {'P_bar': 350.0, 'flujo_kgph': (1.78, 2.22), 'motor_kW': 16.0, 'masa_kg': (1066.0, 1814.0)}
                           #         Contraste: compresor industrial más grande (Bauer H15.2-20), fuentes.BAUER [S]
TANQUE_AGUA    = {'capacidad_gal': 3000.0, 'lado_vacio_ft': 14.0, 'caras': 2, 'tela_ozyd2': (30.0, 62.0)}
                           #         Lastre de agua (D, 2026-10-13): el agua retenida se guarda en un tanque de tela
                           #         recubierta colapsable de agua potable, MIL-T-53029C (DoD, 1994; fuentes.MILT53029)
                           #         [N]: el de 3000 gal es una almohada de 14 x 14 ft vacía (fig. 1, dos caras) y su
                           #         tela pesa de 30 a 62 oz/yd2 (tabla II, ítem 1). Solo la tela: los accesorios
                           #         (válvulas, asas, parches) no se cuentan, cota optimista como las demás cajas
SGT800 = {'P_MWe': (47.5, 50.5, 54.0, 57.0), 'PR': (20.1, 21.0, 21.4, 21.8),
          'm_escape_kgps': (132.8, 134.2, 135.5, 136.6), 'heat_rate_kJkWh': (9547.0, 9389.0, 9206.0, 8970.0),
          'masa_kg': 285000.0}
                           #         Turbina de gas industrial Siemens SGT-800, ficha de 2017, p. 2 (fuentes.SGT800) [S]:
                           #         las cuatro potencias con su relación de presión, caudal de escape y heat rate (ISO,
                           #         gas natural); masa aprox. del paquete clásico, con el generador. Escala de la masa
                           #         del compresor del escape a 14 atm (críticas externas, 2026-10-11): todo el paquete
                           #         se atribuye a su compresor y se escala por el trabajo de compresión
PCI_CH4        = 50.0e6    # [J/kg]  PCI del metano (802.3 kJ/mol; CODATA para CO2 y H2O, JANAF para el CH4): solo para
                           #         restar el combustible del caudal de escape del SGT-800 y obtener su aire (un gas
                           #         natural de 45-50 MJ/kg cambia ese aire en menos de 0.3 %)
T_ISO_C        = 15.0      # [°C]    Condiciones ISO de las turbinas de gas (ISO 3977-2: 15 °C, 101.3 kPa), las de la
                           #         ficha del SGT-800 ("ISO ambient conditions", p. 2)
eta_motor_H2   = (0.30, 0.45)  # [—] Eficiencia de motores de H2: 0.30 conservador; 0.40-0.45 pico al freno de
                           #         motores de combustión de H2 (Verhelst & Wallner 2009) [P]
P_ref_masa     = 10e6      # [W]     Potencia de rechazo de referencia del balance de masa (la que cabe)
n_tau_arranque = 3.0       # [—]     El condensador arranca 3 tau antes de la quema (95 % de su capacidad)

# =============================================================================
# AIRE SECO (composición para la temperatura de llama de derivaciones.py, con janaf.py)
# =============================================================================
# Aire seco, fracciones molares (US Standard Atmosphere 1976, tabla 3); los gases traza
# (Ne, He, Kr, Xe, CH4, H2: 2.7e-5) se suman al Ar (monoatómico casi todo).
AIRE_X = {'N2': 0.78084, 'O2': X_O2_aire, 'CO2': 0.000314}
AIRE_X['Ar'] = 1.0 - sum(AIRE_X.values())

# =============================================================================
# CONTROL — controlador PD de altitud
# =============================================================================

Kp              = 15      # [m³/m]      Ganancia proporcional
# Kd NO es un parámetro fijo: se calcula en cada paso como
#   Kd_auto = 2·zeta·sqrt(Kp·M_eff/((rho_aire - rho_H2)·g))   (eq:kd_app, appendix B)
# con la rigidez real del lazo, para que zeta sea su amortiguamiento (2026-10-09; antes rho_aire·g y 1.35,
# que daban 1.302). Rango real ~0.2–1.6e3 según plataforma/fase. No se
# usa ningún Kd fijo en la simulación.
K_mass_control  = 60       # [kg/s/m³]   Conversión volumen→flujo másico (fase emergencia)
K_mass_trim     = 0.5      # [kg/s/m³]   Misma conversión en la fase trim (asentamiento suave)
zeta            = 1.30     # [—]         Coef. de amortiguamiento del PD adaptativo (sobreamortiguado); ANTES 1.35

max_flujo_H2    = 30      # [kg/s]   Límite del controlador para simulaciones base y multi-drop.
                            # Ref: el barrido paramétrico (t3_barrido_reactivo.py) evalúa 1–160 kg/s.
                            #   ṁmin resultante: 5 kg/s (Hindenburg) a 40 kg/s (Zeppelin NT).
                            #   Rango industrialmente viable con quemadores H2 comerciales.


# =============================================================================
# MANIOBRA PROACTIVA
# =============================================================================

dm_H2_proactivo = -10      # [kg/s]   Consumo de H2 antes de la liberación (esquema estándar)
t_anticip_estandar = 1.0   # [s]      Ventana del pulso proactivo (esquema estándar)
D_lim_oper      = 0.10     # [m]      Profundidad máxima de precondicionamiento (operacional)


# =============================================================================
# UMBRALES DE LA MISIÓN
# =============================================================================

h_max_precision = 1.0      # [m]   |Δh| máximo para acople de precisión
h_max_seguridad = 5.0      # [m]   |Δh| máximo de seguridad (carga externa)
intervalo_multidrop = 30.0 # [s]   Tiempo entre entregas en misiones multi-drop


# =============================================================================
# SIMULACIÓN
# =============================================================================

dt_max          = 0.05     # [s]   Paso máximo del integrador ODE
t_post_liberacion = 400    # [s]   Duración de la sim. después de liberar la carga


# =============================================================================
# ENERGÍA SOLAR — irradiancia MEDIDA (SENAMHI 2003) + factor de forma del casco
# (irradiancia_casco.py). Reemplaza al viejo rendimiento_efectivo = 2.0 kWh/m2/día, que no
# salía de SENAMHI (decisión del autor, 2026-09-25; fuentes.SEN03).
# =============================================================================

SENAMHI_GHI_ANUAL = 4.06      # [kWh/m²/día] Energía solar diaria, promedio anual, estación A. von Humboldt
                               #   (La Molina, Lima; piranómetro 1968-1999), SENAMHI 2003, tabla 2, p. 16
SENAMHI_LAT       = -12.08    # [°]  Latitud de esa estación (tabla 1)
SENAMHI_KT_MES    = (0.50, 0.50, 0.48, 0.47, 0.40, 0.295, 0.295, 0.287, 0.327, 0.40, 0.43, 0.52)
                               # [—]  Transmisividad H/H0 mensual (ene-dic) de esa estación, anexo 6, p. 29,
                               #      LEÍDA de la figura (+-0.01). Con H0 da 4.09 de promedio anual (tabla: 4.06)
G_SC              = 1367.0    # [W/m²] Constante solar de los atlas de esa época (WMO 1981); la actual (1361,
                               #      Kopp & Lean 2011) bajaría H0 un 0.4 %
albedo_suelo      = 0.2       # [—]  Reflectancia del suelo, valor estándar del modelo isotrópico (SUPUESTO;
                               #      aporta 0.2 x 0.18 = 3.6 % del GHI a la mitad superior del casco)

A_paneles_fraccion   = 0.45   # [—]   Fracción del área superficial cubierta con paneles
                               # SUPUESTO: mitad superior del casco sin las puntas (Sun 2024 da la geometría,
                               # no la fracción; con su casco sale 0.449, fuentes.SUN24)
                               # Rango típico: 0.30 – 0.60

eta_panel            = 0.20   # [—]   Eficiencia de conversión del panel fotovoltaico
                               # CIGS thin-film aeroespacial: 0.17–0.22
                               # Amorfos Si (Pande 2018):    0.07–0.10


# =============================================================================
# REGENERACIÓN — electrólisis
# =============================================================================

eta_electrolisis     = 0.63   # [—]       Eficiencia de la electrólisis sobre el PCI: promedio mundial
                              #           2024 (IEA, Global Hydrogen Review 2025, anexo de supuestos, p. 6).
                              #           Única fuente: la energía por kg sale de aquí.
E_H2_electrolisis    = Q_comb_H2 / 3.6e6 / eta_electrolisis   # [kWh/kg] ≈ 52.9
                              #           Piso físico: 39.4 kWh/kg (termoneutral, PCS). El artículo usaba
                              #           45 kWh/kg (= 74 % del PCI, extremo optimista del rango 45-60).
