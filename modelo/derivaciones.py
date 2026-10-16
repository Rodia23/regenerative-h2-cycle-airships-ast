# -*- coding: utf-8 -*-
# [EN] Intermediate steps of the equations of the article, recomputed independently (Appendices A-C).
# =============================================================================
# derivaciones.py — pasos intermedios de cada ecuación del paper 1, para el
# documento de derivaciones (Apéndices A-C del paper). Recalcula cada
# paso a partir de parametros.py con las fórmulas escritas tal cual en el
# documento; verificacion.py (grupo L) los compara con el modelo y con
# comprobaciones numéricas independientes. Salida: datos/derivaciones.json con
# la lista de macros 'Der...' que exporta exportar_numeros.py.
# =============================================================================

import math
import parametros as p
import rutas
import minimax as mm
import inercia
from exp05_condensacion_tiempo_real import p_sat, q_rech, N_INERTE, AIRE_ESTEQ, T_rocio
from exp06_requisitos_minimos import presion_para_phi, T_para_phi

PLAT = p.nombre_dirigible             # caso base del artículo
M_CARGA = 1000.0


def main():
    v = {}                              # valores (para verificacion.py)
    macros = []                         # (nombre, valor, decimales, comentario)

    def m(nombre, valor, dec, comentario):
        v[nombre] = valor
        macros.append((nombre, valor, dec, comentario))

    # ── Entradas ────────────────────────────────────────────────────────────
    rho_a = p.rho_aire_std
    rho_h = p.P_atmos / (p.R_H2 * p.T1)
    k_rho = rho_a / rho_h
    r = 1.0 + p.ratio_O2_H2
    m('DerG', p.g, 5, 'g [m/s2]')
    m('DerRuniv', p.R_univ, 4, 'R universal [J/mol/K]')
    m('DerMHdos', p.M_H2 * 1e3, 5, 'M_H2 [g/mol]')
    m('DerMHdosO', p.M_H2O * 1e3, 5, 'M_H2O [g/mol]')
    m('DerRHdos', p.R_H2, 1, 'R_H2 [J/kg/K]')
    m('DerPatm', p.P_atmos, 0, 'P_atm [Pa]')
    m('DerTuno', p.T1, 2, 'T_1 del H2 [K]')
    m('DerRhoAire', rho_a, 3, 'rho_aire [kg/m3]')
    m('DerRhoHdos', rho_h, 5, 'rho_H2 = P/(R_H2 T_1) [kg/m3]')
    m('DerRazonRho', k_rho, 2, 'rho_aire / rho_H2')
    m('DerR', r, 3, 'r = M_H2O/M_H2')
    m('DerPCI', p.Q_comb_H2 / 1e6, 2, 'PCI [MJ/kg]')
    m('DerLv', p.L_v_agua / 1e3, 1, 'L_v a 25 C [kJ/kg]')

    # ── Altura de escala ────────────────────────────────────────────────────
    m('DerRaire', p.R_aire, 2, 'R del aire [J/kg/K]')
    m('DerTatm', p.T_ref_atm, 2, 'T ISA nivel del mar [K]')
    m('DerGamma', p.gradiente_ISA * 1e3, 1, 'gradiente ISA [K/km]')
    m('DerHsIsoterma', p.R_aire * p.T_ref_atm / p.g, 0, 'R T / g [m]')
    m('DerHs', 1.0 / (p.g / (p.R_aire * p.T_ref_atm) - p.gradiente_ISA / p.T_ref_atm), 0, 'H_s local [m]')

    # ── L(phi) y hidrógeno por tonelada ─────────────────────────────────────
    L1, L0 = k_rho - 1.0 + r, k_rho - 1.0
    L_tanques = k_rho                       # O2 de tanques: la masa no cambia, solo se pierde flotación
    m('DerLuno', L1, 2, 'L(1)')
    m('DerLcero', L0, 2, 'L(0)')
    m('DerLtanques', L_tanques, 2, 'L con O2 de tanques')
    m('DerBuno', 1000.0 / L1, 1, 'H2 por tonelada, L(1) [kg]')
    m('DerBtanques', 1000.0 / L_tanques, 1, 'H2 por tonelada, O2 de tanques [kg]')
    m('DerReduccion', 100.0 * (1.0 - L_tanques / L1), 1, 'reducción por usar O2 del aire [%]')

    # ── q(eta, phi): componentes a eta = 0, phi objetivo ─────────────────────
    phi_o = p.phi_objetivo
    aire = AIRE_ESTEQ
    t_in = p.T_aire_exterior - p.T_ref_comb
    t_out = p.T_salida_cond - p.T_ref_comb
    entra = (p.cp_H2 + aire * p.cp_aire) * t_in
    sale = ((aire - p.ratio_O2_H2) * p.cp_N2 + r * (phi_o * p.cp_agua_liq + (1 - phi_o) * p.cp_vapor)) * t_out
    m('DerAire', aire, 2, 'aire estequiométrico [kg/kg H2]')
    m('DerWOdos', p.w_O2_aire, 4, 'fracción másica de O2 en el aire')
    m('DerCpHdos', p.cp_H2, 0, 'cp H2')
    m('DerCpAire', p.cp_aire, 0, 'cp aire')
    m('DerCpNdos', p.cp_N2, 0, 'cp N2')
    m('DerCpLiq', p.cp_agua_liq, 0, 'cp agua líquida')
    m('DerCpVapor', p.cp_vapor, 0, 'cp vapor')
    m('DerQPCI', p.Q_comb_H2 / 1e6, 2, 'término (1-eta) PCI, eta = 0 [MJ/kg]')
    m('DerQLatente', phi_o * r * p.L_v_agua / 1e6, 2, 'término phi r L_v [MJ/kg]')
    m('DerQEntra', entra / 1e6, 3, 'calor sensible de los reactivos [MJ/kg]')
    m('DerQSale', sale / 1e6, 3, 'calor sensible del escape que sale [MJ/kg]')
    q_o = p.Q_comb_H2 + phi_o * r * p.L_v_agua + entra - sale
    m('DerQTotal', q_o / 1e6, 2, 'q(0, phi objetivo) [MJ/kg]')
    v['q_rech_modelo'] = q_rech(0.0, phi_o)

    # ── phi_max, punto de rocío, presión ────────────────────────────────────
    ps40 = p_sat(p.T_salida_cond)
    m('DerXOdos', p.X_O2_aire, 4, 'fracción molar de O2 en el aire')
    m('DerNmolar', 0.5 * (1.0 - p.X_O2_aire) / p.X_O2_aire, 3, '0.5 (1 - X_O2)/X_O2')
    m('DerN', N_INERTE, 3, 'mol N2+Ar por mol H2O')
    m('DerPsatCuarenta', ps40 / 1e3, 3, 'p_sat a 40 C [kPa]')
    phi40 = 1.0 - N_INERTE * ps40 / (p.P_atmos - ps40)
    m('DerPhiCuarenta', phi40, 3, 'phi_max a 40 C y 1 atm')
    T_r, x_w = T_rocio()
    m('DerXw', 100.0 * x_w, 1, 'fracción molar de agua [%]')
    m('DerTRocio', T_r, 1, 'punto de rocío a 1 atm [C]')
    P_req = presion_para_phi(phi_o, p.T_salida_cond)
    m('DerPReq', P_req / 1e3, 0, 'presión para phi objetivo a 40 C [kPa]')
    m('DerPReqAtm', P_req / p.P_atmos, 1, 'presión para phi objetivo a 40 C [atm]')
    m('DerTUnaAtm', T_para_phi(phi_o), 1, 'temperatura para phi objetivo a 1 atm [C]')
    # mínimo (M-1, auditoría 3): a 40 C el escape ya está saturado a 1 atm (p_s (1 + N) < P_atm), así que solo se
    # comprime el inerte con el vapor en p_s: W = N R T ln((P - p_s) / (P_atm - p_s)); lo condensado sale líquido
    assert ps40 * (1.0 + N_INERTE) < p.P_atmos
    W = N_INERTE / p.M_H2 * p.R_univ * (p.T_salida_cond + 273.15) * math.log((P_req - ps40) / (p.P_atmos - ps40))
    m('DerWComp', W / 1e6, 2, 'trabajo isotérmico mínimo de compresión, solo la fase gaseosa [MJ/kg H2]')
    m('DerMolesGas', (1.0 + N_INERTE) / p.M_H2, 0, 'mol de escape (H2O + N2 + Ar) por kg de H2')

    # ── Masa e inercia del caso base ────────────────────────────────────────
    d = p.PLATAFORMAS[PLAT]
    eq = mm.equilibrio(PLAT, M_CARGA)
    m_H2_eq = eq['m_gas']
    M_tot = eq['M_total']
    M_I = inercia.inercial_casco(d['V_casco'], d['k'])
    m('DerVcasco', d['V_casco'], 0, 'volumen total del casco, caso base [m3]')
    m('DerMaireCasco', p.rho_aire_std * d['V_casco'] / 1e3, 1, 'aire del volumen del casco, rho V_casco [t]')
    m('DerMvacio', d['m_vacio'] / 1e3, 1, 'peso vacío del caso base [t]')
    m('DerMHdosEq', m_H2_eq / 1e3, 2, 'H2 de equilibrio del caso base con 1 t [t]')
    m('DerVgas', eq['V_gas'], 0, 'volumen de gas en equilibrio con 1 t [m3]')
    m('DerLlenado', eq['llenado'], 2, 'fracción del casco con gas, con 1 t')
    m('DerMtotal', M_tot / 1e3, 2, 'M_total del caso base con 1 t [t]')
    m('DerK', d['k'], 2, 'k del caso base (Lamb, su esbeltez)')
    m('DerMI', M_I / 1e3, 1, 'M_I = (1+k) rho V_casco [t]')
    m('DerMIpost', (M_I - M_CARGA) / 1e3, 1, 'M_I - m_load [t]')

    # ── Precondicionamiento (antes de soltar) ───────────────────────────────
    D = p.D_lim_oper
    K_pre = math.sqrt(p.g) * M_CARGA ** 1.5 / math.sqrt(6.0 * M_I * L1 ** 2 * D)
    m('DerKpre', K_pre, 2, 'K_pre [kg/s]')

    # ── Quema limpia (hallazgo C1, 2026-09-27; bitacora/MINIMAX_ANALITICO_2026-09-27.txt) ──
    # Una sola quema a ṁ* que cruza la liberación: antes de soltar hunde D_lim y da v(0) = -3 D/t_ant;
    # después la nave (inercia M_I - m) llega al reposo justo con el H2 de equilibrio m/L.
    # Impulso: M_I |v(0)| = L g ṁ t_ant²/2 = L g ṁ t_r²/2 · M_I/(M_I - m)  ->  t_r = s t_ant, s = √(1 - m/M_I)
    # Masa: ṁ (t_ant + t_r) = m/L  ->  α* = t_ant/(t_ant + t_r) = 1/(1 + s),  ṁ* = K_pre α*^{3/2}
    s_ = math.sqrt(1.0 - M_CARGA / M_I)
    alpha_ex = 1.0 / (1.0 + s_)
    dm_star = K_pre * alpha_ex ** 1.5
    t_ant = (6.0 * M_I * D / (L1 * p.g * dm_star)) ** (1 / 3)
    t_r = s_ * t_ant
    m('DerS', s_, 5, 's = sqrt(1 - m/M_I)')
    m('DerAlphaAprox', 0.5, 3, 'alpha* de la ley general (s = 1)')
    m('DerAlphaExacta', alpha_ex, 4, 'alpha* = 1/(1+s)')
    m('DerDmStar', dm_star, 2, 'ṁ* = K_pre alpha*^1.5 [kg/s]')
    m('DerDmLey', K_pre * 0.5 ** 1.5, 2, 'ṁ* de la ley general, K_pre/2^1.5 [kg/s]')
    m('DerTAnticip', t_ant, 2, 't_anticip [s]')
    m('DerTReactiva', t_r, 2, 'quema después de soltar, s t_anticip [s]')
    m('DerVcero', 3.0 * D / t_ant, 4, 'velocidad hacia abajo al soltar, 3 D/t_ant [m/s]')
    m('DerHdosQuemado', dm_star * (t_ant + t_r), 1, 'H2 quemado ṁ*(t_ant + t_r) = m/L [kg]')
    m('DerHundTotal', (1.0 + s_) * D, 3, 'bajada total (1+s) D_lim [m]')
    # bajada después de soltar (revisión 2026-09-30, R-07): fuerza L g ṁ (t_r - t) sobre M_I - m desde -|v0|;
    # integrando dos veces hasta t_r, |v0| t_r / 3 (= s D_lim)
    m('DerHundReactiva', 3.0 * D / t_ant * t_r / 3.0, 4, 'bajada después de soltar, |v0| t_r / 3 [m]')
    # R-03: la tasa con phi = 1 frente a phi objetivo (ṁ* ∝ 1/L)
    m('DerTasaPhiUnoPct', 100.0 * (1.0 - (k_rho - 1.0 + r * phi_o) / L1), 2,
      'cuánto baja la tasa al dimensionar con phi = 1 en vez de phi objetivo, 1 - L(phi)/L(1) [%]')
    v['resolver'] = mm.resolver(d, M_CARGA, p.h_max_precision, D)

    # ── Cota con reinyección (quema máxima hasta el pico, pico = h_max) ─────
    h = p.h_max_precision
    mu = M_I / (M_I - M_CARGA)
    w = 2.0 * math.cosh(math.acosh(math.sqrt(mu) * (1.0 + h / D)) / 3.0)
    u = w / math.sqrt(mu)
    G = 6.0 + 3.0 * (w + 1.0 / w) / math.sqrt(mu)
    dm_cota = 6.0 / G ** 1.5 * math.sqrt(p.g * M_CARGA ** 3 / (L1 ** 2 * M_I * D))
    m('DerCotaW', w, 4, 'w = 2 cosh(arccosh(sqrt(mu)(1 + h/D))/3)')
    m('DerCotaU', u, 3, 'u = t_pico/t_ant')
    m('DerCotaG', G, 3, 'G = g m t_ant^2/(M_I D)')
    m('DerCotaDm', dm_cota, 2, 'ṁ de la cota con reinyección, 1 m [kg/s]')
    m('DerCotaSobrequema', 100.0 * (6.0 * (1.0 + u) / G - 1.0), 0, 'H2 quemado de más al pico [%]')

    # ── Controlador ─────────────────────────────────────────────────────────
    # (autor, 2026-10-09) K_d y omega_n con la rigidez real (rho_a - rho_H2) g; ANTES con rho_a g y zeta 1.35
    rho_h = rho_a / k_rho
    Kd = 2.0 * p.zeta * math.sqrt(p.Kp * M_I / ((rho_a - rho_h) * p.g))
    omega = math.sqrt((rho_a - rho_h) * p.g * p.Kp / M_I)
    m('DerKd', Kd, 1, 'K_d [m3 s/m]')
    m('DerOmega', omega, 4, 'omega_n [rad/s]')
    # periodo natural con la rigidez real (rho_aire - rho_H2) g K_p (R-11, revisión 2026-09-30)
    m('DerPeriodo', 2.0 * math.pi / omega / 60.0, 1,
      'periodo natural del lazo de altitud con la rigidez real [min]')

    # ── Calor por liberación, potencia, tiempo, aire de enfriamiento ────────
    L_o = k_rho - 1.0 + r * phi_o
    Q_t = 1000.0 / L_o * q_o
    m('DerLobjetivo', L_o, 2, 'L(phi objetivo)')
    m('DerQPorT', Q_t / 1e9, 2, 'Q por tonelada [GJ]')
    m('DerTDiezMW', Q_t / 10e6 / 60.0, 1, 't = Q/P con 10 MW [min]')
    m('DerPrej', v['resolver']['dm_star'] * q_o / 1e6, 0, 'P_rej = ṁ* q, caso base 1 t, 1 m [MW]')
    aire_enfr = q_o / (p.cp_aire * p.dT_aire_enfr)
    m('DerAireEnfr', aire_enfr, 0, 'aire de enfriamiento [kg por kg H2]')
    m('DerRazonAire', aire_enfr / AIRE_ESTEQ, 0, 'aire de enfriamiento / aire de combustión')
    m('DerAguaBuffer', Q_t / (p.cp_agua_liq * p.dT_buffer_agua) / 1e3, 1, 'agua para guardar el calor [t]')
    absorbe = 199.18e3 + p.cp_N2 * (p.T_aire_exterior + 273.15 - 77.355)
    m('DerLNdosAbsorbe', absorbe / 1e3, 0, 'calor que absorbe 1 kg de LN2 [kJ/kg]')
    m('DerLNdos', Q_t / absorbe / 1e3, 1, 'LN2 por tonelada [t]')

    # ── eta mínima ──────────────────────────────────────────────────────────
    q0_85 = q_rech(0.0, phi40)
    m('DerQCeroAtm', q0_85 / 1e6, 1, 'q(0, phi a 1 atm) [MJ/kg]')
    for f, suf in zip(p.perdida_rad, ('Hoy', 'Futuro')):
        m('DerEtaMin' + suf, f * q0_85 / ((1.0 + f) * p.Q_comb_H2), 3, f'eta_min a 1 atm, f = {f}')

    # ── Temperatura de llama (escape que entra a la unidad de rechazo si eta = 0) ──
    import janaf as jn
    M_aire = sum(x * jn.M(e) for e, x in p.AIRE_X.items())      # moles de cada especie en 1 kg de aire seco
    aire = {e: x / M_aire for e, x in p.AIRE_X.items()}
    k_aire = 0.5 / aire['O2']                          # aire justo por mol de H2
    reac = {e: v * k_aire for e, v in aire.items()}
    T0 = p.T_aire_exterior + 273.15
    H_reac = jn.h('H2', T0) + jn.h_mez(reac, T0)
    prod = dict(reac)
    prod['O2'] -= 0.5
    prod['H2O'] = prod.get('H2O', 0.0) + 1.0
    T_ad = jn.T_de_h(prod, H_reac, T_ini=2000.0)
    m('DerTad', T_ad - 273.15, 0, 'temperatura de llama adiabática, aire justo a T exterior [°C]')

    # ── Casco como intercambiador pasivo: W/m2 y por nave ────────────────────
    w_m2 = [hc * p.dT_casco for hc in p.h_conv_casco]
    rango_m = lambda nombre, a, b, dec, com: macros.append((nombre, (a, b), dec, com)) or v.update({nombre: (a, b)})
    rango_m('DerCascoWm', w_m2[0], w_m2[1], 0, 'casco como intercambiador pasivo [W/m2]')
    for n_, d_ in p.PLATAFORMAS.items():
        rango_m('DerCasco' + d_['corto'], d_['A_total'] * w_m2[0] / 1e6, d_['A_total'] * w_m2[1] / 1e6, 1,
                f'casco como intercambiador pasivo [MW] — {n_}')

    # ── Escalón máximo (aprox. m^{3/2}; la quema limpia no depende de h_max) ──
    P10 = 10e6
    m('DerEscalonAprox', 1000.0 * (P10 / (v['resolver']['dm_star'] * q_o)) ** (2 / 3), 0,
      'escalón con 10 MW por la ley m^{3/2} [kg]')

    # ── Sensibilidad a la inercia (ṁ ∝ M_I^{-1/2}) ──────────────────────────
    for ratio_m, suf in ((0.77, 'Restringido'), (1.20, 'Libre')):
        m('DerInercia' + suf, 100.0 * (ratio_m ** -0.5 - 1.0), 1, f'cambio de flujo si M_I x {ratio_m} [%]')

    # ── Versión 8 del paper (2026-09-27): O2, 1 atm, banda fría, masas, balance neto, helio ──
    m('DerOdos', p.ratio_O2_H2, 2, 'O2 por kg de H2 = r - 1 (lastre si se guardara) [kg/kg]')
    B_atm = 1000.0 / (k_rho - 1.0 + r * phi40)                 # H2 por tonelada con phi_max(40 C, 1 atm)
    m('DerBatm', B_atm, 1, 'H2 por tonelada a 1 atm, phi_max(40 C) [kg/t]')
    m('DerRatm', B_atm * phi40, 1, 'H2 regenerable por tonelada a 1 atm = B phi [kg/t]')
    # banda fría (EXP-14): kW/kg del radiador proporcional a (T - T_aire), calibrado con FlyZero
    import exp14_calor_por_temperatura as e14m
    kf = [kw / (T - p.T_aire_exterior) for kw, T in zip((10.0, 20.0), e14m.T_RECHAZO_FLYZERO)]
    m('DerKFriaHoy', kf[0], 3, 'k de la banda fría, 2026 [kW/kg por K]')
    m('DerKFriaFuturo', kf[1], 3, 'k de la banda fría, 2035 [kW/kg por K]')
    phi14 = e14m.phi_T(p.T_salida_cond, p.P_atmos)
    B14 = 1000.0 / (k_rho - 1.0 + r * phi14)
    Q_frio = B14 * e14m.banda_fria(p.P_atmos, phi14)['Q_banda']
    m('DerQFrioT', Q_frio / 1e9, 2, 'calor de la banda fría por tonelada, 1 atm [GJ/t]')
    m('DerTina', Q_frio / (p.cp_agua_liq * 10.0) / 1e3, 1, 'agua para guardar la banda fría subiendo 10 K [t/t]')
    # caja de conversión con la quema de referencia (10 MW de calor, eta de referencia)
    q_ref = q_rech(p.rend_gen, phi40)
    m('DerQRef', q_ref / 1e6, 1, 'q(eta ref, phi 1 atm) [MJ/kg]')
    m('DerQuemaRef', p.P_ref_masa / q_ref, 3, 'quema que produce 10 MW de calor [kg/s]')
    # regeneración, compresión y tanques (arquitectura (b): todo el H2 regenerado a los tanques)
    m('DerPpem', B_atm * phi40 * p.E_H2_electrolisis / 24.0, 0, 'electrolizador para 1 t en 24 h [kW]')
    m('DerEcomp', p.W_iso_20_350 / p.eta_compresion_H2, 2, 'electricidad de compresión típica [kWh/kg]')
    m('DerCompFlujo', B_atm * phi40 / p.t_regeneracion_h, 2, 'compresor al ritmo del electrolizador [kg/h]')
    m('DerTanquePorKg', (1.0 - p.wt_COPV) / p.wt_COPV, 1, 'kg de tanque por kg de H2 guardado')
    # helio: trabajo isotérmico de gas ideal (contraste del de gas real de EXP-12)
    m('DerWHeIdeal', p.R_univ / p.M_He * p.T1 * math.log(20e6 / p.P_atmos) / 1e6, 2,
      'trabajo isotérmico de gas ideal, He de 1 atm a 20 MPa [MJ/kg]')

    # ── Casco como intercambiador pasivo (sec. 5.3): W/m2 = h dT, por la superficie del casco ──
    m('DerHConvCalma', p.h_conv_casco[0], 0, 'convección exterior sobre el casco, calma [W/m2K]')
    m('DerHConvViento', p.h_conv_casco[1], 0, 'convección exterior sobre el casco, con viento [W/m2K]')
    m('DerDTCasco', p.dT_casco, 0, 'diferencia media pared-aire [K]')
    m('DerAtotalBase', d['A_total'], 0, 'superficie del casco del caso base [m2]')

    # ── Regeneración con paneles (sec. 6): días por tonelada = R E_el / (G_ef ef A_pv), caso base, a-Si, invierno ──
    ir = rutas.leer_datos('irradiancia_casco')['estaciones']['invierno']
    A_pv = p.A_paneles_fraccion * d['A_total']
    G_ef = ir['GHI'] * ir['F']                                  # sobre la mitad superior del casco
    ef_asi, kgm2_asi = p.PV_TEC['a-Si']
    E_el_t = B_atm * phi40 * p.E_H2_electrolisis              # kWh por tonelada
    m('DerApvBase', A_pv, 0, 'área de paneles del caso base [m2]')
    m('DerGefInv', G_ef, 2, 'irradiación sobre la mitad superior del casco, invierno [kWh/m2/día]')
    m('DerEdiaBase', A_pv * G_ef * ef_asi, 0, 'electricidad diaria de los paneles a-Si del caso base, invierno [kWh]')
    m('DerElecT', E_el_t, 0, 'electrólisis por tonelada a 1 atm [kWh]')
    m('DerDiasBase', E_el_t / (A_pv * G_ef * ef_asi), 2, 'días por tonelada, caso base, a-Si, invierno')
    m('DerMasaPvBase', A_pv * kgm2_asi / 1e3, 1, 'masa de los paneles a-Si del caso base [t]')

    # ── Operación en altura (EXP-15): día caluroso con el gradiente ISA, 1 km ──
    T0 = p.T_aire_exterior + 273.15
    T_h = T0 - p.gradiente_ISA * 1000.0
    expo = p.g / (p.R_aire * p.gradiente_ISA)
    P_h = p.P_atmos * (T_h / T0) ** expo
    rho_rel = (P_h / p.P_atmos) * (T0 / T_h)
    T_sal_h = T_h - 273.15 + (p.T_salida_cond - p.T_aire_exterior)
    ps_h = p_sat(T_sal_h)
    m('DerTAireAlt', T_h - 273.15, 1, 'aire a 1 km, día caluroso + gradiente ISA [C]')
    m('DerExpAlt', expo, 2, 'exponente g/(R_aire Gamma)')
    m('DerPAlt', P_h / 1e3, 1, 'presión a 1 km [kPa]')
    m('DerRhoRelAlt', rho_rel, 3, 'densidad del aire a 1 km / nivel del mar')
    m('DerTSalAlt', T_sal_h, 1, 'salida del condensador a 1 km [C]')
    m('DerPhiAlt', 1.0 - N_INERTE * ps_h / (P_h - ps_h), 3, 'phi_max a 1 km')

    # ── Condensar a presión (discusión): la banda fría a 5 atm ──
    P5 = 5.0 * p.P_atmos
    T_r5 = e14m.T_sat(P5 / (1.0 + N_INERTE))
    T_f5 = e14m.T_sat((1.0 - phi14) / (1.0 - phi14 + N_INERTE) * P5)
    m('DerTRocioCinco', T_r5, 0, 'rocío a 5 atm [C]')
    m('DerTFinCinco', T_f5, 0, 'fin de la condensación a 5 atm, misma phi [C]')
    # mínimo, solo la fase gaseosa (M-1): a T_f5 el escape no está saturado a 1 atm (T_f5 > rocío a 1 atm), así que
    # de 1 atm a la presión de rocío P_r = p_s (1 + N) se comprimen los 1 + N mol, y de P_r a 5 atm solo el inerte
    ps_f5 = p_sat(T_f5)
    P_r5 = ps_f5 * (1.0 + N_INERTE)
    assert p.P_atmos < P_r5 < P5
    RT5 = p.R_univ * (T_f5 + 273.15) / p.M_H2
    W5 = (1.0 + N_INERTE) * RT5 * math.log(P_r5 / p.P_atmos) + N_INERTE * RT5 * math.log((P5 - ps_f5) / (P_r5 - ps_f5))
    m('DerWCinco', W5 / 1e6, 1, 'trabajo isotérmico mínimo para comprimir el escape a 5 atm, solo la fase gaseosa [MJ/kg H2]')
    v['phi14'] = phi14

    # ── Gas isotérmico (discusión): sustentación por K de sobrecalentamiento del gas ──
    m('DerKgPorK', M_tot / p.T1, 0, 'cambio de sustentación por K del gas, caso base con 1 t [kg/K]')

    rutas.guardar_datos('derivaciones', {'valores': v, 'macros': macros})
    for nombre, valor, dec, com in macros:
        txt = f'{valor[0]:.{dec}f}--{valor[1]:.{dec}f}' if isinstance(valor, (list, tuple)) else f'{valor:.{dec}f}'
        print(f'{nombre:22s} {txt:>14s}  {com}')


if __name__ == '__main__':
    main()
