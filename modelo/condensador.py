# -*- coding: utf-8 -*-
# [EN] Heat-rejection unit: installed power, first-order lag and condensed fraction (Section 3.4).
# =============================================================================
# condensador.py — dinámica mínima del condensador para las simulaciones del
# artículo (EXP-09). Un solo lugar; lo usa drop_flujo_optimo_perdrop.correr().
#
# El combustor sigue siendo una caja negra (quema H2 con aire justo). Por cada kg
# de H2 quemado salen 8.94 kg de vapor (M_H2O/M_H2); como mucho condensa la fracción de equilibrio
# phi a la temperatura de salida (dilución con N2, EXP-05); el resto sale siempre
# con el N2. Condensar exige rechazar h_agua = q / (r phi) por kg de agua, con
# q = q_rech(0, phi) de EXP-05 (balance de energía, sin exportar trabajo). El condensador:
#   - tiene una capacidad P_cap [W] (None = sin límite);
#   - llega a su potencia con un retardo térmico de primer orden tau [s];
#   - el vapor condensable que no alcanza a condensar:
#       modo 'ventear': sale con el N2 (se pierde el lastre);
#       modo 'pulmon':  queda como gas dentro de la envolvente y se condensa
#                       después; mientras tanto ocupa volumen a T_pulmon y 1 atm,
#                       así que empuja hacia arriba: (rho_aire V_vapor - m_vapor) g.
# Sin condensador (cond=None) las simulaciones hacen lo del artículo: r kg de agua
# a bordo al instante por kg de H2 (condensación instantánea: la cota ideal).
# =============================================================================

import math
import parametros as p
from exp05_condensacion_tiempo_real import phi_max, q_rech

R_AGUA = 1.0 + p.ratio_O2_H2        # kg de agua por kg de H2 (M_H2O/M_H2)


def config(P_cap_W=None, tau_s=None, modo='ventear', T_salida=None, prearranque=True):
    """Parámetros del condensador para una corrida.
    prearranque = True (la OPERACIÓN del paper, sec. 2.3): la unidad recibe la orden de arrancar
    n_tau_arranque·tau antes de la quema previa y su CAPACIDAD sube hacia P_cap con el retardo tau;
    condensa lo que esa capacidad permite. La simulación fija el instante de la orden ('t_arranque').
    prearranque = False (solo sensibilidad): la potencia SIGUE a la quema desde cero, con retardo."""
    phi = phi_max(p.T_salida_cond if T_salida is None else T_salida)
    q = q_rech(0.0, phi)
    return {'P_cap': float('inf') if P_cap_W is None else float(P_cap_W),
            'tau': p.tau_condensador if tau_s is None else float(tau_s),
            'modo': modo, 'phi': phi, 'q': q, 'h_agua': q / (R_AGUA * phi),
            'prearranque': bool(prearranque), 't_arranque': None}


def potencia_inicial(c):
    """Capacidad al comienzo de la simulación (t = 0) si la orden de arranque fue antes de t = 0."""
    if not c['prearranque'] or c['t_arranque'] is None or not math.isfinite(c['P_cap']) or c['t_arranque'] > 0.0:
        return 0.0
    if c['tau'] <= 0.0:
        return c['P_cap']
    return c['P_cap'] * (1.0 - math.exp(c['t_arranque'] / c['tau']))


def volumen_pulmon(B_v):
    """Volumen [m3] del vapor retenido como gas (T_pulmon, 1 atm)."""
    return B_v * p.R_vapor * (p.T_pulmon + 273.15) / p.P_atmos


def derivadas(quema, P_c, B_v, c, t=None):
    """quema [kg/s de H2]; P_c potencia (o capacidad, con prearranque) del condensador [W]; B_v vapor en
    el pulmón [kg]; t [s] instante (para la orden de arranque).
    Devuelve (dP_c, agua condensada [kg/s], agua venteada [kg/s], dB_v [kg/s])."""
    condensable = R_AGUA * c['phi'] * quema
    siempre_sale = R_AGUA * (1.0 - c['phi']) * quema
    hay_pulmon = c['modo'] == 'pulmon' and B_v > 1e-9
    if c['prearranque'] and c['t_arranque'] is not None:
        encendido = t is not None and t >= c['t_arranque']
        if not math.isfinite(c['P_cap']):              # capacidad sin límite ya encendida: condensa todo
            dP, tope = 0.0, (float('inf') if encendido else 0.0)
        elif c['tau'] > 0.0:
            dP = ((c['P_cap'] if encendido else 0.0) - P_c) / c['tau']
            tope = max(P_c, 0.0) / c['h_agua']
        else:
            dP, tope = 0.0, (c['P_cap'] if encendido else 0.0) / c['h_agua']
        if c['modo'] == 'pulmon':
            cond = min(tope, condensable + B_v / 1.0) if hay_pulmon else min(condensable, tope)
            return dP, cond, siempre_sale, condensable - cond
        cond = min(condensable, tope)
        return dP, cond, siempre_sale + (condensable - cond), 0.0
    if c['modo'] == 'pulmon':
        demanda = c['P_cap'] if (hay_pulmon or quema > 0.0) else 0.0
    else:
        demanda = min(c['P_cap'], quema * c['q'])
    if not math.isfinite(demanda):                 # sin límite: sigue a la quema (y vacía el pulmón)
        demanda = quema * c['q'] + (B_v * c['h_agua'] / max(c['tau'], 1.0) if hay_pulmon else 0.0)
    if c['tau'] > 0.0:
        dP = (demanda - P_c) / c['tau']
        tope = max(P_c, 0.0) / c['h_agua']
    else:
        dP = 0.0
        tope = demanda / c['h_agua']
    if c['modo'] == 'pulmon':
        # con vapor esperando, condensa a toda la potencia disponible, sin vaciar más de lo que hay
        cond = min(tope, condensable + B_v / 1.0) if hay_pulmon else min(condensable, tope)
        return dP, cond, siempre_sale, condensable - cond
    cond = min(condensable, tope)
    return dP, cond, siempre_sale + (condensable - cond), 0.0
