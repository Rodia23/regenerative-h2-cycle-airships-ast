# -*- coding: utf-8 -*-
# [EN] Closed-form sizing of the clean burn (Equation 7, Appendix A).
# =============================================================================
# minimax.py — forma cerrada del dimensionamiento (Apéndice A) en UN solo lugar.
#
# Maniobra que se dimensiona (quema LIMPIA, hallazgo C1 de la auditoría, 2026-09-27):
# una sola quema a la tasa instalada ṁ*, continua a través de la liberación. Antes
# de soltar hunde la nave D_lim (desplazamiento cúbico) y le da velocidad hacia
# abajo; después de soltar la lleva al reposo justo cuando el H2 quemado es el de
# equilibrio, m_load/L. No se quema de más ni se devuelve H2 de los tanques: la nave
# nunca sube por encima del punto de liberación y termina en reposo (1 + s)·D_lim
# abajo, con s = √(1 − m_load/M_I).
#   α*         = 1/(1 + s)                      (fracción compensada antes de soltar, ≈ 1/2)
#   ṁ*         = K_pre·α*^(3/2),  K_pre = √(g m³ / (6 M_I L² D_lim))
#   t_anticip  = (6 M_I D_lim / (L g ṁ*))^(1/3),  t_reactiva = s·t_anticip
# h_max NO dimensiona esta maniobra (se cumple sola si h_max ≥ (1 + s)·D_lim). La
# ley general del artículo toma s = 1 (misma inercia antes y después de soltar):
# ṁ* = √(g m³ / (48 L² M_I D_lim)), α* = 1/2 para todas las naves.
#
# cota_reinyeccion(): el mínimo ABSOLUTO del flujo instalado (quema máxima hasta el
# pico, pico = h_max). Exige quemar de más y devolver ese H2 desde los tanques; el
# artículo lo menciona como cota, no dimensiona con él.
#
# Derivación completa y verificación: bitacora/MINIMAX_ANALITICO_2026-09-27.txt.
# Antes de 2026-09-27 este módulo tenía un minimax entre una fase de precondicio-
# namiento (D_lim) y una reactiva (h_max) que acreditaba solo rho_aire/rho_H2 y
# suponía v(0) = 0; coincidía con la quema limpia a h_max = 1 m por compensación
# de errores (C ≈ 0.98) y no a 5 m.
# =============================================================================

import numpy as np
import parametros as p
import inercia

rho_aire = p.rho_aire_std
rho_H2   = p.P_atmos / (p.R_H2 * p.T1)
denom    = rho_aire / rho_H2 - 1               # rho_aire/rho_H2 - 1 (equilibrio de flotación)
ratio    = rho_H2 / rho_aire                   # rho_H2/rho_aire
L_eff    = rho_aire / rho_H2 + p.ratio_O2_H2   # rho_aire/rho_H2 + 8


def nave(n):
    """Datos de una nave: acepta el nombre (clave de parametros.PLATAFORMAS) o el dict mismo."""
    return p.PLATAFORMAS[n] if isinstance(n, str) else n


def equilibrio(n, m_load, rho_gas=None):
    """Estado de equilibrio de la nave n con m_load [kg] a bordo, en el estado de referencia:
    gas justo para flotar, V_gas = (m_vacio + m_load)/(rho_aire - rho_gas) <= V_casco; el resto del
    casco es aire. La inercia es (1 + k)·rho_aire·V_casco, pese lo que pese la nave (inercia.py)."""
    d = nave(n)
    rho_gas = rho_H2 if rho_gas is None else rho_gas
    V_gas = (d['m_vacio'] + m_load) / (rho_aire - rho_gas)
    if V_gas > d['V_casco'] * (1.0 + 1e-9):
        raise ValueError(f"{d['corto']}: {m_load:.0f} kg no flotan ({V_gas:.0f} m3 de gas > casco {d['V_casco']:.0f} m3)")
    m_gas = rho_gas * V_gas
    M = d['m_vacio'] + m_load + m_gas
    return {'V_gas': V_gas, 'm_gas': m_gas, 'M_total': M, 'llenado': V_gas / d['V_casco'],
            'M_I': inercia.inercial(M, rho_aire, d['V_casco'], V_gas, d['k'])}


def limpio(M_I, m_load, D_lim, L_gas=None):
    """Quema limpia con la inercia M_I al soltar (con masa añadida); después de soltar acelera M_I - m_load.
    Es la ÚNICA fuente del flujo instalado: la usan resolver() y las simulaciones (drop_flujo_optimo*.py)."""
    L_gas = L_eff if L_gas is None else L_gas
    s         = np.sqrt(1.0 - m_load / M_I)
    K_pre     = np.sqrt(p.g * m_load**3 / (6.0 * M_I * L_gas**2 * D_lim))
    alpha     = 1.0 / (1.0 + s)
    dm_star   = K_pre * alpha ** 1.5
    t_anticip = (6.0 * M_I * D_lim / (L_gas * p.g * dm_star)) ** (1.0 / 3.0)
    return {'s': s, 'K_pre': K_pre, 'alpha': alpha, 'dm_star': dm_star, 't_anticip': t_anticip,
            't_reactiva': s * t_anticip, 'hundimiento_total': (1.0 + s) * D_lim}


def cota_reinyeccion(M_I, m_load, h_max, D_lim, L_gas=None):
    """Mínimo absoluto del flujo instalado: quema máxima desde la quema previa (que hunde D_lim) hasta el
    pico, con el pico en h_max. u = (tiempo de soltar al pico)/t_anticip; sobrequema = H2 quemado al pico
    sobre el de equilibrio, menos 1 (hay que devolverlo desde los tanques)."""
    L_gas = L_eff if L_gas is None else L_gas
    mu = M_I / (M_I - m_load)
    w  = 2.0 * np.cosh(np.arccosh(np.sqrt(mu) * (1.0 + h_max / D_lim)) / 3.0)
    u  = w / np.sqrt(mu)
    G  = 6.0 + 3.0 * (w + 1.0 / w) / np.sqrt(mu)
    dm = 6.0 / G ** 1.5 * np.sqrt(p.g * m_load**3 / (L_gas**2 * M_I * D_lim))
    return {'dm': dm, 'u': u, 'sobrequema': 6.0 * (1.0 + u) / G - 1.0}


def resolver(n, m_load, h_max, D_lim, rho_gas=None, L_gas=None, k_ma=None):
    """Dimensionamiento de una liberación de m_load [kg] desde la nave n (nombre o dict de
    parametros.PLATAFORMAS) con hundimiento previo D_lim [m].

    h_max [m] solo entra en la cota con reinyección (no dimensiona la quema limpia).
    Por defecto el gas es H2 con el agua retenida (L_eff). Para otro gas (EXP-12, helio) se pasan
    rho_gas [kg/m3] y L_gas = kg de pesadez ganados por kg de gas procesado (He comprimido que
    queda a bordo: rho_aire/rho_He). Misma nave (casco, peso vacío y carga); el gas de equilibrio
    cambia, la inercia no (depende del casco).

    La inercia incluye la masa añadida (inercia.py): (1 + k)·rho_aire·V_casco antes de soltar y
    eso menos m_load después. k_ma reemplaza el k de la nave (k_ma = 0: sin masa añadida)."""
    if rho_gas is None and L_gas is not None:
        raise ValueError('L_gas sin rho_gas se ignoraría: pasar los dos (para H2, rho_gas=minimax.rho_H2)')
    L_g = L_eff if rho_gas is None else L_gas
    d = nave(n)
    eq  = equilibrio(d, m_load, rho_gas)
    k   = d['k'] if k_ma is None else k_ma
    M_I = inercia.inercial_casco(d['V_casco'], k, rho_aire)   # inercia antes de soltar
    ex      = limpio(M_I, m_load, D_lim, L_g)
    cota    = cota_reinyeccion(M_I, m_load, h_max, D_lim, L_g)
    return {
        'M_total': eq['M_total'], 'M_inercial': M_I, 'epsilon': m_load / (rho_aire * d['V_casco']),
        'V_gas': eq['V_gas'], 'llenado': eq['llenado'], 'k': k,
        **ex,
        # ley general (s = 1: misma inercia antes y después de soltar)
        'alpha_aprox': 0.5, 'dm_star_aprox': ex['K_pre'] * 0.5 ** 1.5,
        'dm_cota': cota['dm'], 'u_cota': cota['u'], 'sobrequema_cota': cota['sobrequema'],
    }
