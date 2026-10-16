# -*- coding: utf-8 -*-
# [EN] Vertical inertia with added mass (Equation 2, Section 3.1).
# =============================================================================
# inercia.py — masa inercial VERTICAL del dirigible: una sola fuente para todas las
# simulaciones y formas cerradas (decisión del autor 2026-09-24: la masa añadida entra
# a la dinámica, no como comentario).
#
# Al acelerar en vertical, la nave mueve más que su propia masa:
#   - su masa (estructura, carga, gas de sustentación, agua retenida);
#   - el aire que va dentro del casco (balonets y espacio entre celdas): suma inercia,
#     no peso (su peso lo compensa el aire que desplaza);
#   - la masa añadida del aire exterior, k·rho_aire·V_casco, con k el coeficiente
#     transversal de un esferoide prolato de la esbeltez de la nave (Lamb, Hydrodynamics
#     §373; Tuckerman, NACA Report 210, 1926; parametros.k_lamb). Cada nave usa SU k.
# V_casco es el volumen TOTAL del casco, un dato de la nave (parametros.PLATAFORMAS), y no
# cambia en la maniobra. La masa añadida no entra en el peso ni en la flotación: solo en la
# inercia, y por eso no cambia ningún resultado estático.
#
# En CUALQUIER equilibrio (fuerza neta cero: M = rho_aire·V_gas) la inercia es
#     M_I = M + rho_aire·(V_casco - V_gas) + k·rho_aire·V_casco = (1 + k)·rho_aire·V_casco,
# pese lo que pese la nave: lo que no es gas es aire interior. Quemar B kg de H2 suma L·B
# (sale gas, entra aire, queda el agua) y soltar la carga resta m_load; en el nuevo
# equilibrio (B = m_load/L) vuelve al mismo valor (autor, 2026-09-28; hasta entonces el
# casco se suponía lleno de gas y se escribía (1 + k)·M_total, que es lo mismo solo así).
# =============================================================================

import parametros as p


def inercial(M_nave, rho_aire, V_casco, V_gas, k):
    """Masa inercial vertical [kg] en la simulación.
    M_nave: masa propia (estructura, carga, gas, agua); V_gas: volumen actual de gas en el casco;
    k: masa añadida de la nave (parametros.PLATAFORMAS[...]['k'])."""
    return M_nave + rho_aire * max(V_casco - V_gas, 0.0) + k * rho_aire * V_casco


def inercial_casco(V_casco, k, rho_aire=None):
    """Forma cerrada: inercia en equilibrio, (1 + k)·rho_aire·V_casco [kg]."""
    rho_aire = p.rho_aire_std if rho_aire is None else rho_aire
    return (1.0 + k) * rho_aire * V_casco
