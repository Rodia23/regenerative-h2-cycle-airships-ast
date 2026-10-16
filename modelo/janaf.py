# -*- coding: utf-8 -*-
# [EN] Ideal-gas properties from NIST-JANAF (adiabatic flame temperature, Section 9).
# =============================================================================
# janaf.py — propiedades de gas ideal de N2, O2, Ar, H2O, CO2 y H2 (temperatura de llama).
# DATOS COPIADOS de las tablas NIST-JANAF (Chase, M. W., NIST-JANAF Thermochemical
# Tables, 4th ed., J. Phys. Chem. Ref. Data Monograph 9, 1998), texto plano de
# https://janaf.nist.gov/tables/<id>.txt, consulta 2026-09-24. Columnas por fila:
#   T [K], Cp [J/mol K], S° [J/mol K] (p° = 0.1 MPa), H - H(298.15) [kJ/mol].
# Entre nodos: h(T) y S°(T) por Hermite cúbico con pendientes Cp y Cp/T, que
# reproduce EXACTO los valores tabulados en cada nodo. La verificación (grupo K)
# compara además con los ajustes de Shomate de NIST WebBook (misma fuente).
# Los números están en janaf_datos.py, GENERADO por generar_janaf.py; no editarlos a mano.
# =============================================================================

import bisect
from janaf_datos import DATOS   # generado por generar_janaf.py (tablas en datos_fuente/janaf/)

R = 8.314462618        # [J/mol K] (CODATA 2018, el mismo de parametros.py); S° es a p° = 0.1 MPa

_T = {e: [f[0] for f in d[3]] for e, d in DATOS.items()}


def _tramo(esp, T):
    t = _T[esp]
    if not t[0] <= T <= t[-1]:
        raise ValueError(f'{esp}: T = {T:.1f} K fuera de la tabla ({t[0]:g}-{t[-1]:g} K)')
    i = min(max(bisect.bisect_right(t, T) - 1, 0), len(t) - 2)
    return DATOS[esp][3][i], DATOS[esp][3][i + 1]


def _hermite(T, T0, T1, y0, y1, d0, d1):
    """Valor y derivada del Hermite cúbico con valores y pendientes en T0 y T1."""
    dx = T1 - T0
    u = (T - T0) / dx
    h00, h10, h01, h11 = 2*u**3 - 3*u**2 + 1, u**3 - 2*u**2 + u, -2*u**3 + 3*u**2, u**3 - u**2
    y = h00 * y0 + h10 * dx * d0 + h01 * y1 + h11 * dx * d1
    dy = ((6*u**2 - 6*u) * y0 + (3*u**2 - 4*u + 1) * dx * d0 + (-6*u**2 + 6*u) * y1 + (3*u**2 - 2*u) * dx * d1) / dx
    return y, dy


def h(esp, T):
    """Entalpía molar absoluta [J/mol] = dfH°(298.15) + [H - H(298.15)](T)."""
    a, b = _tramo(esp, T)
    y, _ = _hermite(T, a[0], b[0], a[3] * 1e3, b[3] * 1e3, a[1], b[1])
    return DATOS[esp][2] * 1e3 + y


def cp(esp, T):
    """Cp molar [J/mol K] = dh/dT del mismo interpolante (coherente con h)."""
    a, b = _tramo(esp, T)
    return _hermite(T, a[0], b[0], a[3] * 1e3, b[3] * 1e3, a[1], b[1])[1]


def s0(esp, T):
    """Entropía molar estándar S°(T) [J/mol K] a p° = 0.1 MPa."""
    a, b = _tramo(esp, T)
    return _hermite(T, a[0], b[0], a[2], b[2], a[1] / a[0], b[1] / b[0])[0]


def M(esp):
    return DATOS[esp][1]


# ── Mezclas de composición fija (moles por unidad de base) ─────────────────
def h_mez(n, T):
    return sum(ni * h(e, T) for e, ni in n.items() if ni)


def s0_mez(n, T):
    return sum(ni * s0(e, T) for e, ni in n.items() if ni)


def masa(n):
    return sum(ni * M(e) for e, ni in n.items() if ni)


def moles(n):
    return sum(n.values())


def T_de_h(n, H, T_ini=800.0):
    """Temperatura con entalpía H [J] (Newton con cp de la mezcla)."""
    T = T_ini
    for _ in range(60):
        f = h_mez(n, T) - H
        c = sum(ni * cp(e, T) for e, ni in n.items() if ni)
        dT = -f / c
        T = min(max(T + dT, _T['N2'][0]), _T['N2'][-1])
        if abs(dT) < 1e-9:
            return T
    raise RuntimeError('T_de_h no converge')
