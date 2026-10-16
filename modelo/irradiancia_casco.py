# -*- coding: utf-8 -*-
# [EN] Solar irradiance on the hull panels, measured in Lima (Section 6).
# =============================================================================
# irradiancia_casco.py — irradiancia sobre los paneles del casco (regeneración solar, paper 1).
#
# Datos medidos: SENAMHI 2003, estación A. von Humboldt (Lima): transmisividad H/H0 de cada
# mes (anexo 6) y energía solar diaria anual (tabla 2). H = K_T H0 con H0 del día
# representativo de cada mes (Klein 1977).
#
# Método (cielo isotrópico de Liu & Jordan; Duffie & Beckman, caps. 1-2):
#   - fracción difusa de las medias mensuales: Page (1964, p. 383), H_d/H = 1.00 - 1.13 K_T;
#   - los paneles cubren la mitad superior del casco (ángulo central -90 a +90 grados, como en
#     Sun et al. 2024) sin las puntas; ese cuerpo central se toma como un cilindro horizontal.
#     Para un elemento de normal n, lo que recibe de la directa es max(0, n·s). Su promedio
#     sobre la mitad superior es (|s_p| + s_z)/pi, con s_p la componente del sol en el plano de
#     la sección (deducción: integrar cos(theta - psi) donde es positivo). Se integra en el día
#     con el peso de la radiación extraterrestre (método de Liu & Jordan / Klein) -> R_b, y se
#     promedia sobre el rumbo del eje (la nave se orienta con el viento); N-S y E-W dan el rango;
#   - difusa isotrópica: el promedio de (1 + cos b)/2 sobre la mitad superior es 1/2 + 1/pi;
#   - reflejada por el suelo: albedo x (1/2 - 1/pi).
#   F = (1 - f_d) R_b + f_d (1/2 + 1/pi) + albedo (1/2 - 1/pi);  G_ef = F H  por m2 de panel.
# Salida: tabla mensual; invierno (jun-ago), verano (dic-feb) y anual (promedios por día).
# =============================================================================

import math
import numpy as np
import parametros as p
import rutas

PAGE_B = 1.13                     # Page (1964, p. 383, ec. [6]): H_d/H = 1.00 - PAGE_B K_T (medias mensuales)
F_DIFUSA = 0.5 + 1.0 / math.pi    # difusa isotrópica sobre la mitad superior de un cilindro
F_SUELO = 0.5 - 1.0 / math.pi     # reflejada por el suelo sobre la misma mitad
DIA_REPR = (17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344)   # Klein 1977
DIAS_MES = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
MESES = ('ene', 'feb', 'mar', 'abr', 'may', 'jun', 'jul', 'ago', 'set', 'oct', 'nov', 'dic')
ESTACIONES = {'invierno': (5, 6, 7), 'verano': (11, 0, 1), 'anual': tuple(range(12))}
RUMBOS = np.radians(np.arange(0.0, 180.0, 5.0))   # rumbo del eje del casco (0 = N-S, 90 = E-W)


def declinacion(n):
    """Declinación solar [rad] (Cooper 1969)."""
    return math.radians(23.45) * math.sin(2.0 * math.pi * (284 + n) / 365.0)


def H0_dia(lat_deg, n, G_sc=None):
    """Irradiación extraterrestre diaria sobre un plano horizontal [J/m2]."""
    G_sc = p.G_SC if G_sc is None else G_sc
    phi, d = math.radians(lat_deg), declinacion(n)
    ws = math.acos(-math.tan(phi) * math.tan(d))
    return (24 * 3600 * G_sc / math.pi * (1 + 0.033 * math.cos(2 * math.pi * n / 365.0))
            * (math.cos(phi) * math.cos(d) * math.sin(ws) + ws * math.sin(phi) * math.sin(d)))


def directa_media_semicilindro(s, eje):
    """Promedio de max(0, n·s) sobre la mitad superior de un cilindro horizontal de eje unitario 'eje'
    (forma cerrada; exige s_z >= 0)."""
    s_a = float(np.dot(s, eje))
    return (math.sqrt(max(0.0, 1.0 - s_a * s_a)) + s[2]) / math.pi


def Rb_dia(lat_deg, n, rumbo):
    """Razón diaria directa sobre la mitad superior del casco / directa sobre el horizontal,
    con el peso de la extraterrestre (Liu & Jordan). rumbo: ángulo del eje desde el N [rad]."""
    phi, d = math.radians(lat_deg), declinacion(n)
    ws = math.acos(-math.tan(phi) * math.tan(d))
    w = np.linspace(-ws, ws, 2001)
    sz = math.sin(phi) * math.sin(d) + math.cos(phi) * math.cos(d) * np.cos(w)
    se = -math.cos(d) * np.sin(w)
    sn = math.cos(phi) * math.sin(d) - math.sin(phi) * math.cos(d) * np.cos(w)
    sa = sn * math.cos(rumbo) + se * math.sin(rumbo)
    curvo = (np.sqrt(np.clip(1.0 - sa * sa, 0.0, 1.0)) + sz) / math.pi
    return float(np.trapezoid(curvo, w) / np.trapezoid(sz, w))


def factor_forma(fd, Rb):
    return (1.0 - fd) * Rb + fd * F_DIFUSA + p.albedo_suelo * F_SUELO


def main():
    lat = p.SENAMHI_LAT
    meses = []
    for i, (kt, n) in enumerate(zip(p.SENAMHI_KT_MES, DIA_REPR)):
        H = kt * H0_dia(lat, n) / 3.6e6                        # [kWh/m2/día]
        fd = 1.0 - PAGE_B * kt
        rb = [Rb_dia(lat, n, r) for r in RUMBOS]
        rb_ns, rb_ew, rb_m = Rb_dia(lat, n, 0.0), Rb_dia(lat, n, math.pi / 2), float(np.mean(rb))
        F = factor_forma(fd, rb_m)
        meses.append({'mes': MESES[i], 'dia': n, 'KT': kt, 'H0_kWh': H0_dia(lat, n) / 3.6e6, 'GHI': H, 'fd': fd,
                      'Rb_NS': rb_ns, 'Rb_EW': rb_ew, 'Rb': rb_m, 'F': F,
                      'F_min': factor_forma(fd, min(rb)), 'F_max': factor_forma(fd, max(rb)),
                      'F_fd0': factor_forma(0.0, rb_m), 'F_fd1': factor_forma(1.0, rb_m), 'efectiva': F * H})
    est = {}
    for nom, idx in ESTACIONES.items():
        w = [DIAS_MES[i] for i in idx]
        prom = lambda k: sum(meses[i][k] * wi for i, wi in zip(idx, w)) / sum(w)
        ghi, ef = prom('GHI'), prom('efectiva')
        est[nom] = {'GHI': ghi, 'efectiva': ef, 'F': ef / ghi, 'fd': prom('fd'),
                    'F_min': min(meses[i]['F_min'] for i in idx), 'F_max': max(meses[i]['F_max'] for i in idx),
                    'F_fd0': min(meses[i]['F_fd0'] for i in idx), 'F_fd1': max(meses[i]['F_fd1'] for i in idx)}
    rutas.guardar_datos('irradiancia_casco', {
        'fuente': 'SENAMHI 2003, estación A. von Humboldt (Lima)', 'lat': lat, 'anual_tabla': p.SENAMHI_GHI_ANUAL,
        'anual_reconstruido': est['anual']['GHI'], 'F_difusa': F_DIFUSA, 'F_suelo': F_SUELO,
        'albedo': p.albedo_suelo, 'page_b': PAGE_B, 'meses': meses, 'estaciones': est})
    print(f"GHI anual rearmado {est['anual']['GHI']:.2f} kWh/m2/día (tabla 2: {p.SENAMHI_GHI_ANUAL})")
    print('mes   KT    GHI   f_d   Rb(NS/EW/prom)      F     G_ef')
    for m in meses:
        print(f"{m['mes']}  {m['KT']:.3f} {m['GHI']:5.2f} {m['fd']:5.2f}  {m['Rb_NS']:.3f}/{m['Rb_EW']:.3f}/{m['Rb']:.3f}"
              f"  {m['F']:.3f}  {m['efectiva']:5.2f}")
    for nom, e in est.items():
        print(f"{nom:<9} GHI {e['GHI']:.2f}  F {e['F']:.3f} (rumbo {e['F_min']:.3f}-{e['F_max']:.3f}; "
              f"f_d 0-1: {e['F_fd0']:.3f}-{e['F_fd1']:.3f})  sobre el casco {e['efectiva']:.2f} kWh/m2/día")


if __name__ == '__main__':
    main()
