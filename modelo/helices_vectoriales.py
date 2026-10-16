# -*- coding: utf-8 -*-
# [EN] Propeller power and response ceiling of the hold-down device (Section 5.2).
# =============================================================================
# helices_vectoriales.py — discusión del paper 1, "Braking" (revisión 2026-10-01). FUERA de la premisa del paper 1
# (ningún propulsor frena la nave): ¿qué pasaría si las hélices vectoriales que la nave ya lleva para maniobrar la
# sostuvieran mientras se quema al ritmo que el rechazo instalado puede absorber? Su masa ya está en la masa vacía;
# la energía, de la planta de potencia de la propia nave. Solo para despejar la duda de un revisor y mostrar qué
# queda por resolver:
#   - la quema dura t = m q / (L P_rej) (la misma de la sec. 5.3);
#   - la fuerza a sostener empieza en m g y baja linealmente a cero mientras la quema quita la sustentación;
#   - hélices ideales (teoría de cantidad de movimiento) con área total de discos A y figura de mérito FM:
#     P = F^(3/2) / (FM raiz(2 rho A)); con F lineal a cero, E = P_max t / 2.5 (integral de (1 - u)^1.5);
#   - la energía, comparada con la química del H2 que se quema por tonelada (m/L · PCI);
#   - la unidad de rechazo sigue haciendo falta y debe condensar a la temperatura del ciclo: su masa con las clases de
#     la tabla 7 (aviones, a 80-160 °C; haces de condensador de central, a la temperatura del ciclo), al rechazo de
#     referencia y a uno menor (más liviana, pero sostiene más tiempo);
#   - la potencia crece como F^1.5: escalones grandes;
#   - cota "todas las hélices verticales" (para cuando haya datos de propulsión por nave): escalón máximo
#     m_max = (FM raiz(2 rho A) P_inst)^(2/3) / g, en una grilla paramétrica (no va al artículo todavía);
#   - cualquier freno (hélices u otro) como caja negra: debe dar m g y tomarla antes de que la nave, soltada desde el
#     reposo, suba h_max: tau = raiz(2 h_max (M_I - m) / (m g)), por nave (techo de respuesta, con fuerza infinita);
#   - con hélices fijas la energía escala como m^(5/2) / P_rej (P ~ m^1.5 y la quema dura m q / (L P_rej)).
# Salida: resultados/datos/helices_vectoriales.json (exportar_numeros.py, verificacion.py grupo Z)
# =============================================================================

import io
import sys
import parametros as p
import rutas
import minimax as mm
from exp05_condensacion_tiempo_real import q_rech
from exp02_quema_centrada import L_ef

P_INST_GRILLA = (0.5e6, 1e6, 2e6, 5e6, 10e6)    # [W] potencia de propulsión instalada, solo grilla exploratoria
NAVES_RESPUESTA = ('Flying Whales LCA60T', 'Zeppelin NT')   # techo de respuesta del freno, en el artículo
A_GRILLA = (20.0, 50.0, 100.0, 200.0)            # [m2] área total de discos, solo grilla exploratoria


def potencia(F, A, rho=None):
    """Potencia ideal [W] de hélices de área total A [m2] que dan la fuerza F [N] en vuelo estacionario."""
    rho = p.rho_aire_std if rho is None else rho
    return F ** 1.5 / (p.freno_FM * (2.0 * rho * A) ** 0.5)


def escalon_maximo(P, A, rho=None):
    """Mayor escalón [kg] que sostienen hélices de área A [m2] con toda la potencia P [W] en vertical."""
    rho = p.rho_aire_std if rho is None else rho
    return (p.freno_FM * (2.0 * rho * A) ** 0.5 * P) ** (2.0 / 3.0) / p.g


def respuesta(nom, m):
    """Techo de respuesta [s] de cualquier freno para soltar m [kg] desde el reposo: el tiempo en que la nave libre
    (después de soltar acelera M_I - m) sube h_max."""
    M_I = mm.equilibrio(nom, m)['M_I']
    a = m * p.g / (M_I - m)
    return {'m_kg': m, 'M_I_kg': M_I, 'a_mps2': a, 'tau_techo_s': (2.0 * p.h_max_precision / a) ** 0.5}


def masas_rechazo(P, cand):
    """Masa [t] de la unidad de rechazo de P [W] con las clases de la tabla 7."""
    return {'avion_t': (P / 1e3 / p.kWkg_radiador[1] / 1e3, P / 1e3 / p.kWkg_radiador[0] / 1e3),
            'acc_haces_t': P / 1e3 / cand['ACC_haces']['kWkg'][0] / 1e3,
            'acc_sistema_t': P / 1e3 / cand['ACC_sistema']['kWkg'][0] / 1e3}


def sostener(m, P_rej, q, L):
    """Quema de m [kg] al rechazo P_rej [W] con las hélices sosteniendo: duración, potencia y energía por área."""
    t = m * q / (L * P_rej)
    E_quim = m / L * p.Q_comb_H2
    P_quema = P_rej / q * p.Q_comb_H2                  # potencia química de la quema
    filas = []
    for A in p.freno_A_discos:
        P_max = potencia(m * p.g, A)
        E = P_max * t / 2.5
        filas.append({'A_m2': A, 'P_max_W': P_max, 'E_J': E, 'E_frac_quim': E / E_quim,
                      'P_frac_quema': P_max / P_quema})
    return {'P_rej_W': P_rej, 't_quema_s': t, 'E_quim_J': E_quim, 'helices': filas}


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    m = float(p.m_carga)
    q = q_rech(0.0, p.phi_objetivo)
    L = L_ef(p.phi_objetivo)
    cand = {c['clave']: c for c in rutas.leer_datos('req_masa_rechazo')['candidatos']}
    ref = sostener(m, p.P_ref_masa, q, L)
    bajo = sostener(m, p.freno_P_bajo, q, L)
    ref['rechazo'] = masas_rechazo(p.P_ref_masa, cand)
    bajo['rechazo'] = masas_rechazo(p.freno_P_bajo, cand)
    mg = p.freno_escalon_grande
    t_g = mg * q / (L * p.P_ref_masa)
    grande = {'m_kg': mg, 'P_W': [potencia(mg * p.g, A) for A in p.freno_A_discos], 't_quema_s': t_g,
              'E_J': [potencia(mg * p.g, A) * t_g / 2.5 for A in p.freno_A_discos]}
    resp = {nom: respuesta(nom, m) for nom in NAVES_RESPUESTA}
    cota = [{'P_inst_W': Pi, 'A_m2': A, 'm_max_kg': escalon_maximo(Pi, A)} for Pi in P_INST_GRILLA for A in A_GRILLA]
    rutas.guardar_datos('helices_vectoriales', {'m_kg': m, 'F_N': m * p.g, 'q_Jkg': q, 'L': L,
                                                'rho_aire': p.rho_aire_std, 'FM': p.freno_FM,
                                                'referencia': ref, 'bajo': bajo, 'escalon_grande': grande,
                                                'respuesta': resp,
                                                'cota_todas_verticales': cota})
    for caso in (ref, bajo):
        print(f"Rechazo {caso['P_rej_W']/1e6:.0f} MW: quema de 1 t en {caso['t_quema_s']/60:.1f} min; energía química "
              f"del H2 {caso['E_quim_J']/3.6e6:.0f} kWh")
        for f in caso['helices']:
            print(f"  discos {f['A_m2']:.0f} m2: P_max {f['P_max_W']/1e3:.0f} kW "
                  f"({100*f['P_frac_quema']:.1f} % de la potencia química de la quema), {f['E_J']/3.6e6:.1f} kWh "
                  f"({100*f['E_frac_quim']:.2f} % de la química)")
        r = caso['rechazo']
        print(f"  unidad de rechazo: aviones {r['avion_t'][0]:.2f}-{r['avion_t'][1]:.2f} t; haces de ACC "
              f"{r['acc_haces_t']:.1f} t; ACC completo {r['acc_sistema_t']:.0f} t")
    print(f"Escalón de {mg/1e3:.0f} t: {min(grande['P_W'])/1e6:.1f}-{max(grande['P_W'])/1e6:.1f} MW de hélices")
    print(f"  energía al rechazo de referencia: {min(grande['E_J'])/3.6e6:.0f}-{max(grande['E_J'])/3.6e6:.0f} kWh "
          f"(escala como m^2.5 / P_rej)")
    for nom, r in resp.items():
        print(f"Freno, {p.PLATAFORMAS[nom]['corto']}, {m/1e3:.0f} t: tomar m g antes de {r['tau_techo_s']:.2f} s "
              f"(M_I {r['M_I_kg']/1e3:.1f} t)")
    print('Cota "todas verticales" (exploratoria), escalón máximo [t]:')
    for Pi in P_INST_GRILLA:
        print(f"  {Pi/1e6:4.1f} MW: " + '  '.join(f"A {x['A_m2']:.0f} m2 {x['m_max_kg']/1e3:5.1f}"
                                                  for x in cota if x['P_inst_W'] == Pi))


if __name__ == '__main__':
    main()
