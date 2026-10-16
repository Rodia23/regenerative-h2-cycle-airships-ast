# -*- coding: utf-8 -*-
# [EN] Numerical verification: exact solution and convergence (Section 8).
# =============================================================================
# ver02_numerica.py — VER-02: verificación NUMÉRICA de la simulación del artículo
# (PIRT, tabla S2 del suplementario: integración numérica).
#
# A. Verificación del CÓDIGO contra una solución exacta. La misma función de la
#    simulación del artículo (drop_flujo_optimo_perdrop.correr), con el control y el
#    precondicionamiento apagados (ṁ* = 1e-12 kg/s: el recorte a +-ṁ* anula la quema)
#    y densidad constante, reduce la dinámica a
#        M_I dv/dt = m_d g - (1/2) rho Cd A v|v|,   M_I = M_ini - m_d + k M_ini,
#    cuya solución exacta desde el reposo es v = v_t tanh(t/t_v), z = v_t t_v ln cosh(t/t_v),
#    con v_t = sqrt(2 m_d g / (rho Cd A)) y t_v = M_I v_t / (m_d g). Se compara con la
#    tolerancia anterior del artículo y con una estricta.
# B. Verificación de la SOLUCIÓN (convergencia): el caso base del artículo (PF, 1 t,
#    minimax) y uno con condensador real (100 MW, tau = 5 s, ventear) con paso y
#    tolerancia cada vez más finos; la diferencia entre la configuración ADOPTADA
#    (rtol 1e-6, la de la simulación del artículo desde VER-02) y la más fina es la
#    incertidumbre numérica u_num de cada resultado; la anterior (rtol 1e-3) queda como registro.
# Salidas: resultados/datos/ver02_numerica.json, resultados/tablas/ver02_numerica.tex
# =============================================================================

import math
import numpy as np
import parametros as p
import rutas
import condensador
import drop_flujo_optimo_perdrop as sim
import inercia

# (nombre, max_step, rtol, atol). atol = None: la del artículo (1e-6, un solo número); si no, atol POR ESTADO
# = rtol x escala del estado (un atol único en J y en kg a la vez no deja ajustar la tolerancia: los estados
# de energía y de agua parten de cero y su derivada salta ~1e8 J/s al empezar el pulso)
CONFIGS = (('anterior (rtol 1e-3)', 0.05, 1e-3, None), ('adoptada (rtol 1e-6)', 0.05, 1e-6, 'escalada'),
           ('paso /2, rtol 1e-8', 0.025, 1e-8, 'escalada'), ('paso /4, rtol 1e-10', 0.0125, 1e-10, 'escalada'))
# escala de cada estado: la misma de la simulación (una sola fuente)
ESCALAS = sim.ESCALAS_ESTADO
ESCALAS_COND = sim.ESCALAS_COND


def atol_de(rt, at, con_cond):
    if at is None:
        return 1e-6
    return np.array(ESCALAS + (ESCALAS_COND if con_cond else ())) * rt
T_POST = 130.0


def exacta(t, m_d, M_I, rho, Cd, A):
    v_t = math.sqrt(2.0 * m_d * p.g / (rho * Cd * A))
    t_v = M_I * v_t / (m_d * p.g)
    return v_t * np.tanh(t / t_v), v_t * t_v * np.log(np.cosh(t / t_v)), v_t, t_v


def parte_A():
    H_orig = p.H_escala
    p.H_escala = 1e15                                   # densidad constante (solo para esta prueba)
    try:
        out = []
        for nom, ms, rt, at in (CONFIGS[0], CONFIGS[2]):
            r = sim.correr(graficar=False, verbose=False, dm_opt=1e-12, t_post=T_POST, devolver_serie=True,
                           max_step=ms, rtol=rt, atol=atol_de(rt, at, False))
            s = r['serie']
            t_d = s['t_drops'][0]
            m_d = p.m_carga
            M_I = inercia.inercial_casco(p.NAVE_BASE['V_casco'], p.NAVE_BASE['k']) - m_d   # tras soltar
            post = (s['t'] > t_d + 1.0)
            tt = s['t'][post] - t_d
            v_ex, z_ex, v_t, t_v = exacta(tt, m_d, M_I, p.rho_aire_std, p.Cd, sim.DIRIGIBLE['Area'])
            e_v = float(np.max(np.abs(s['v'][post] - v_ex) / np.maximum(np.abs(v_ex), 1e-9)))
            e_z = float(np.max(np.abs((s['h'][post] - p.h0) - z_ex) / np.maximum(np.abs(z_ex), 1e-9)))
            antes = s['t'] < t_d
            out.append({'config': nom, 'max_step': ms, 'rtol': rt, 'err_rel_v': e_v, 'err_rel_z': e_z,
                        'v_terminal_ms': v_t, 't_v_s': t_v, 'z_final_m': float(z_ex[-1]),
                        'reposo_antes_m': float(np.max(np.abs(s['h'][antes] - p.h0)))})
        return out
    finally:
        p.H_escala = H_orig


def corrida(caso, ms, rt, at):
    cond = None if caso == 'base' else condensador.config(100e6, tau_s=p.tau_condensador, modo='ventear')
    r = sim.correr(graficar=False, verbose=False, cond=cond, devolver_serie=True, max_step=ms, rtol=rt,
                   atol=atol_de(rt, at, cond is not None))
    s = r['serie']
    return {'error_max_m': r['error_max_m'], 'dh_abajo_m': r['dh_abajo_m'], 'dh_arriba_m': r['dh_arriba_m'],
            'H2_quemado_kg': float(s['m_H2'][0] - np.min(s['m_H2'])), 'agua_kg': float(s['agua'][-1])}


def parte_B():
    out = {}
    for caso in ('base', 'condensador_100MW_tau5'):
        filas = [dict(config=n, max_step=ms, rtol=rt, **corrida(caso, ms, rt, at)) for n, ms, rt, at in CONFIGS]
        ref = filas[-1]
        for f in filas:
            f['u_rel'] = {k: abs(f[k] - ref[k]) / max(abs(ref[k]), 1e-12)
                          for k in ('error_max_m', 'dh_abajo_m', 'H2_quemado_kg', 'agua_kg')}
        out[caso] = filas
    return out


def tabla(A, B):
    L = [r'\begin{tabular}{l l r r r r}', r'\toprule',
         r'Prueba & Configuración & \multicolumn{4}{c}{Error o diferencia relativa} \\',
         r'\midrule']
    for a in A:
        L.append(f"Solución exacta (control apagado) & {a['config']} & $v$: {a['err_rel_v']:.1e} & $z$: "
                 f"{a['err_rel_z']:.1e} & & \\\\")
    for caso, filas in B.items():
        L.append(r'\midrule')
        nombre = 'Caso base (minimax)' if caso == 'base' else 'Condensador 100 MW, $\\tau$ = 5 s'
        for f in filas[:-1]:
            u = f['u_rel']
            L.append(f"{nombre} & {f['config']} & exc.: {u['error_max_m']:.1e} & baja: {u['dh_abajo_m']:.1e} & "
                     f"H$_2$: {u['H2_quemado_kg']:.1e} & agua: {u['agua_kg']:.1e} \\\\")
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('ver02_numerica'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/ver02_numerica.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')


def main():
    A = parte_A()
    B = parte_B()
    rutas.guardar_datos('ver02_numerica', {'exacta': A, 'convergencia': B, 'configs': [list(c) for c in CONFIGS]})
    tabla(A, B)
    print('A. Solución exacta (control apagado, densidad constante):')
    for a in A:
        print(f"   {a['config']:22s} err v {a['err_rel_v']:.2e}  err z {a['err_rel_z']:.2e}  (v_t {a['v_terminal_ms']:.3f} m/s, "
              f"t_v {a['t_v_s']:.1f} s, z final {a['z_final_m']:.1f} m, reposo previo {a['reposo_antes_m']:.1e} m)")
    print('B. Convergencia (diferencia contra la configuración más fina):')
    for caso, filas in B.items():
        for f in filas:
            print(f"   {caso:24s} {f['config']:22s} exc {f['error_max_m']:.6f} m  H2 {f['H2_quemado_kg']:.4f} kg  "
                  f"agua {f['agua_kg']:.3f} kg  u_rel {max(f['u_rel'].values()):.1e}")


if __name__ == '__main__':
    main()
