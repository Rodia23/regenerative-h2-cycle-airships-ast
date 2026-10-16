# -*- coding: utf-8 -*-
# [EN] Requirement map over release sizes and thresholds, by simulation (check of the closed form).
# =============================================================================
# mapa_requisitos.py — Paso 1 del plan (opción B): mapa de requisitos VERIFICADO
#
# Pregunta: para soltar m_load toneladas con una tolerancia h_max, ¿qué flujo de
# H2 (y qué potencia térmica, aire, agua...) hace falta, de verdad?
#
# Para cada plataforma x carga (1 t ... su capacidad, máx. 30 t ~ contenedor ISO
# de 40 pies) x h_max {1, 3.5, 5, 7} m x D_lim {0.10 m, h/10, h/5, h/3}:
#   1) forma cerrada (minimax.py, quema limpia) -> ṁ* de partida;
#   2) simulación completa (drop_flujo_optimo_perdrop.correr) SIN reinyección
#      (el controlador solo quema, la premisa de la quema limpia) + bisección:
#      el MENOR ṁ instalado con el que la nave queda dentro de ±h_max y se
#      asienta (antes, 2026-09-27: solo la excursión hacia arriba y con el PD
#      devolviendo H2 de los tanques, lo que bajaba el flujo; hallazgo C1);
#   3) por (plataforma, carga, h_max): caso CONSERVADOR = D_lim 0.10 m y MEJOR
#      caso = el D_lim factible con menor ṁ verificado.
# Además, por (plataforma, carga, h_max) se simula D_lim = h_max con el ṁ* de la
# forma cerrada, para cuantificar cuánto se pasa del umbral: la quema limpia baja
# (1+s)·D_lim en total, así que se pasa ~2x (hallazgo Paso 1, ahora explicado).
#
# Salidas: datos/mapa_requisitos_puntos.csv (cada D_lim), datos/mapa_requisitos.csv
# (conservador y mejor caso + requisitos derivados), datos/mapa_requisitos.json,
# figuras/mapa_requisitos.png, tablas/mapa_requisitos.tex.
# Es el paso más lento del pipeline (~6 min con 12 núcleos): correr_todo --rapido lo omite.
# =============================================================================

import os
import io
import sys
import csv
import time
import math
from concurrent.futures import ProcessPoolExecutor

import parametros as p
import rutas
import minimax as mm

H_MAX_LIST       = [1.0, 3.5, 5.0, 7.0]
CARGA_MAX_MAPA_T = 30.0          # ~ masa bruta máx. de un contenedor ISO de 40 pies
T_POST           = 250.0         # s simulados después de la liberación
TOL_REL          = 0.01          # precisión de la bisección en ṁ (1 %)
H_TABLA          = p.h_max_seguridad   # h_max de la figura y la tabla del artículo
CARGAS_TABLA     = [1, 5, 10, 12.5, 20, 30]

AIRE_POR_KG_H2 = p.ratio_O2_H2 / p.w_O2_aire      # kg de aire estequiométrico por kg de H2


def cargas(nombre):
    cap = min(p.PLATAFORMAS[nombre]['carga_max'] / 1000.0, CARGA_MAX_MAPA_T)
    lista = [float(x) for x in range(1, int(cap) + 1)]
    if cap > int(cap):
        lista.append(cap)
    return lista


def candidatos_D(h):
    return sorted({round(x, 4) for x in (p.D_lim_oper, h / 10, h / 5, h / 3)})


# ── Simulación de un punto (corre en procesos separados) ────────────────────
def _simular(nombre, m_t, h, D, dm):
    import drop_flujo_optimo_perdrop as sim
    sim.usar_nave(nombre, m_t * 1000.0)
    M = mm.equilibrio(nombre, m_t * 1000.0)['M_I']             # inercia en equilibrio (inercia.py)
    t_ant = (6.0 * M * D / (mm.L_eff * p.g * dm)) ** (1.0 / 3.0)
    t_lib = max(p.intervalo_multidrop, 1.2 * t_ant + 5.0)     # espacio para la anticipación
    r = sim.correr(1, m_t, graficar=False, verbose=False, h_max=h, D_lim=D, dm_opt=dm,
                   intervalo=t_lib, t_post=T_POST, inyeccion=False)
    return {k: r[k] for k in ('dh_arriba_m', 'dh_abajo_m', 'error_max_m', 'v_max_ms', 'a_max_ms2',
                              't_estab_s', 't_anticipo', 'Energia_kWh', 'm_H2_neto')}


def verificar_punto(args):
    """Menor ṁ (±1 %) con el que la nave, solo quemando, queda dentro de ±h y se asienta."""
    nombre, m_t, h, D = args
    s = mm.resolver(nombre, m_t * 1000.0, h, D)
    dm_cf, alpha = s['dm_star'], s['alpha']      # alpha: solo informativo en la salida
    n = 0

    def sim(dm):
        nonlocal n
        n += 1
        return _simular(nombre, m_t, h, D, dm)

    cumple = lambda r: r['error_max_m'] <= h * (1 + 1e-9) and r['t_estab_s'] == r['t_estab_s']
    r0 = sim(dm_cf)
    if cumple(r0):                      # buscar hacia abajo
        hi, r_hi, lo = dm_cf, r0, None
        x = dm_cf
        while lo is None and x > 1e-3 * dm_cf:
            x /= 1.5
            rx = sim(x)
            if cumple(rx):
                hi, r_hi = x, rx
            else:
                lo = x
        if lo is None:
            lo = x
    else:                                      # buscar hacia arriba
        lo, hi, r_hi = dm_cf, None, None
        x = dm_cf
        while hi is None and x < 50 * dm_cf:
            x *= 1.5
            rx = sim(x)
            if cumple(rx):
                hi, r_hi = x, rx
            else:
                lo = x
    if hi is None:                             # ni con 50x la forma cerrada
        return {'plataforma': nombre, 'carga_t': m_t, 'h_max_m': h, 'D_lim_m': D,
                'epsilon': s['epsilon'], 'alpha': alpha, 'dm_cerrada_kgps': dm_cf,
                'dm_verificado_kgps': float('nan'), 'factible': False, 'n_sims': n}
    while hi / lo > 1.0 + TOL_REL:
        mid = math.sqrt(hi * lo)
        rm = sim(mid)
        if cumple(rm):
            hi, r_hi = mid, rm
        else:
            lo = mid
    factible = r_hi['error_max_m'] <= h * (1 + 1e-9)
    return {'plataforma': nombre, 'carga_t': m_t, 'h_max_m': h, 'D_lim_m': D,
            'epsilon': s['epsilon'], 'alpha': alpha, 'dm_cerrada_kgps': dm_cf,
            'dm_verificado_kgps': hi, 'factible': factible, 'n_sims': n, **r_hi}


def sobrepaso_D_igual_h(args):
    """D_lim = h_max con el ṁ* de la forma cerrada: cuánto se pasa del umbral."""
    nombre, m_t, h = args
    s = mm.resolver(nombre, m_t * 1000.0, h, h)
    r = _simular(nombre, m_t, h, h, s['dm_star'])
    return (nombre, m_t, h, r['error_max_m'] / h)


# ── Requisitos derivados (sin elegir equipo) ────────────────────────────────
def requisitos(m_t, dm):
    m_H2 = m_t * 1000.0 / mm.L_eff                          # H2 por entrega [kg]
    agua = (1 + p.ratio_O2_H2) * m_H2                       # agua por entrega [kg]
    return {
        'potencia_termica_MW':  dm * p.Q_comb_H2 / 1e6 if dm == dm else float('nan'),
        'aire_min_kgps':        dm * AIRE_POR_KG_H2 if dm == dm else float('nan'),
        'H2_por_entrega_kg':    m_H2,
        'energia_termica_GJ':   m_H2 * p.Q_comb_H2 / 1e9,
        'agua_por_entrega_kg':  agua,
        'Q_cond_latente_kW':    agua * p.L_v_agua / 1e3 / (p.t_cond_objetivo * 60.0),
    }


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    t0 = time.time()
    tareas = [(n, m_t, h, D) for n in p.PLATAFORMAS for m_t in cargas(n)
              for h in H_MAX_LIST for D in candidatos_D(h)]
    tareas_D = [(n, m_t, h) for n in p.PLATAFORMAS for m_t in cargas(n) for h in H_MAX_LIST]
    print(f'Puntos a verificar: {len(tareas)} (+ {len(tareas_D)} con D_lim = h_max)')
    with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 2) - 1)) as ex:
        puntos = list(ex.map(verificar_punto, tareas, chunksize=4))
        sobrep = list(ex.map(sobrepaso_D_igual_h, tareas_D, chunksize=4))
    print(f'Simulaciones: {sum(r["n_sims"] for r in puntos) + len(tareas_D)}  ({time.time() - t0:.0f} s)')

    # ── CSV con todos los D_lim ──
    campos = ['plataforma', 'carga_t', 'h_max_m', 'D_lim_m', 'epsilon', 'alpha', 'dm_cerrada_kgps',
              'dm_verificado_kgps', 'factible', 'n_sims', 'dh_arriba_m', 'dh_abajo_m', 'error_max_m',
              'v_max_ms', 'a_max_ms2', 't_estab_s', 't_anticipo', 'Energia_kWh', 'm_H2_neto']
    with open(rutas.dato('mapa_requisitos_puntos', 'csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for r in puntos:
            w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items() if k in campos})

    # ── Resumen por (plataforma, carga, h_max): conservador y mejor caso ──
    sob = {(n, m_t, h): x for n, m_t, h, x in sobrep}
    filas = []
    for n in p.PLATAFORMAS:
        for m_t in cargas(n):
            for h in H_MAX_LIST:
                grupo = [r for r in puntos if r['plataforma'] == n and r['carga_t'] == m_t and r['h_max_m'] == h]
                cons = next(r for r in grupo if abs(r['D_lim_m'] - p.D_lim_oper) < 1e-9)
                fact = [r for r in grupo if r['factible']]
                mejor = min(fact, key=lambda r: r['dm_verificado_kgps']) if fact else None
                fila = {'plataforma': n, 'corto': p.PLATAFORMAS[n]['corto'], 'carga_t': m_t, 'h_max_m': h,
                        'epsilon': cons['epsilon'],
                        'cons_dm_cerrada_kgps': cons['dm_cerrada_kgps'],
                        'cons_dm_verificado_kgps': cons['dm_verificado_kgps'], 'cons_factible': cons['factible'],
                        'cons_error_max_m': cons.get('error_max_m', float('nan')),
                        'cons_t_estab_s': cons.get('t_estab_s', float('nan')),
                        'mejor_D_lim_m': mejor['D_lim_m'] if mejor else float('nan'),
                        'mejor_dm_verificado_kgps': mejor['dm_verificado_kgps'] if mejor else float('nan'),
                        'sobrepaso_D_igual_h': sob[(n, m_t, h)]}
                ref = mejor if mejor else cons
                for k in ('dh_arriba_m', 'dh_abajo_m', 'error_max_m', 'v_max_ms', 'a_max_ms2', 't_estab_s',
                          't_anticipo'):
                    fila['mejor_' + k] = ref.get(k, float('nan')) if mejor else float('nan')
                for k, v in requisitos(m_t, fila['mejor_dm_verificado_kgps']).items():
                    fila['mejor_' + k] = v
                fila['cons_potencia_termica_MW'] = requisitos(m_t, fila['cons_dm_verificado_kgps'])['potencia_termica_MW']
                filas.append(fila)
    with open(rutas.dato('mapa_requisitos', 'csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        for r in filas:
            w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})

    # ── JSON para numeros.tex ──
    en_capacidad = {}
    for n in p.PLATAFORMAS:
        cap = cargas(n)[-1]
        r = next(f for f in filas if f['plataforma'] == n and f['carga_t'] == cap and f['h_max_m'] == H_TABLA)
        en_capacidad[n] = {k: r[k] for k in ('corto', 'carga_t', 'epsilon', 'cons_dm_verificado_kgps',
                                             'mejor_dm_verificado_kgps', 'mejor_D_lim_m',
                                             'mejor_potencia_termica_MW', 'mejor_aire_min_kgps',
                                             'mejor_t_estab_s', 'mejor_error_max_m')}
    vals_sob = [f['sobrepaso_D_igual_h'] for f in filas]
    cons_malos = [f for f in filas if f['cons_dm_verificado_kgps'] > f['cons_dm_cerrada_kgps'] * (1 + TOL_REL)]
    rutas.guardar_datos('mapa_requisitos', {
        'h_max_tabla_m': H_TABLA, 'aire_por_kg_H2': AIRE_POR_KG_H2,
        'puntos': len(puntos), 'filas': len(filas),
        'sobrepaso_D_igual_h_min': min(vals_sob), 'sobrepaso_D_igual_h_max': max(vals_sob),
        'conservador_forma_cerrada_insuficiente': len(cons_malos),
        'conservador_factor_max': max(f['cons_dm_verificado_kgps'] / f['cons_dm_cerrada_kgps'] for f in filas),
        'mejor_sin_solucion': sum(1 for f in filas if f['mejor_dm_verificado_kgps'] != f['mejor_dm_verificado_kgps']),
        'en_capacidad_h5': en_capacidad,
    })

    # ── Consola ──
    print(f'\nD_lim = h_max con la forma cerrada: se pasa del umbral {min(vals_sob):.2f}x a {max(vals_sob):.2f}x h_max')
    print(f'Conservador (D_lim = {p.D_lim_oper} m): la forma cerrada se queda corta en {len(cons_malos)} de {len(filas)} casos')
    print(f'\n{"Plataforma":<22}{"carga":>6}{"h":>5}{"eps":>7}{"ṁ cons":>9}{"ṁ mejor":>9}{"D_lim":>7}'
          f'{"P_th MW":>9}{"aire kg/s":>10}{"|Δh| m":>8}{"t_est s":>8}')
    for f in filas:
        if f['h_max_m'] == H_TABLA and f['carga_t'] in CARGAS_TABLA:
            print(f'{f["plataforma"]:<22}{f["carga_t"]:>6g}{f["h_max_m"]:>5g}{f["epsilon"]:>7.3f}'
                  f'{f["cons_dm_verificado_kgps"]:>9.1f}{f["mejor_dm_verificado_kgps"]:>9.1f}{f["mejor_D_lim_m"]:>7.2f}'
                  f'{f["mejor_potencia_termica_MW"]:>9.0f}{f["mejor_aire_min_kgps"]:>10.0f}'
                  f'{f["mejor_error_max_m"]:>8.2f}{f["mejor_t_estab_s"]:>8.1f}')

    graficar(filas)
    tabla(filas)
    print(f'\nTotal: {time.time() - t0:.0f} s')


def graficar(filas):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import plot_style  # noqa
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for n, d in p.PLATAFORMAS.items():
        fs = sorted((f for f in filas if f['plataforma'] == n and f['h_max_m'] == H_TABLA), key=lambda f: f['carga_t'])
        x = [f['carga_t'] for f in fs]
        mk = 'o' if len(x) == 1 else None
        ax.plot(x, [f['cons_dm_verificado_kgps'] for f in fs], color=d['color'], lw=2.2, ls='-', marker=mk, ms=7)
        ax.plot(x, [f['mejor_dm_verificado_kgps'] for f in fs], color=d['color'], lw=2.0, ls='--',
                marker=('s' if mk else None), ms=6, mfc='white')
    ax.set_yscale('log')
    ax.axvline(CARGA_MAX_MAPA_T, color='0.45', ls=':', lw=1.2)
    ax.text(CARGA_MAX_MAPA_T - 0.4, 0.35, 'ISO 40-ft container (~30 t)', rotation=90, ha='right',
            va='center', fontsize=8, color='0.35', transform=ax.get_xaxis_transform())
    ax.set_xlim(0, 31)
    ax.set_xlabel(r'Payload released  $m_\mathrm{load}$  [t]')
    ax.set_ylabel(r'Verified installed H$_2$ flow  $\dot{m}_{H_2}$  [kg/s]')
    ax.grid(True, which='both', alpha=0.2)
    leg1 = ax.legend([Line2D([0], [0], color=d['color'], lw=2.2) for d in p.PLATAFORMAS.values()],
                     list(p.PLATAFORMAS), loc='upper left', title='Airship (up to its payload capacity)',
                     fontsize=8.5, title_fontsize=8.5)
    ax.add_artist(leg1)
    ax.legend([Line2D([0], [0], color='k', lw=2.2, ls='-'), Line2D([0], [0], color='k', lw=2.0, ls='--')],
              [rf'Conservative ($D_\mathrm{{lim}}$ = {p.D_lim_oper:.2f} m)', r'Best verified $D_\mathrm{lim}$'],
              loc='lower right', fontsize=8.5, title=rf'$h_\mathrm{{max}}$ = {H_TABLA:g} m', title_fontsize=8.5)
    fig.tight_layout()
    fig.savefig(rutas.figura('mapa_requisitos'), dpi=200, bbox_inches='tight')
    plt.close(fig)


def tabla(filas):
    def num(x, dec):
        if x != x:
            return '--'
        s = f'{x:,.{dec}f}' if abs(x) >= 10000 else f'{x:.{dec}f}'
        return s.replace(',', '{,}')
    L = []
    L.append(r'\begin{table}[H]')
    L.append(r'\caption{Verified requirement map at the safety threshold $h_\text{max} = ' + f'{H_TABLA:g}' +
             r'$~m (single release). $\dot{m}$: minimum installed H\textsubscript{2} flow for which the full '
             r'nonlinear simulation keeps $|\Delta h| \le h_\text{max}$, conservative ($D_\text{lim} = ' +
             f'{p.D_lim_oper:.2f}' + r'$~m) and with the best verified preconditioning depth. $P_\text{th}$, air '
             r'intake (stoichiometric), $|\Delta h|_\text{max}$ and $t_\text{stab}$ refer to the best case. '
             r'Generated by \texttt{modelo/mapa\_requisitos.py}.}')
    L.append(r'\label{tab:mapa_requisitos}')
    L.append(r'\centering\small')
    L.append(r'\begin{tabular}{@{} l r r r r r r r r r @{}}')
    L.append(r'\toprule')
    L.append(r'\textbf{Platform} & \makecell{$m_\text{load}$\\(t)} & $\varepsilon$ & '
             r'\makecell{$\dot{m}$ cons.\\(kg/s)} & \makecell{$\dot{m}$ best\\(kg/s)} & '
             r'\makecell{$D_\text{lim}$\\(m)} & \makecell{$P_\text{th}$\\(MW)} & \makecell{Air\\(kg/s)} & '
             r'\makecell{$|\Delta h|_\text{max}$\\(m)} & \makecell{$t_\text{stab}$\\(s)} \\')
    L.append(r'\midrule')
    nombres = {'LZ 129 Hindenburg': r'LZ~129 Hindenburg', 'Zeppelin NT': r'Zeppelin~NT'}
    for i, n in enumerate(p.PLATAFORMAS):
        fs = [f for f in filas if f['plataforma'] == n and f['h_max_m'] == H_TABLA and f['carga_t'] in CARGAS_TABLA]
        for j, f in enumerate(sorted(fs, key=lambda f: f['carga_t'])):
            L.append(f'{nombres.get(n, n) if j == 0 else ""} & {f["carga_t"]:g} & {num(f["epsilon"], 3)} & '
                     f'{num(f["cons_dm_verificado_kgps"], 1)} & {num(f["mejor_dm_verificado_kgps"], 1)} & '
                     f'{num(f["mejor_D_lim_m"], 2)} & {num(f["mejor_potencia_termica_MW"], 0)} & '
                     f'{num(f["mejor_aire_min_kgps"], 0)} & {num(f["mejor_error_max_m"], 2)} & '
                     f'{num(f["mejor_t_estab_s"], 0)} \\\\')
        if i < len(p.PLATAFORMAS) - 1:
            L.append(r'\midrule')
    L.append(r'\bottomrule')
    L.append(r'\end{tabular}')
    L.append(r'\end{table}')
    with open(rutas.tabla('mapa_requisitos'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')


def solo_figuras():
    """Redibuja figura y tabla desde datos/mapa_requisitos.csv, sin volver a simular."""
    with open(rutas.dato('mapa_requisitos', 'csv'), encoding='utf-8') as f:
        filas = []
        for r in csv.DictReader(f):
            for k, v in r.items():
                if k not in ('plataforma', 'corto'):
                    r[k] = (v == 'True') if v in ('True', 'False') else float(v)
            filas.append(r)
    graficar(filas)
    tabla(filas)
    print('Figura y tabla regeneradas desde el csv (sin simular).')


if __name__ == '__main__':
    if '--solo-figuras' in sys.argv:
        solo_figuras()
    else:
        main()
