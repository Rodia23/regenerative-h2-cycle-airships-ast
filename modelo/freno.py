# -*- coding: utf-8 -*-
# [EN] Hold-down device: sharing of the release, recovery time and energy (Section 5.2, Table S10).
# =============================================================================
# freno.py — paper 1, sección 5.2 ("With a hold-down device"). El dispositivo de sujeción (freno) es la 4.ª caja
# negra del ciclo: cualquier cosa que sostenga la nave hacia abajo mientras se quema (las hélices vectoriales que la
# nave ya tiene, alimentadas por el motor del ciclo o, más allá de su alcance, por la planta de la nave, o un cable a
# tierra si la nave está amarrada).
# Con un rechazo instalado P (parametros.freno_P_grilla):
#   - el ciclo compensa EN TIEMPO REAL solo el escalón m_e(P) cuya quema limpia pide P. La ecuación 7 con la inercia
#     de la liberación completa (s de m) da ṁ*(m_e) = ṁ*(m) (m_e/m)^1.5, así que m_e = m (P / P_rt)^(2/3), con
#     P_rt = ṁ*(m) q(0, phi) el rechazo en tiempo real de todo el escalón (sección 5.3);
#   - el freno toma el resto, F = (m - m_e) g, que baja linealmente a cero mientras la quema quita la sustentación;
#   - toda la quema pasa por la misma unidad: dura t = (m / L(phi)) q / P, el TIEMPO DE RECUPERACIÓN (mismo L y q
#     que la sección 5.3: 10.6 min por tonelada con 10 MW);
#   - hélices ideales (freno_A_discos, freno_FM): P_hél = F^1.5 / (FM raiz(2 rho A)), E = P_max t / 2.5;
#   - masa de la unidad de rechazo a P con las clases de la tabla 7 (haces de ACC: condensan a la temperatura del
#     ciclo; intercambiadores de avión: 80-160 °C), y la P con la que los haces igualan la carga útil del LCA60T;
#   - respuesta: el techo (helices_vectoriales.respuesta) y, con fuerza finita (1 + mu) m g (mu del margen del autor)
#     y el mejor control (a tope desde la liberación), tau <= raiz(h_max / (a G(mu))),
#     G = u (1 - mu u / 2), (1 + mu)(1 - e^-u) = mu u;
#   - la liviandad certificada del Zeppelin NT (FAA TCDS AS1CE p. 2) frente a su carga útil;
#   - la energía de regenerar el H2 del escalón (EXP-16, 1 atm), para comparar con la del freno;
#   - chequeo con la simulación completa (drop_flujo_optimo_perdrop.correr(helices=...)): 1 t con 10 MW en el LCA60T
#     y en el NT, con quema previa y el freno de fuerza (1 + mu) m g.
# Salidas: resultados/datos/freno.json, resultados/tablas/art_freno.tex, resultados/figuras/freno_recuperacion.png
# (verificacion.py, grupo Z; exportar_numeros.py, macros Freno*).
# =============================================================================

import io
import math
import sys
import numpy as np
import parametros as p
import rutas
import minimax as mm
from exp05_condensacion_tiempo_real import q_rech
from exp02_quema_centrada import L_ef
from helices_vectoriales import potencia, masas_rechazo, respuesta
import drop_flujo_optimo_perdrop as sim   # aquí arriba: al importarse envuelve sys.stdout

Q = q_rech(0.0, p.phi_objetivo)                  # calor por kg de H2 (sección 5.3)
L_CALOR = L_ef(p.phi_objetivo)                   # H2 por tonelada para el calor y el tiempo (sección 5.3)
L_TASA = mm.L_eff                                # la tasa de la ecuación 7 se dimensiona con phi = 1
MU = 1.0 / (1.0 - p.hel_margen_autor[0]) - 1.0   # fuerza máxima del freno (1 + mu) m g
KWH = 3.6e6


def u_retardo(mu):
    """Raíz positiva de (1 + mu)(1 - e^-u) = mu u (instante, en unidades de tau, en que la subida se detiene)."""
    lo, hi = 1e-12, (1.0 + mu) / mu
    for _ in range(200):
        u = 0.5 * (lo + hi)
        if (1.0 + mu) * (1.0 - math.exp(-u)) - mu * u > 0.0:
            lo = u
        else:
            hi = u
    return 0.5 * (lo + hi)


def g_retardo(mu):
    """Subida mínima, en unidades de a tau^2, con retardo de primer orden y fuerza (1 + mu) m g mandada a tope."""
    u = u_retardo(mu)
    return u * (1.0 - mu * u / 2.0)


def tiempo_real(nom, m):
    """Rechazo en tiempo real de todo el escalón m [kg] (sección 5.3) y la inercia con la carga."""
    M_I = mm.equilibrio(nom, m)['M_I']
    lim = mm.limpio(M_I, m, p.D_lim_oper)
    return {'M_I_kg': M_I, 's': lim['s'], 'dm_star': lim['dm_star'], 'P_rt_W': lim['dm_star'] * Q,
            't_nat_s': (m / L_TASA) / lim['dm_star']}


def reparto(nom, m, P, cand):
    """Con el rechazo instalado P [W]: lo que el ciclo toma en tiempo real, lo que toma el freno, cuánto dura la
    recuperación, la potencia y la energía del freno y la masa de la unidad de rechazo."""
    rt = tiempo_real(nom, m)
    m_e = m * min(1.0, P / rt['P_rt_W']) ** (2.0 / 3.0)
    F = (m - m_e) * p.g
    t = (m / L_CALOR) * Q / P
    P_hel = [potencia(F, A) for A in p.freno_A_discos]
    mr = masas_rechazo(P, cand)
    return {'nave': nom, 'corto': p.PLATAFORMAS[nom]['corto'], 'm_kg': m, 'P_W': P, 'P_rt_W': rt['P_rt_W'],
            's': rt['s'], 'M_I_kg': rt['M_I_kg'], 'dm_star_kgps': rt['dm_star'], 't_nat_s': rt['t_nat_s'],
            'm_e_kg': m_e, 'frac_freno': (m - m_e) / m, 'F_N': F, 't_rec_s': t,
            'P_hel_W': P_hel, 'E_hel_J': [x * t / 2.5 for x in P_hel],
            'masa_acc_haces_t': mr['acc_haces_t'], 'masa_avion_t': mr['avion_t']}


def simular(nom, m, P):
    """Chequeo con la simulación completa: quema previa y quema a P / q, con el freno de fuerza (1 + mu) m g, sus
    discos menores (la energía con los mayores escala como raiz(A1/A2))."""
    sim.usar_nave(nom, m)
    dm = P / Q
    t_quema = m / (L_TASA * dm)
    hel = {'F_max': (1.0 + MU) * m * p.g, 'tau': p.hel_tau, 'wn': p.hel_wn, 'zeta': p.hel_zeta,
           'A': p.freno_A_discos[0], 'FM': p.freno_FM, 'bidireccional': True}
    r = sim.correr(graficar=False, verbose=False, toneladas_por_drop=m / 1000.0, dm_opt=dm, helices=hel,
                   t_post=1.3 * t_quema + 120.0, max_step=0.5)
    sim.usar_nave(p.nombre_dirigible)
    return {'nave': nom, 'm_kg': m, 'P_W': P, 'dm_kgps': dm, 'sube_m': r['dh_arriba_m'], 'baja_m': r['dh_abajo_m'],
            'H2_kg': r['m_H2_bruto'], 'E_hel_J': r['E_hel_J'], 'P_hel_max_W': r['P_hel_max_W'],
            'satura': r['satura'], 'A_m2': p.freno_A_discos[0]}


def anticipacion(nom, m, P):
    """El freno ANTICIPA como la quema limpia (revisión 2026-10-06, ronda 2). La deducción de la ecuación 7 vale con la
    tasa de fuerza F' en lugar de L g ṁ: F'* = L g ṁ*(m), la fuerza llega al escalón en t_F = m g / F'* = t_nat, la
    fracción 1/(1+s) antes de soltar, y la nave baja (1+s) D_lim sin subir. La quema a P (sin quema previa) solo
    agrega fuerza hacia abajo después de soltar (P/q es un 2 % de ṁ* con 10 MW); el lazo la absorbe."""
    rt = tiempo_real(nom, m)
    Fp = L_TASA * p.g * rt['dm_star']
    t_F = m * p.g / Fp
    t_previa = t_F / (1.0 + rt['s'])
    return {'nave': nom, 'm_kg': m, 'P_W': P, 'rampa_Nps': Fp, 't_nat_s': rt['t_nat_s'],
            't_F_s': t_F, 't_previa_s': t_previa, 's': rt['s'], 'M_I_kg': rt['M_I_kg'],
            'bajada_m': (1.0 + rt['s']) * p.D_lim_oper, 'v0_mps': 3.0 * p.D_lim_oper / t_previa}


def _corre_freno(nom, m, P, tau, anticipa):
    """Simulación con el freno de fuerza (1 + mu) m g y retardo tau, SIN quema previa; con anticipa, la orden sube en
    rampa (adelantada un tau). Devuelve la subida y la bajada [m]."""
    a = anticipacion(nom, m, P)
    hel = {'F_max': (1.0 + MU) * m * p.g, 'tau': tau, 'wn': p.hel_wn, 'zeta': p.hel_zeta,
           'A': p.freno_A_discos[0], 'FM': p.freno_FM, 'bidireccional': True}
    if anticipa:
        hel['anticipar'] = {'rampa': a['rampa_Nps'], 't_previa': a['t_previa_s'], 'adelanto': tau}
    dm = P / Q
    sim.usar_nave(nom, m)
    r = sim.correr(graficar=False, verbose=False, toneladas_por_drop=m / 1000.0, dm_opt=dm, helices=hel,
                   precond=False, t_post=1.3 * m / (L_TASA * dm) + 120.0, max_step=0.5)
    sim.usar_nave(p.nombre_dirigible)
    return r['dh_arriba_m'], -r['dh_abajo_m']


def tau_max(nom, m, P, anticipa, lo=0.2, hi=6.0, n=9):
    """Retardo máximo del freno con |Δh| <= h_max hacia arriba y hacia abajo (bisección, sin quema previa)."""
    def cumple(tau):
        sube, baja = _corre_freno(nom, m, P, tau, anticipa)
        return max(sube, baja) <= p.h_max_precision
    if not cumple(lo):
        return lo
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if cumple(mid) else (lo, mid)
    return lo


def _fmt_miles(x, dec=0):
    s = f'{x:,.{dec}f}'
    return s.replace(',', '{,}')


def _rango(a, b, dec=0):
    ta, tb = _fmt_miles(a, dec), _fmt_miles(b, dec)
    return ta if ta == tb else f'{ta}--{tb}'


def _tiempo(t_s):
    return f'{t_s / 3600.0:.1f}~h' if t_s >= 7200.0 else f'{t_s / 60.0:.1f}~min'


def tabla(filas):
    """resultados/tablas/art_freno.tex (inglés)."""
    nombres = {'FW': 'LCA60T', 'ZNT': 'Zeppelin~NT'}
    L = ['% GENERADO por modelo/freno.py. NO EDITAR A MANO.',
         r'\begin{tabular}{@{} l r r r r r r r @{}}', r'\toprule',
         # (autor, 2026-10-08) encabezado en tres filas (nombre, resto del nombre, unidad) para angostar las columnas
         # (autor, 2026-10-09) ANTES: Real-time step; m_e es el escalón que absorbe la quema sola, no un flujo
         r'\textbf{Airship,} & \textbf{$P_\text{cap}$} & \textbf{Absorbed} & \textbf{Device} & '
         r'\textbf{Recovery} & \textbf{Device} & \textbf{Device} & \textbf{Rejection} \\',
         r'\textbf{step} & \textbf{(MW)} & \textbf{step $m_e$} & \textbf{force} & \textbf{time} & '
         r'\textbf{power} & \textbf{energy} & \textbf{unit} \\',
         r' & & \textbf{(kg)} & \textbf{(kN)} & & \textbf{(MW)} & \textbf{(kWh)} & \textbf{(t)} \\', r'\midrule']
    previo = None
    for f in filas:
        grupo = (f['corto'], f['m_kg'])
        if previo is not None and grupo != previo:
            L.append(r'\addlinespace')
        etiqueta = f"{nombres[f['corto']]}, {f['m_kg'] / 1000:.0f}~t" if grupo != previo else ''
        previo = grupo
        dec_p = 2 if max(f['P_hel_W']) < 1e6 else 1
        L.append(f"{etiqueta} & {f['P_W'] / 1e6:.0f} & {_fmt_miles(f['m_e_kg'])} & {f['F_N'] / 1e3:.1f} & "
                 f"{_tiempo(f['t_rec_s'])} & {_rango(min(f['P_hel_W']) / 1e6, max(f['P_hel_W']) / 1e6, dec_p)} & "
                 f"{_rango(min(f['E_hel_J']) / KWH, max(f['E_hel_J']) / KWH)} & {f['masa_acc_haces_t']:.1f}"
                 f"{'' if f['cabe'] else chr(92) + 'textsuperscript{a}'} \\\\")
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('art_freno'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')


def tabla_rev(filas):
    """resultados/tablas/art_freno_rev.tex: propuesta de la sección 5 fusionada (autor, 2026-10-09). Solo el LCA60T
    (el NT va en el texto), el tiempo de recuperación antes del freno, fuerza y potencia juntas bajo "Hold-down
    device" y sin la energía de las hélices. (autor, 2026-10-09) Sin la columna de la unidad de rechazo: es
    P_cap / ReqACCHacesKWkg, va en una frase del texto."""
    for f in filas:
        if f['corto'] == 'FW' and not f['cabe']:
            sys.exit('tabla_rev: una unidad de rechazo del LCA60T no cabe en la carga; la tabla ya no la marca')
    L = ['% GENERADO por modelo/freno.py. NO EDITAR A MANO.',
         r'\begin{tabular}{@{} l r r r r r @{}}', r'\toprule',
         # (autor, 2026-10-09) "Hold-down device" solo en la fila de arriba; los demás encabezados empiezan debajo
         r' & & & & \multicolumn{2}{c}{\textbf{Hold-down device}} \\', r'\cmidrule(l){5-6}',
         r'\textbf{Release} & \textbf{$P_\text{cap}$} & \textbf{Absorbed} & \textbf{Recovery} & \textbf{Force} & '
         r'\textbf{Power} \\',
         r' & \textbf{(MW)} & \textbf{by the burn,} & \textbf{time} & \textbf{(kN)} & \textbf{(MW)} \\',
         r' & & \textbf{$m_e$ (kg)} & & & \\', r'\midrule']
    previo = None
    for f in (f for f in filas if f['corto'] == 'FW'):
        if previo is not None and f['m_kg'] != previo:
            L.append(r'\addlinespace')
        etiqueta = f"{f['m_kg'] / 1000:.0f}~t" if f['m_kg'] != previo else ''
        previo = f['m_kg']
        dec_p = 2 if max(f['P_hel_W']) < 1e6 else 1
        L.append(f"{etiqueta} & {f['P_W'] / 1e6:.0f} & {_fmt_miles(f['m_e_kg'])} & {_tiempo(f['t_rec_s'])} & "
                 f"{f['F_N'] / 1e3:.1f} & {_rango(min(f['P_hel_W']) / 1e6, max(f['P_hel_W']) / 1e6, dec_p)} \\\\")
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('art_freno_rev'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')


def figura(cand):
    """resultados/figuras/freno_recuperacion.png: rechazo frente al tiempo de recuperación, LCA60T."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401
    nom = 'Flying Whales LCA60T'
    colores = ('#2a78d6', '#eb6834', '#1baf7a')       # paleta categórica de referencia, validada (3 series)
    estilos = ('-', '--', '-.')
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    kwkg = cand['ACC_haces']['kWkg'][0]
    P_lim = p.PLATAFORMAS[nom]['carga_max'] / 1000.0 * kwkg       # MW con los que los haces pesan la carga útil
    curvas = []
    for (m, c, ls) in zip(p.freno_escalones[nom], colores, estilos):
        rt = tiempo_real(nom, m)
        Qt = (m / L_CALOR) * Q
        t = np.logspace(math.log10(rt['t_nat_s']), math.log10(p.freno_t_fig_h * 3600.0), 200)
        ax.loglog(t / 60.0, Qt / t / 1e6, color=c, ls=ls, label=f'{m / 1000:.0f} t')
        ax.plot(rt['t_nat_s'] / 60.0, Qt / rt['t_nat_s'] / 1e6, 'o', ms=8, color=c, mec='white', mew=2, zorder=5)
        ax.annotate(f'{m / 1000:.0f} t', xy=(t[-1] / 60.0, Qt / t[-1] / 1e6), xytext=(4, 0),
                    textcoords='offset points', va='center', fontsize=10, color='#333333')
        curvas.append((rt['t_nat_s'] / 60.0, Qt / rt['t_nat_s'] / 1e6))
    ax.axhline(P_lim, color='#52514e', lw=1.0, ls=':')
    ax.annotate('condenser bundles = LCA60T payload', xy=(0.035, P_lim), xytext=(0, -11),
                textcoords='offset points', fontsize=9, color='#333333')
    ax.axhline(p.P_ref_masa / 1e6, color='#52514e', lw=0.8, ls='-', alpha=0.5)
    ax.annotate(f'{p.P_ref_masa / 1e6:.0f} MW (reference unit)', xy=(0.035, p.P_ref_masa / 1e6), xytext=(0, -11),
                textcoords='offset points', fontsize=9, color='#333333')
    from matplotlib.lines import Line2D
    marca = Line2D([], [], ls='none', marker='o', ms=8, color='#52514e', mec='white', mew=2,
                   label='buoyancy control alone\n(real time, no device)')
    ax.set_xlabel('Recovery time (min)')
    ax.set_ylabel(r'Installed heat rejection $P_\mathrm{cap}$ (MW)')
    ax.set_xlim(0.03, p.freno_t_fig_h * 60.0 * 2.5)
    h, _ = ax.get_legend_handles_labels()
    ax.legend(handles=h + [marca], title='Step released', loc='upper right', frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(rutas.figura('freno_recuperacion'), dpi=180, bbox_inches='tight')
    plt.close(fig)
    return {'P_lim_haces_W': P_lim * 1e6, 'kWkg_haces': kwkg}


def main():
    if (sys.stdout.encoding or '').lower() != 'utf-8':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    cand = {c['clave']: c for c in rutas.leer_datos('req_masa_rechazo')['candidatos']}
    # Auditoría 2026-10-06 (M1, M2): la grilla termina en la potencia con la que los haces pesan la carga útil del
    # LCA60T, y cada fila dice si la unidad de rechazo cabe en la carga útil de SU nave
    kwkg = cand['ACC_haces']['kWkg'][0]
    limite = {nom: {'P_lim_W': p.PLATAFORMAS[nom]['carga_max'] / 1000.0 * kwkg * 1e6} for nom in p.freno_escalones}
    for nom, lc in limite.items():
        lc['t_rec_1t_s'] = (1000.0 / L_CALOR) * Q / lc['P_lim_W']
    grilla = list(p.freno_P_grilla) + [limite['Flying Whales LCA60T']['P_lim_W']]
    filas = [reparto(nom, m, P, cand) for nom, escalones in p.freno_escalones.items() for m in escalones
             for P in grilla]
    for f in filas:
        f['cabe'] = f['masa_acc_haces_t'] * 1e3 <= p.PLATAFORMAS[f['nave']]['carga_max'] * (1.0 + 1e-9)
    resp = {}
    a_G = g_retardo(MU)
    for nom in p.freno_escalones:
        r = respuesta(nom, 1000.0)
        resp[nom] = dict(r, mu=MU, G=a_G, u=u_retardo(MU),
                         tau_primer_orden_s=r['tau_techo_s'] / math.sqrt(2.0 * a_G))
    nt = p.HELICES_NAVES['Zeppelin NT']
    b16 = {f['caso']: f for f in rutas.leer_datos('exp16_balance_neto')['filas']}['B']
    regen_kwh_t = [b16['E_electrolisis'] + x for x in b16['E_compresion']]     # por tonelada liberada, 1 atm
    sims = [simular(nom, 1000.0, p.P_ref_masa) for nom in ('Flying Whales LCA60T', 'Zeppelin NT')]
    for s in sims:
        f = next(x for x in filas if x['nave'] == s['nave'] and x['m_kg'] == s['m_kg'] and x['P_W'] == s['P_W'])
        s['E_cerrada_J'] = f['E_hel_J'][0]
        s['razon_E'] = s['E_hel_J'] / f['E_hel_J'][0]
        s['H2_esperado_kg'] = s['m_kg'] / L_TASA
    # Ronda 2 (2026-10-06): anticipar con el freno y prescindir de la quema previa
    ant = {nom: anticipacion(nom, 1000.0, p.P_ref_masa) for nom in p.freno_escalones}
    for f in filas:                       # sin quema previa el freno toma el escalón entero (m_e = 0)
        f['factor_P_sin_previa'] = (f['m_kg'] / (f['m_kg'] - f['m_e_kg'])) ** 1.5
    nt_nom = 'Zeppelin NT'
    taus = {'espera_s': tau_max(nt_nom, 1000.0, p.P_ref_masa, False),
            'anticipa_s': tau_max(nt_nom, 1000.0, p.P_ref_masa, True)}
    sube_e, baja_e = _corre_freno(nt_nom, 1000.0, p.P_ref_masa, taus['anticipa_s'], False)
    sube_a, baja_a = _corre_freno(nt_nom, 1000.0, p.P_ref_masa, taus['anticipa_s'], True)
    taus.update({'nave': nt_nom, 'espera_con_tau_ant': [sube_e, baja_e], 'anticipa_con_tau_ant': [sube_a, baja_a]})
    sube_0, baja_0 = _corre_freno(nt_nom, 1000.0, p.P_ref_masa, p.hel_tau, False)
    taus['sin_previa_tau_base'] = [sube_0, baja_0]
    fig = figura(cand)
    tabla(filas)
    tabla_rev(filas)
    rutas.guardar_datos('freno', {'Q_Jkg': Q, 'L_calor': L_CALOR, 'L_tasa': L_TASA, 'mu': MU,
                                  'margen': p.hel_margen_autor[0], 'filas': filas, 'respuesta': resp,
                                  'NT': {'liviandad_cert_kg': nt['liviandad_cert_kg'],
                                         'carga_kg': p.PLATAFORMAS['Zeppelin NT']['carga_max']},
                                  'regen_kWh_por_t': regen_kwh_t, 'simulacion': sims, 'figura': fig,
                                  'anticipacion': ant, 'tau_nt': taus, 'limite_carga': limite})
    for f in filas:
        print(f"{f['corto']:4s} {f['m_kg'] / 1e3:3.0f} t, {f['P_W'] / 1e6:3.0f} MW: m_e {f['m_e_kg']:6.1f} kg "
              f"({100 * (1 - f['frac_freno']):5.2f} %), freno {f['F_N'] / 1e3:6.1f} kN, recuperación "
              f"{f['t_rec_s'] / 60:7.1f} min, hélices {min(f['P_hel_W']) / 1e6:6.2f}-{max(f['P_hel_W']) / 1e6:6.2f} MW "
              f"y {min(f['E_hel_J']) / KWH:8.0f}-{max(f['E_hel_J']) / KWH:8.0f} kWh; haces {f['masa_acc_haces_t']:.1f} t")
    for nom, r in resp.items():
        print(f"{p.PLATAFORMAS[nom]['corto']}: techo {r['tau_techo_s']:.2f} s; con fuerza (1+{MU:.2f}) m g y el mejor "
              f"control {r['tau_primer_orden_s']:.2f} s (G = {r['G']:.3f})")
    print(f"Haces de ACC iguales a la carga útil del LCA60T con {fig['P_lim_haces_W'] / 1e6:.0f} MW")
    for nom, a in ant.items():
        print(f"Anticipar {p.PLATAFORMAS[nom]['corto']}: rampa {a['rampa_Nps'] / 1e3:.2f} kN/s, escalón en "
              f"{a['t_F_s']:.1f} s ({a['t_previa_s']:.1f} s antes de soltar), bajada {a['bajada_m']:.3f} m")
    print(f"NT, retardo máximo con |dh| <= {p.h_max_precision} m: esperando {taus['espera_s']:.2f} s, anticipando "
          f"{taus['anticipa_s']:.2f} s; con este último, esperando sube/baja {sube_e:.2f}/{baja_e:.2f} m")
    for s in sims:
        print(f"Simulación {p.PLATAFORMAS[s['nave']]['corto']}, 1 t a {s['P_W'] / 1e6:.0f} MW: sube {s['sube_m']:.3f} m, "
              f"H2 {s['H2_kg']:.1f} kg (esperado {s['H2_esperado_kg']:.1f}), energía sim / cerrada "
              f"{s['razon_E']:.3f}, satura {s['satura']}")


if __name__ == '__main__':
    main()
