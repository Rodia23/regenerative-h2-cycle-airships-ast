# [EN] Delivery in consecutive releases (Section 5.3, Figures S5-S6).
# =============================================================================
# exp20_reparto_sueltas.py — EXP-20: repartir la carga en sueltas con el freno (autor, 2026-10-12). En el paper: la
# sección 5.3 (figura 2x2 de las cuatro naves); en el suplementario, la familia universal. Está en correr_todo.py.
#
# Con el freno, la nave lleva una carga de negocio N y la suelta en sueltas de tamaño s, reusando el mismo equipo.
# Para que cada suelta se recupere en t hace falta el rechazo P = Q s / t, y el equipo pesa
#     compresor de H2 (unidades de catálogo: no se achica con la suelta, decisión del autor; crece si una no alcanza)
#     + s (b' + a Q / t),   b' = b - compresor por tonelada (tanques y PEM, que sí escalan con la suelta);
# la carga útil C lo contiene junto con la carga de negocio:  N = C - compresor - s (b' + a Q / t),  con N >= s.
# (D, 2026-10-13) El agua retenida de cada suelta queda a bordo hasta la regeneración: el tanque de agua
# (fig_presupuesto_masa, w t por tonelada liberada) guarda la de TODA la carga entregada, no la de una suelta, así
# que sale de b' y va con N:  N (1 + w) = C - compresor - s (b' + a Q / t). En la suelta única (N = s) vuelve a
# ser b' + w, la b de la carta de diseño.
# El borde N = s es la suelta única: la recuperación más corta para ese tamaño de suelta. Con s = 1 t el compresor fijo
# y el proporcional pesan lo mismo, así que ahí coincide con la carta de diseño (EXP-18) y la barra HDM.
# La familia universal del suplementario toma todo proporcional, como la Fig. 4: N / C = 1 - f (b + a Q / t), f = s / C,
# la misma para toda nave (cada una entra por su f; el compresor fijo la baja en compresor / C).
# (2026-10-13) El compresor de H2 va en unidades de catálogo (fig_presupuesto_masa.compresor_t): una por suelta mientras
# alcance a comprimir, al ritmo del electrolizador, el H2 que la suelta regenera; más si la suelta lo pide. La familia
# universal usa el compresor proporcional a su capacidad (b lineal), el límite de las naves grandes.
# Hélices como freno: al soltar sujetan s g (sin restar lo que la quema absorbe sola: conservador) con la potencia ideal
# de helices_vectoriales.potencia; el motor entrega a la planta k P, k = [eta PCI - f_v (q0 - eta PCI)] / q0, con los
# eta y f_v de la potencia neta de la tabla de rutas (macro DoceNetaCiclo), ya descontado el compresor del escape
# (2026-10-12: k = [eta PCI - f_v (q0 - eta PCI + W) - W] / q0). Alcanzan si P_hel(s) <= k P, o sea
# t <= t_max(s) = k Q s / P_hel(s): dos bordes, pesimista (k mínimo, discos de 20 m2) y optimista (k máximo, 100 m2).
# Más lento que eso el motor, que va con el rechazo instalado, es muy chico: amarre (o una batería, no contada).
# Reparto: se carga todo en la base y se suelta en varios puntos; entre sueltas no se devuelve sustentación (el agua
# queda como lastre; la regeneración va después). Carga divisible.
# Salidas: resultados/figuras/exp20_reparto_naves.png (2x2), exp20_reparto_universal.png;
#          resultados/datos/exp20_reparto_sueltas.json; resultados/tablas/exp20_reparto_sueltas.txt
# =============================================================================
import sys

import numpy as np

import parametros as p
import rutas

ORDEN = ('Hind', 'FW', 'ZR', 'ZNT')                       # como las Figs. 3 y 5
SUELTAS = {'Hind': (0.25, 0.5, 1.0, 2.0, 5.0), 'FW': (0.5, 1.0, 2.0, 5.0, 10.0, 20.0),
           'ZR': (0.1, 0.25, 0.5, 1.0, 2.0), 'ZNT': (0.025, 0.05, 0.1, 0.2)}     # [t]; el NT, hasta lo que cabe con el compresor fijo (0.24 t)
FRACS = (0.01, 0.02, 0.05, 0.1, 0.2)                      # familia universal: suelta / carga útil
T_TABLA = (600.0, 1800.0, 3600.0, 86400.0)                # [s]
MARGEN_REF = {'nave': 'FW', 's_t': 1.0, 't_s': 600.0, 'x': 0.10}   # B2: margen de diseño del 10 % sobre la entrega de referencia
T_EJE = (60.0, 86400.0)                                   # [s] de 1 min a 24 h
ZONA = ('#ebe8e2', '#d3cfc6')                      # fondo de la Fig. 6: gris claro, gris medio
GRIS = ('#c4c1bb', '#a8a49d', '#8a867e', '#6b675f', '#4a4740', '#28261f')


def datos():
    import fig_presupuesto_masa as fpm
    import fig_esferoides as fe
    from balance_masa import f_caida_presion
    import helices_vectoriales as hv
    d = fpm.datos_apilado()
    a, b, Q = d['t_por_MW'], d['fijas_t'], d['calor_t_J']                      # [t/MW], [t/t], [J/t]
    aQ = a * Q / 1e6                                                            # [s]
    # (2026-10-13) el compresor de H2 en unidades de catálogo (fpm.compresor_t): una por suelta hasta lo que una
    # comprime al ritmo del electrolizador, más si la suelta lo pide; comp es una unidad
    cu, w = d['compresor'], d['agua_t']                                         # w: tanque de agua [t por t]
    b1 = d['fijas_sin_comp_t'] - w                                              # b': tanques y PEM [t por t]
    comp = cu['kg'] / 1e3                                                       # [t] una unidad
    q0 = next(r for r in rutas.leer_datos('exp12_helio')['rutas'] if 'phi = 0.99' in r['ruta'])['q_Jkg']
    # (2026-10-12) el motor también paga el compresor del escape (14 atm), que sale como calor
    Wc = rutas.leer_datos('exp06_requisitos_minimos')['W_comp_MJkg'] * 1e6
    ks = [(eta * p.Q_comb_H2 - fv * (q0 - eta * p.Q_comb_H2 + Wc) - Wc) / q0
          for eta in p.eta_motor_H2 for fv in (min(f_caida_presion()), max(p.perdida_rad))]
    k_min, k_max = min(ks), max(ks)
    A_min, A_max = min(p.freno_A_discos), max(p.freno_A_discos)

    def t_hel(s, k, A):
        return k * Q * s / hv.potencia(s * 1000.0 * p.g, A)                   # [s]

    naves = {}
    for c in ORDEN:
        nom = next(n for n, v in p.PLATAFORMAS.items() if v['corto'] == c)
        C = p.PLATAFORMAS[nom]['carga_max'] / 1e3
        filas = []
        for s in SUELTAS[c]:
            cs = fpm.compresor_t(s, cu)                                         # [t] las unidades de esta suelta
            disp = C - cs - s * (1.0 + b1 + w)
            t_min = s * aQ / disp if disp > 0 else None                         # una sola suelta (N = s)
            fila = {'s_t': s, 'comp_t': cs, 't_min_s': t_min, 't_hel_pes_s': t_hel(s, k_min, A_min),
                    't_hel_opt_s': t_hel(s, k_max, A_max), 'P_hel_pes_W': hv.potencia(s * 1000.0 * p.g, A_min),
                    'N_t': {f'{x:g}': (C - cs - s * (b1 + aQ / x)) / (1.0 + w) for x in T_TABLA},
                    'N_max_t': (C - cs - s * b1) / (1.0 + w)}
            filas.append(fila)
        naves[c] = {'nombre': fe.CORTO[c], 'carga_t': C, 'color': p.PLATAFORMAS[nom]['color'], 'sueltas': filas,
                    'unica_max_t': fpm.suelta_unica_max(C, b1 + w, 0.0, cu)}     # la suelta única más grande (t muy largo)
    # (2026-10-14, autor: la cota pesimista de la entrega) la misma entrega en sueltas de 1 t con el sistema ACC entero en
    # vez de sus haces (fpm.datos_apilado, t_por_MW_acc): solo cambia a; b', b_w, Q y el compresor de H2 son los mismos
    aQ_acc = d['t_por_MW_acc'] * Q / 1e6                                       # [s]
    acc = {'a_t_por_MW': d['t_por_MW_acc'], 'aQ_s': aQ_acc, 'naves': {}}
    for c in ORDEN:
        C = naves[c]['carga_t']
        cs = fpm.compresor_t(1.0, cu)
        disp = C - cs - (1.0 + b1 + w)
        acc['naves'][c] = {'t_min_s': aQ_acc / disp if disp > 0 else None,
                           'N_t': {f'{x:g}': (C - cs - (b1 + aQ_acc / x)) / (1.0 + w) for x in T_TABLA}}
    # B2 (2026-10-13; exploración en exploracion/B2_margen_masa): un margen de diseño x sobre el caso de referencia
    # (LCA60T, sueltas de 1 t de 10 min). Masa: todo el equipo x (1 + x) (a, b', b_w y el compresor); calor: solo
    # a Q x (1 + x) (más calor por tonelada con el mismo rechazo instalado, o rechazo instalado de sobra)
    C_fw = p.PLATAFORMAS[next(n for n, v in p.PLATAFORMAS.items() if v['corto'] == 'FW')]['carga_max'] / 1e3
    cs1 = fpm.compresor_t(MARGEN_REF['s_t'], cu)

    def N_margen(fm, fq):
        s_, t_ = MARGEN_REF['s_t'], MARGEN_REF['t_s']
        return (C_fw - fm * cs1 - s_ * (fm * b1 + fm * fq * aQ / t_)) / (1.0 + fm * w)
    margen = {**MARGEN_REF, 'N_t': N_margen(1.0, 1.0),
              'N_masa_t': N_margen(1.0 + MARGEN_REF['x'], 1.0), 'N_calor_t': N_margen(1.0, 1.0 + MARGEN_REF['x'])}
    # (2026-10-13) la entrega a 1 atm, sin el compresor del escape (auditoría: los dos casos como dato): a y Q de 1 atm
    # de fpm.datos_apilado; b', b_w y el compresor de H2 son los mismos (cajas fijas por tonelada liberada; b_w es el
    # de la clase liviana, la diferencia de agua entre los dos casos le cambia menos de 0.2 kg por tonelada)
    u1 = d['una_atm']
    una_atm = {'a_t_por_MW': u1['t_por_MW'], 'Q_J_por_t': u1['calor_t_J'], 'aQ_s': u1['t_por_MW'] * u1['calor_t_J'] / 1e6}
    return {'una_atm': una_atm, 'acc_entero': acc,
            'a_t_por_MW': a, 'b_t_por_t': b, 'b_sin_comp_t_por_t': b1, 'w_t_por_t': w, 'comp_t': comp, 'compresor': cu,
            'b_lineal_t_por_t': b1 + comp * cu['kgph_por_t'] / cu['cap_kgph'], 'Q_J_por_t': Q, 'aQ_s': aQ,
            'k': [k_min, k_max], 'A_m2': [A_min, A_max], 'naves': naves, 'margen': margen}


def _ejes_x(ax, angosto=False):
    from matplotlib.ticker import FixedLocator, NullFormatter
    ax.set_xscale('log')
    ax.set_xlim(T_EJE[0] / 60, T_EJE[1] / 60)
    pos = [1, 10, 60, 600, 1440] if not angosto else [1, 10, 60, 1440]
    lab = ['1 min', '10 min', '1 h', '10 h', '24 h'] if not angosto else ['1 min', '10 min', '1 h', '24 h']
    ax.xaxis.set_major_locator(FixedLocator(pos))
    ax.set_xticklabels(lab)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(labelsize=7.5)
    ax.grid(True, which='major', alpha=0.25)


def _tono(color, w):
    """El color de la nave aclarado (w > 0, hacia blanco) u oscurecido (w < 0, hacia negro)."""
    from matplotlib.colors import to_rgb
    r = np.array(to_rgb(color))
    return tuple(r + (1.0 - r) * w) if w >= 0 else tuple(r * (1.0 + w))


def _rot(s):
    return f'{s:g} t'                                      # en toneladas (autor, 2026-10-12)


def figura_naves(d):
    """Fig. 6 del paper (autor, 2026-10-12): curvas continuas del mismo grosor en el color de la nave; el freno como
    zonas de fondo: cada punto (t, N) sobre la suelta única define un tamaño de suelta s = (C - comp - N) / (b' + aQ/t),
    que se clasifica con el t máximo de las hélices para ese s (blanco: hélices o amarre; gris claro: amarre, o hélices
    solo con discos grandes y motor eficiente; gris medio: solo amarre o batería); la curva de 1 t, con sus valores a
    10 y 30 min."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    import helices_vectoriales as hv

    import fig_presupuesto_masa as fpm
    aQ, b1, comp, Q = d['aQ_s'], d['b_sin_comp_t_por_t'], d['comp_t'], d['Q_J_por_t']
    w = d['w_t_por_t']
    cu = d['compresor']
    cap_t = cu['cap_kgph'] / cu['kgph_por_t']                                   # [t] lo que comprime una unidad
    (k0, k1), (A0, A1) = d['k'], d['A_m2']
    t = np.geomspace(*T_EJE, 600)
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 6.2))
    for ax, c in zip(axs.flat, ORDEN):
        n = d['naves'][c]
        C, col = n['carga_t'], n['color']
        T, Ny = np.meshgrid(t, np.linspace(0.0, 1.03 * C, 400))
        # (2026-10-13) el tamaño de suelta con el compresor en unidades: con n unidades, s = (C - n comp - N) / (b' +
        # aQ/T), válido si cae en el tramo de n unidades
        s = np.full(T.shape, np.nan)
        for nu in range(1, int(np.ceil(C / cap_t)) + 2):
            sn = (C - nu * comp - Ny * (1.0 + w)) / (b1 + aQ / T)
            tramo = np.isnan(s) & (sn <= nu * cap_t) & ((sn > (nu - 1) * cap_t) | (nu == 1))
            s = np.where(tramo, sn, s)
        s = np.nan_to_num(s, nan=-1.0)
        valido = (s > 0) & (Ny >= s)
        sc = np.clip(s, 1e-9, None) * 1000.0 * p.g
        tp = k0 * Q * s / hv.potencia(sc, A0)
        to = k1 * Q * s / hv.potencia(sc, A1)
        clase = np.where(~valido, np.nan, np.where(T <= tp, 0, np.where(T <= to, 1, 2)))
        cs = ax.contourf(T / 60, Ny, clase, levels=[-0.5, 0.5, 1.5, 2.5], colors=['#ffffff', ZONA[0], ZONA[1]],
                         zorder=0)
        cs.set_edgecolor('face')                                 # sin rendijas entre los rellenos en el PDF
        for f in n['sueltas']:
            sv = f['s_t']
            N = (C - fpm.compresor_t(sv, cu) - sv * (b1 + aQ / t)) / (1.0 + w)
            ok = N >= sv
            if not ok.any():
                continue
            ax.plot(t[ok] / 60, N[ok], color=col, lw=1.7, zorder=3, solid_capstyle='round')
            i0 = int(np.argmax(ok))
            ax.annotate(_rot(sv), xy=(t[i0] / 60, N[i0]), xytext=(-2, 3), textcoords='offset points', fontsize=7,
                        ha='right', va='bottom', color='#222222')
            if abs(sv - 1.0) < 1e-9:                      # el hallazgo: la suelta de 1 t a 10 y 30 min
                for tm in (600.0, 1800.0):
                    Nm = (C - fpm.compresor_t(sv, cu) - sv * (b1 + aQ / tm)) / (1.0 + w)
                    if Nm >= sv:
                        ax.plot(tm / 60, Nm, 'o', ms=6, color=col, mec='white', mew=1.2, zorder=5)
                        ax.annotate(f'{Nm:.0f} t' if Nm >= 10 else f'{Nm:.1f} t', xy=(tm / 60, Nm), xytext=(6, -10),
                                    textcoords='offset points', fontsize=7.5, fontweight='bold', color='#222222',
                                    bbox=dict(facecolor='white', edgecolor='none', pad=0.8, alpha=0.85), zorder=6)
        ax.plot(t / 60, [fpm.suelta_unica_max(C, b1 + w, aQ / x, cu) for x in t], color='#c0392b', lw=1.3,
                ls=(0, (4, 2)), zorder=4)
        ax.set_ylim(0, 1.03 * C)
        ax.set_title(n['nombre'], fontsize=9.5, fontweight='bold', loc='left', pad=3, color='#111111')
        ax.text(0.98, 0.04, f'payload {C:.3g} t', transform=ax.transAxes, ha='right', fontsize=7.5, color='#555555')
        _ejes_x(ax, angosto=True)
        ax.set_axisbelow(False)
    for ax in axs[1]:
        ax.set_xlabel('Recovery time of each release', fontsize=8.5)
    for ax in axs[:, 0]:
        ax.set_ylabel('Delivered cargo (t)', fontsize=8.5)
    fig.legend(handles=[Patch(facecolor='#ffffff', edgecolor='#999999', lw=0.6, label='propellers or tether'),
                        Patch(facecolor=ZONA[0], edgecolor='#999999', lw=0.6,
                              label='tether, or propellers with large disks and an efficient engine'),
                        Patch(facecolor=ZONA[1], edgecolor='#999999', lw=0.6, label="tether, or power from the airship's own plant or a battery"),
                        Line2D([], [], color='#c0392b', lw=1.3, ls=(0, (4, 2)), label='single release')],
               loc='lower center', ncol=2, fontsize=7.5, frameon=False, handlelength=1.6)
    plot_style.normalizar_paper(fig)
    fig.tight_layout(rect=(0, 0.07, 1, 1), h_pad=1.0, w_pad=1.2)
    plot_style.guardar_paper(fig, rutas.figura('exp20_reparto_naves'))
    plt.close(fig)


def figura_universal(d):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401
    from matplotlib.lines import Line2D

    # (2026-10-13) con el compresor de H2 en unidades de catálogo la familia ya no es exacta en toda nave: va con el
    # compresor proporcional a su capacidad (b lineal), el límite de las naves grandes, donde las sueltas piden varias
    # unidades; en las chicas la unidad entera pesa más
    aQ, b, w = d['aQ_s'], d['b_lineal_t_por_t'], d['w_t_por_t']
    t = np.geomspace(*T_EJE, 1200)
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    tonos = ['#0072B2'] * len(FRACS)
    for f, g in zip(FRACS, tonos):
        y = 100 * (1 - f * (b + aQ / t)) / (1.0 + w)
        ok = y >= 100 * f
        ax.plot(t[ok] / 60, y[ok], color=g, lw=1.6, zorder=3)
        i0 = int(np.argmax(ok))
        ax.annotate(f'{100 * f:g} %', xy=(t[i0] / 60, y[i0]), xytext=(-3, 5), textcoords='offset points',
                    fontsize=7.5, ha='right', va='bottom', color='#222222')
        for n, mk in ((2, 's'), (5, 'D'), (10, '^')):
            den = 1 - f * (b + n * (1.0 + w))
            if den > 0 and T_EJE[0] <= f * aQ / den <= T_EJE[1]:
                ax.plot(f * aQ / den / 60, 100 * n * f, mk, ms=4.5, color=g, mec='#222222', mew=0.6, zorder=5)
    ax.plot(t / 60, 100 / (1 + b + w + aQ / t), color='#c0392b', lw=1.2, ls=(0, (4, 2)), zorder=4)
    _ejes_x(ax)
    ax.set_ylim(0, 100)
    ax.set_xlabel('Recovery time of each release', fontsize=9)
    ax.set_ylabel('Delivered cargo / payload (%)', fontsize=9)
    ax.legend(handles=[Line2D([], [], color='#0072B2', lw=1.6, label='releases of the size marked (% of the payload)'),
                       Line2D([], [], color='#c0392b', lw=1.2, ls=(0, (4, 2)), label='single release (Eq. 10)'),
                       Line2D([], [], ls='none', marker='s', ms=4.5, color=_tono('#0072B2', 0.3), mec='#222222', label='2 releases'),
                       Line2D([], [], ls='none', marker='D', ms=4.5, color=_tono('#0072B2', 0.3), mec='#222222', label='5 releases'),
                       Line2D([], [], ls='none', marker='^', ms=4.5, color=_tono('#0072B2', 0.3), mec='#222222', label='10 releases')],
              loc='lower right', fontsize=7, frameon=False)
    fig.tight_layout()
    plot_style.guardar_paper(fig, rutas.figura('exp20_reparto_universal'))
    plt.close(fig)


def informe(d):
    hm = lambda s: f'{s / 60:.1f} min' if s < 3600 else f'{s / 3600:.1f} h'
    L = ['EXP-20 — REPARTIR LA CARGA EN SUELTAS CON EL FRENO (generado por modelo/exp20_reparto_sueltas.py)',
         f"N (1 + w) = C - compresor(s) - s (b' + a Q / t); a Q = {d['aQ_s'] / 3600:.2f} h, b' = {d['b_sin_comp_t_por_t']:.3f} t/t, "
         f"w = {1e3 * d['w_t_por_t']:.2f} kg/t (tanque de agua), "
         f"compresor de H2 en unidades de {d['comp_t']:.3f} t, una por "
         f"{d['compresor']['cap_kgph'] / d['compresor']['kgph_por_t']:.2f} t de suelta; k = {d['k'][0]:.3f}-{d['k'][1]:.3f}",
         '']
    for c in ORDEN:
        n = d['naves'][c]
        L.append(f"{n['nombre']} (carga útil {n['carga_t']:.2f} t)")
        for f in n['sueltas']:
            L.append(f"  suelta {_rot(f['s_t']):>7}: desde {hm(f['t_min_s']) if f['t_min_s'] else '—':>9}; "
                     + ', '.join(f"{hm(float(x))}: {max(v, 0):.2f} t" for x, v in f['N_t'].items())
                     + f"; compresor {f['comp_t']:.2f} t; máx {f['N_max_t']:.2f} t; hélices hasta {hm(f['t_hel_pes_s'])} / {hm(f['t_hel_opt_s'])}")
    return '\n'.join(L)


def main():
    d = datos()
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    figura_naves(d)
    figura_universal(d)
    rutas.guardar_datos('exp20_reparto_sueltas', d)
    txt = informe(d)
    with open(rutas.tabla('exp20_reparto_sueltas').replace('.tex', '.txt'), 'w', encoding='utf-8') as fh:
        fh.write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
