# [EN] Stacked mass budget without and with a hold-down device (Figure 5).
# =============================================================================
# fig_presupuesto_masa.py — El presupuesto de masa (tabla 5 del paper 1) en barras (propuesta, 2026-10-11)
#
# Las mismas cajas y cifras de la tabla (req_masa_rechazo), de mayor a menor. La tabla mezcla unidades (las que
# escalan con el flujo en t por kg/s de H2, las fijas en t por tonelada liberada): aquí todas en toneladas para el caso
# base, 1 t en el LCA60T, las que escalan por la tasa de la cota ideal (ṁ*), para compararlas con su carga útil.
# Dos versiones para elegir (autor, 2026-10-11): en % de la carga útil (barras, eje lineal; llena hasta el mínimo y
# clara hasta el máximo) y en t con eje logarítmico (rangos con puntos, no barras). La unidad de rechazo no se elige:
# es lo que deja el presupuesto con la batería más liviana.
# EN EL PAPER (autor, 2026-10-11): presupuesto_masa_apilado_log.png es la figura de la sección 5.2 (las cuatro naves,
# sin freno, con 10 MW y con el rechazo más grande que cabe; la tabla 5 pasó al suplemento). Macros Pila* en
# exportar_numeros; chequeos en verificacion (grupo I). Las demás versiones quedan como alternativas.
# Con el freno (autor, 2026-10-11): el mismo presupuesto con 10 MW de rechazo instalado (P_ref_masa), otra paleta
# para no confundirlo con la cota ideal. La tasa es P / q (sección 5.2), la misma en todas las naves: cambia solo la
# carga útil, así que los cuatro paneles van con la misma escala. Motor (que mueve los ventiladores: sin batería,
# cambio 7), haces de ACC hasta el sistema ACC entero, compresor del escape a 14 atm (cambio 6), las fijas y el freno
# (las hélices que ya lleva la nave; el cable queda en tierra). Total con los haces: el motor más liviano y el más pesado.
# Salidas: resultados/figuras/presupuesto_masa_pct.png, _log.png, _naves.png (las cuatro naves, en t),
#          _freno_naves.png (con el freno, 10 MW); resultados/datos/fig_presupuesto_masa.json.
# =============================================================================
import math
import sys

import rutas

AZUL, NARANJA = '#2a78d6', '#eb6834'       # paleta categórica de referencia (validada: 2 grupos)


# ── Compresor de H2 en unidades de catálogo (autor, 2026-10-13) ──────────────────────────────────────────────────
# Una unidad (parametros.COMPRESOR_REF, la más chica de 6000 psi que alcanza el caudal de una tonelada) comprime el H2
# que el electrolizador regenera al ritmo de la suelta en t_regeneracion_h. Una suelta de m toneladas lleva las
# unidades que pide su caudal, al menos una: no se achica con la suelta (decisión del autor en EXP-20) y crece cuando
# una no alcanza. Antes: una unidad por tonelada en la carta de diseño (EXP-18) y una sola fija en el reparto (EXP-20,
# EXP-22), que dejaban techos distintos para la misma suelta única (24 y 31 t en el LCA60T). Lo usan la pila de masas,
# EXP-18, EXP-20 y EXP-22.
def compresor_unidad():
    import parametros as p
    return {'kg': p.COMPRESOR_REF['masa_kg'][0], 'cap_kgph': min(p.COMPRESOR_REF['flujo_kgph']),
            'kgph_por_t': rutas.leer_datos('balance_masa')['recuperado_por_t_kg'] / p.t_regeneracion_h}


def unidades_compresor(m_t, cu):
    return max(1, math.ceil(m_t * cu['kgph_por_t'] / cu['cap_kgph'] - 1e-9))


def compresor_t(m_t, cu):
    """Masa del compresor de H2 para una suelta de m_t toneladas [t]."""
    return unidades_compresor(m_t, cu) * cu['kg'] / 1e3


def suelta_unica_max(C, bp, x, cu):
    """La suelta única más grande que cabe en la carga útil C [t]: m (1 + b' + x) + compresor(m) <= C, con x = a Q / t
    (x = 0: recuperación muy larga, el techo) y b' las fijas sin el compresor [t por t]. Con n unidades caben hasta
    n veces lo que una comprime; se toma el mejor n."""
    mu, cap_t = cu['kg'] / 1e3, cu['cap_kgph'] / cu['kgph_por_t']
    mejor, n = 0.0, 1
    while C - n * mu > 0.0:
        m = (C - n * mu) / (1.0 + bp + x)
        mejor = max(mejor, min(m, n * cap_t))
        if m <= n * cap_t:
            break
        n += 1
    return mejor
AGUA, VIOLETA = '#1baf7a', '#4a3aa7'       # la del freno, distinta (validada; el agua pide rótulos: los lleva)


def datos():
    r = rutas.leer_datos('req_masa_rechazo')['requisitos']
    m = r['mdot_kgps']
    fj = r['fijas_kg']
    cajas = [  # (nombre, clase, rango [t], grupo)
        ('Energy conversion unit', 'H$_2$ engine', [x * m for x in r['conversion_t_por_kgps']], 'escala'),
        ('Energy buffer', 'Li-ion pack', [x * m for x in r['buffer_quemador_t_por_kgps']], 'escala'),
        ('Heat-rejection unit', 'budget remainder', [0.0, r['queda_rechazo_t_por_kgps'] * m], 'resto'),
        ('H$_2$ tanks', 'Type IV', [x / 1e3 for x in fj['tanques_por_t']], 'fija'),
        ('H$_2$ compressor', 'electric', [x / 1e3 for x in fj['compresor']], 'fija'),
        ('Regeneration unit', 'PEM stack', [x / 1e3 for x in fj['regeneracion']], 'fija'),
        ('Water tank', 'fabric', [x / 1e3 for x in fj['tanque_agua']], 'fija'),
    ]
    cajas.sort(key=lambda c: -c[2][1])
    return {'mdot_kgps': m, 'carga_t': r['carga_t'], 'fijas_t': r['fijas_t'],
            'cajas': [{'nombre': n, 'clase': k, 't': t, 'grupo': g} for n, k, t, g in cajas]}


def _fmt(t):
    return f'{t:.0f}' if t >= 10 else (f'{t:.1f}' if t >= 1 else f'{t:.2f}')


def _tmin(t):
    """Minutos con un decimal bajo 100 (10.6 min, como en el paper)."""
    return f'{t:.1f}' if t < 100 else f'{t:.0f}'


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import plot_style  # noqa: F401
    return plt


def _ejes_comunes(ax, r, x_carga, txt_carga, leyenda, donde='lower right'):
    ax.set_yticks([-i for i in range(len(r['cajas']))])
    ax.set_yticklabels([f"{c['nombre']}\n({c['clase']})" for c in r['cajas']], fontsize=9)
    ax.tick_params(axis='y', which='both', left=False, right=False)
    ax.tick_params(top=False)
    ax.set_ylim(-len(r['cajas']) + 0.4, 0.6)
    ax.grid(axis='y', visible=False)
    # la carga útil: la referencia contra la que se mide todo
    ax.axvline(x_carga, color='#52514e', lw=1.0, ls=(0, (4, 3)), zorder=0)
    ax.annotate(txt_carga, xy=(x_carga, 0.42), xytext=(4, 0), textcoords='offset points', fontsize=9,
                color='#52514e', va='center', ha='left')
    ax.legend(handles=leyenda, loc=donde, fontsize=8.5, frameon=False)


def _leyenda(r, marcador):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    txt = (f"Scales with the burn rate ($\\times$ {r['mdot_kgps']:.2f} kg/s)", 'Fixed, per tonne released')
    if marcador:
        return [Line2D([], [], marker='o', ls='-', lw=2.5, color=c, label=t) for c, t in zip((AZUL, NARANJA), txt)]
    return [Patch(color=c, label=t) for c, t in zip((AZUL, NARANJA), txt)]


def figura_pct(r):
    """Barras en % de la carga útil (eje lineal): la carga útil es el 100 %."""
    plt = _plt()
    k = 100.0 / r['carga_t']
    fig, ax = plt.subplots(figsize=(7.4, 3.3))
    for i, c in enumerate(r['cajas']):
        y = -i
        lo, hi = c['t'][0] * k, c['t'][1] * k
        col = NARANJA if c['grupo'] == 'fija' else AZUL
        if c['grupo'] == 'resto':
            ax.barh(y, hi, height=0.62, facecolor='white', edgecolor=AZUL, hatch='////', lw=1.0)
            txt = rf"$\leq$ {_fmt(hi)} %"
        else:
            ax.barh(y, lo, height=0.62, color=col, lw=0)
            if hi > lo * 1.001:
                ax.barh(y, hi - lo, left=lo, height=0.62, color=col, alpha=0.32, lw=0)
            if hi < 3.0:                 # las de menos del 3 % no se ven a esta escala: un punto al final
                ax.plot(hi, y, 'o', ms=5, color=col, mec='white', mew=1.0, clip_on=False, zorder=3)
            txt = f'{_fmt(lo)} %' if hi <= lo * 1.001 else f'{_fmt(lo)}–{_fmt(hi)} %'
        ax.text(hi + 10.0, y, txt, va='center', ha='left', fontsize=9, color='#222222')
    _ejes_comunes(ax, r, 100.0, 'LCA60T payload', _leyenda(r, False))
    ax.set_xlim(0.0, 1.2 * max(c['t'][1] for c in r['cajas']) * k)
    ax.set_xlabel('Share of the LCA60T payload, 1-tonne release (%)')
    fig.tight_layout()
    fig.savefig(rutas.figura('presupuesto_masa_pct'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def figura_log(r):
    """Rangos en t sobre eje logarítmico: un segmento del mínimo al máximo con puntos (no barras: en log el largo de
    una barra no es proporcional al valor y su inicio es arbitrario)."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7.4, 3.3))
    for i, c in enumerate(r['cajas']):
        y = -i
        lo, hi = c['t']
        col = NARANJA if c['grupo'] == 'fija' else AZUL
        if c['grupo'] == 'resto':        # "a lo sumo": rombo hueco con flecha hacia abajo
            ax.annotate('', xy=(hi / 4.0, y), xytext=(hi, y),
                        arrowprops=dict(arrowstyle='-|>', color=AZUL, lw=1.2, mutation_scale=9))
            ax.plot(hi, y, 'D', ms=7, mfc='white', mec=AZUL, mew=1.5, zorder=3)
            txt = rf"$\leq$ {_fmt(hi)} t"
        else:
            xs = [lo, hi] if hi > 1.15 * lo else [hi]   # un rango que no se distingue en log: un punto (el texto lo da)
            if len(xs) == 2:
                ax.plot(xs, [y, y], color=col, lw=2.5, alpha=0.45, solid_capstyle='round', zorder=2)
            ax.plot(xs, [y] * len(xs), 'o', ms=7, color=col, mec='white', mew=1.2, zorder=3)
            txt = f'{_fmt(lo)} t' if hi <= lo * 1.001 else f'{_fmt(lo)}–{_fmt(hi)} t'
        ax.annotate(txt, xy=(hi, y), xytext=(9, 0), textcoords='offset points', va='center', ha='left', fontsize=9,
                    color='#222222')
    ax.set_xscale('log')
    _ejes_comunes(ax, r, r['carga_t'], f"LCA60T payload, {r['carga_t']:.0f} t", _leyenda(r, True), 'upper left')
    ax.set_xlim(0.02, 3000.0)
    ax.set_xticks([0.1, 1, 10, 100, 1000])
    ax.set_xticklabels(['0.1', '1', '10', '100', '1000'])
    ax.set_xlabel('Mass for a 1-tonne release on the LCA60T (t)')
    fig.tight_layout()
    fig.savefig(rutas.figura('presupuesto_masa_log'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def datos_naves():
    """El mismo presupuesto para las cuatro naves (autor, 2026-10-11): 1 t liberada en cada una, con su tasa de la cota
    ideal (minimax.resolver, la de la tabla 5) y su carga útil. Las masas por kg/s y las fijas por tonelada son las
    mismas; al condensador le queda la carga útil menos las fijas y la batería más liviana, si queda algo."""
    import parametros as p
    import minimax as mm
    import fig_esferoides as fe
    r = rutas.leer_datos('req_masa_rechazo')['requisitos']
    fj = r['fijas_kg']
    naves = []
    for c in fe.ORDEN:
        nom = next(n for n, d in p.PLATAFORMAS.items() if d['corto'] == c)
        m = mm.resolver(nom, 1000.0, p.h_max_precision, p.D_lim_oper)['dm_star']
        carga = p.PLATAFORMAS[nom]['carga_max'] / 1e3
        bat = [x * m for x in r['buffer_quemador_t_por_kgps']]
        naves.append({'corto': c, 'nombre': fe.CORTO[c], 'carga_t': carga, 'mdot_kgps': m,
                      'queda_t': carga - r['fijas_t'] - bat[0], 'cajas': [
                          ('Conversion (H$_2$ engine)', [x * m for x in r['conversion_t_por_kgps']], 'escala'),
                          ('Buffer (Li-ion pack)', bat, 'escala'),
                          ('Heat rejection (remainder)', None, 'resto'),
                          ('H$_2$ tanks', [x / 1e3 for x in fj['tanques_por_t']], 'fija'),
                          ('H$_2$ compressor', [x / 1e3 for x in fj['compresor']], 'fija'),
                          ('Regeneration (PEM)', [x / 1e3 for x in fj['regeneracion']], 'fija'),
                          ('Water tank', [x / 1e3 for x in fj['tanque_agua']], 'fija')]})
    if abs(naves[fe.ORDEN.index('FW')]['mdot_kgps'] - r['mdot_kgps']) > 1e-9 * r['mdot_kgps']:
        sys.exit('la tasa del LCA60T no es la de la tabla 5')
    return naves


def figura_naves(naves):
    """Una fila por nave, en t y eje lineal como la primera versión; cada fila con su escala y su carga útil. Solo la
    clase más liviana de cada caja (autor, 2026-10-11: el rango hasta la más pesada va al texto)."""
    plt = _plt()
    fig, axes = plt.subplots(len(naves), 1, figsize=(7.4, 10.0))
    for ax, nv in zip(axes, naves):
        x_max = 1.22 * max(max(t[0] for _, t, _ in nv['cajas'] if t), nv['carga_t'])
        for i, (nombre, t, grupo) in enumerate(nv['cajas']):
            y = -i
            col = NARANJA if grupo == 'fija' else AZUL
            if grupo == 'resto':
                if nv['queda_t'] > 0:
                    ax.barh(y, nv['queda_t'], height=0.62, facecolor='white', edgecolor=AZUL, hatch='////', lw=1.0)
                    txt, x_txt = rf"$\leq$ {_fmt(nv['queda_t'])} t", nv['queda_t']
                else:
                    txt, x_txt = f"none left: short by {_fmt(-nv['queda_t'])} t", 0.0
            else:
                lo = t[0]
                ax.barh(y, lo, height=0.62, color=col, lw=0)
                if lo < 0.005 * x_max:   # invisibles a esta escala: un punto al final
                    ax.plot(lo, y, 'o', ms=4.5, color=col, mec='white', mew=0.9, clip_on=False, zorder=3)
                txt, x_txt = f'{_fmt(lo)} t', lo
            ax.annotate(txt, xy=(x_txt, y), xytext=(6, 0), textcoords='offset points', va='center', ha='left',
                        fontsize=8.5, color='#222222', zorder=4,     # fondo blanco: la línea de la carga útil los cruza
                        bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
        ax.set_yticks([-i for i in range(len(nv['cajas']))])
        ax.set_yticklabels([n for n, _, _ in nv['cajas']], fontsize=8.5)
        ax.tick_params(axis='y', which='both', left=False, right=False)
        ax.tick_params(top=False)
        ax.set_ylim(-len(nv['cajas']) + 0.4, 0.75)
        ax.grid(axis='y', visible=False)
        ax.set_xlim(0.0, x_max)
        ax.axvline(nv['carga_t'], color='#52514e', lw=1.0, ls=(0, (4, 3)), zorder=0)
        ax.annotate(f"payload, {_fmt(nv['carga_t'])} t", xy=(nv['carga_t'], 0.55), xytext=(4, 0),
                    textcoords='offset points', fontsize=8.5, color='#52514e', va='center', ha='left')
        ax.set_title(f"{nv['nombre']}: 1-tonne release, {nv['mdot_kgps']:.2f} kg/s", fontsize=10, loc='left')
    axes[-1].set_xlabel('Mass of the lightest class (t)')
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=AZUL, label='Scales with the burn rate'),
                        Patch(color=NARANJA, label='Fixed, per tonne released')],
               loc='upper center', ncol=2, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(rutas.figura('presupuesto_masa_naves'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def _compresor_escape_t(dm):
    """Masa del compresor del escape a la presión de phi objetivo, escalada con el paquete del SGT-800 [t, mín.-máx.]:
    la misma fórmula que exportar_numeros (macro FrenoCompMasaDiez: trabajo isentrópico de una etapa; la eficiencia se
    cancela)."""
    import parametros as p
    from exp05_condensacion_tiempo_real import p_sat, N_INERTE, AIRE_ESTEQ
    r6 = rutas.leer_datos('exp06_requisitos_minimos')
    n = N_INERTE / p.M_H2
    nv = n * p_sat(p.T_salida_cond) / (p.P_atmos - p_sat(p.T_salida_cond))
    C = (AIRE_ESTEQ - p.ratio_O2_H2) * p.cp_N2 + nv * p.M_H2O * p.cp_vapor
    ws = C * (p.T_salida_cond + 273.15) * (r6['P_req_atm'] ** ((n + nv) * p.R_univ / C) - 1.0)
    s8, ka = p.SGT800, (p.gamma_aire - 1.0) / p.gamma_aire
    w8 = [(me - hr * P / 3.6 / p.PCI_CH4 * 1e3) * p.cp_aire * (p.T_ISO_C + 273.15) * (pr ** ka - 1.0)
          for P, pr, me, hr in zip(s8['P_MWe'], s8['PR'], s8['m_escape_kgps'], s8['heat_rate_kJkWh'])]
    return [s8['masa_kg'] * dm * ws / w / 1e3 for w in (max(w8), min(w8))]


def datos_freno():
    """El presupuesto con el freno y P_ref_masa de rechazo instalado, 1 t, en las cuatro naves."""
    import parametros as p
    import fig_esferoides as fe
    from exp05_condensacion_tiempo_real import q_rech
    from exp02_quema_centrada import L_ef
    rq = rutas.leer_datos('req_masa_rechazo')
    r, cand = rq['requisitos'], {c['clave']: c for c in rq['candidatos']}
    P, q = p.P_ref_masa, q_rech(0.0, p.phi_objetivo)
    dm = P / q                                       # la tasa con el freno (sección 5.2): no depende de la nave
    fj = r['fijas_kg']
    hx = [P / 1e3 / cand[k]['kWkg'][0] / 1e3 for k in ('ACC_haces', 'ACC_sistema')]
    motor = [x * dm for x in r['conversion_t_por_kgps']]
    comp = _compresor_escape_t(dm)
    fijas = [[x / 1e3 for x in fj[k]] for k in ('tanques_por_t', 'compresor', 'regeneracion', 'tanque_agua')]
    cajas = [('Heat rejection (ACC bundles)', hx, 'escala'),
             ('Conversion (H$_2$ engine)', motor, 'escala'),
             ('Exhaust compressor (14 atm)', comp, 'escala'),
             ('H$_2$ tanks', fijas[0], 'fija'), ('H$_2$ compressor', fijas[1], 'fija'),
             ('Regeneration (PEM)', fijas[2], 'fija'), ('Water tank', fijas[3], 'fija')]
    cajas.sort(key=lambda c: -c[1][0])               # la figura muestra la clase más liviana: de mayor a menor en ella
    total = [motor[0] + hx[0] + comp[0] + sum(f[0] for f in fijas), motor[1] + hx[0] + comp[1] + sum(f[1] for f in fijas)]
    naves = []
    for c in fe.ORDEN:
        d = next(x for x in p.PLATAFORMAS.values() if x['corto'] == c)
        naves.append({'corto': c, 'nombre': fe.CORTO[c], 'carga_t': d['carga_max'] / 1e3})
    return {'P_MW': P / 1e6, 'mdot_kgps': dm, 't_rec_min': (1000.0 / L_ef(p.phi_objetivo)) * q / P / 60.0,
            'cajas': cajas, 'total_t': total, 'naves': naves}


def figura_freno_naves(d):
    """Una fila por nave, en t y eje lineal; con el freno las barras son las mismas en todas: misma escala. Solo la
    clase más liviana (autor, 2026-10-11): aquí es optimista (el sistema ACC entero pesa 74 t, no 15): va al texto."""
    plt = _plt()
    from matplotlib.patches import Patch
    filas = d['cajas'] + [('Hold-down device', None, 'freno'), ('Total', d['total_t'], 'total')]
    x_max = 1.15 * max(max(c[1][0] for c in d['cajas']), max(nv['carga_t'] for nv in d['naves']))
    fig, axes = plt.subplots(len(d['naves']), 1, figsize=(7.4, 11.0), sharex=True)
    for ax, nv in zip(axes, d['naves']):
        for i, (nombre, t, grupo) in enumerate(filas):
            y = -i
            if grupo == 'freno':
                txt, x_txt = '0 t: the propellers already on board (a tether stays on the ground)', 0.0
            elif grupo == 'total':
                lo = t[0]
                ax.barh(y, lo, height=0.62, facecolor='white', edgecolor='#333333', hatch='////', lw=1.0)
                txt, x_txt = f'{_fmt(lo)} t', lo
            else:
                lo = t[0]
                col = VIOLETA if grupo == 'fija' else AGUA
                ax.barh(y, lo, height=0.62, color=col, lw=0)
                if lo < 0.005 * x_max:   # invisibles a esta escala: un punto al final
                    ax.plot(lo, y, 'o', ms=4.5, color=col, mec='white', mew=0.9, clip_on=False, zorder=3)
                txt, x_txt = f'{_fmt(lo)} t', lo
            ax.annotate(txt, xy=(x_txt, y), xytext=(6, 0), textcoords='offset points', va='center', ha='left',
                        fontsize=8.5, color='#222222', zorder=4,     # fondo blanco: la línea de la carga útil los cruza
                        bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
        ax.axhline(-len(filas) + 1.5, color='#c3c2b7', lw=0.7, zorder=0)       # separa el total
        ax.set_yticks([-i for i in range(len(filas))])
        ax.set_yticklabels([f[0] for f in filas], fontsize=8.5)
        ax.tick_params(axis='y', which='both', left=False, right=False)
        ax.tick_params(top=False)
        ax.set_ylim(-len(filas) + 0.4, 0.75)
        ax.grid(axis='y', visible=False)
        ax.set_xlim(0.0, x_max)
        ax.axvline(nv['carga_t'], color='#52514e', lw=1.0, ls=(0, (4, 3)), zorder=0)
        ax.annotate(f"payload, {_fmt(nv['carga_t'])} t", xy=(nv['carga_t'], 0.55), xytext=(4, 0),
                    textcoords='offset points', fontsize=8.5, color='#52514e', va='center', ha='left')
        ax.set_title(f"{nv['nombre']}: 1-tonne release with a hold-down device, {d['P_MW']:.0f} MW "
                     f"({d['mdot_kgps']:.3f} kg/s, {d['t_rec_min']:.1f} min)", fontsize=10, loc='left')
    axes[-1].set_xlabel('Mass of the lightest class (t)')
    fig.legend(handles=[Patch(color=AGUA, label=f"Scales with the heat rejection ({d['P_MW']:.0f} MW)"),
                        Patch(color=VIOLETA, label='Fixed, per tonne released'),
                        Patch(facecolor='white', edgecolor='#333333', hatch='////', label='Total')],
               loc='upper center', ncol=3, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(rutas.figura('presupuesto_masa_freno_naves'), dpi=300, bbox_inches='tight')
    plt.close(fig)


APILADO = (('Heat rejection (ACC bundles)', '#00869e'), ('Exhaust compressor (14 atm)', '#e06c5f'),
           ('Conversion (H$_2$ engine)', '#4a3aa7'), ('Fixed boxes (tanks, compressor, PEM)', '#cf9b3a'),
           ('Cargo released (1 t)', '#c9c4bb'))
# (autor, 2026-10-11: colores "más bonitos", y que no se confundan con los de las naves en las otras figuras, que
# eran casi los mismos) arena, violeta, coral y verde azulado: cálido y frío alternados entre vecinos de la pila,
# validados (validate_palette.js: todo pasa; la arena pide rótulos, los lleva, y una línea blanca separa los
# segmentos); la carga que se suelta, en gris cálido (no es equipo: también va a bordo y compite por la carga útil)
# ANTES: '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#b5b3ad' (la paleta de referencia)


def _tinta(col):
    """Texto negro o blanco según la luminancia del fondo (los segmentos oscuros llevan el número en blanco)."""
    from matplotlib.colors import to_rgb
    r, g, b = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in to_rgb(col)]
    return 'white' if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.18 else '#0b0b0b'


def datos_apilado(masas_kg=None):
    """Sin y con freno, los MISMOS componentes y la clase más liviana (autor, 2026-10-11): el rechazo con los haces de
    ACC también en la cota ideal (en la tabla 5 es el resto), el motor en las dos (con motor no hay batería) y el
    compresor del escape en las dos (14 atm); las fijas juntas. Sin freno: la tasa de tiempo real de cada nave.
    masas_kg: la carga que suelta cada nave ({corto: kg}; por omisión 1 t en todas). Los tanques y el PEM van por
    tonelada liberada, así que se escalan con ella; el compresor de H2, en unidades de catálogo (compresor_t)."""
    masas_kg = masas_kg or {}
    import parametros as p
    import minimax as mm
    import fig_esferoides as fe
    from exp05_condensacion_tiempo_real import q_rech, phi_max
    from exp02_quema_centrada import L_ef
    rq = rutas.leer_datos('req_masa_rechazo')
    r, cand = rq['requisitos'], {c['clave']: c for c in rq['candidatos']}
    q, kw = q_rech(0.0, p.phi_objetivo), cand['ACC_haces']['kWkg'][0]
    # a 1 atm, sin el compresor del escape (comparación, autor 2026-10-12; la base sigue a 14 atm): se retiene phi_max
    # a la salida del condensador; los haces pesan lo mismo por MW y el motor va a la tasa de q(0, phi_max). La carta
    # de diseño de EXP-18 lee de aquí su curva de 1 atm
    phi1 = phi_max(p.T_salida_cond)
    q1 = q_rech(0.0, phi1)
    cu = compresor_unidad()
    bp = r['fijas_t'] - cu['kg'] / 1e3                       # b': tanques, PEM y tanque de agua (sin el compresor)
    # el tanque de agua por tonelada liberada (D, 2026-10-13): en una suelta única va en b'; en el reparto guarda el
    # agua de todo lo entregado, así que EXP-20 y EXP-22 lo sacan de b' y lo cargan a la carga de negocio
    w = r['fijas_kg']['tanque_agua'][0] / 1e3

    def fijas(m_t):
        return bp * m_t + compresor_t(m_t, cu)

    def pila(dm, m_t=1.0):
        return [dm * q / 1e3 / kw / 1e3, _compresor_escape_t(dm)[0], r['conversion_t_por_kgps'][0] * dm,
                fijas(m_t), m_t]

    # con freno, el rechazo más grande que cabe: lo que escala con la tasa es lineal en P (t por MW), así que
    # P_max = (carga útil - fijas - la tonelada) / (t por MW); la tonelada se recupera en calor por tonelada / P_max
    por_MW = sum(pila(1e6 / q)[:3])
    # lo mismo con el sistema ACC entero en vez de los haces (el texto lo da: los haces son la cota optimista)
    por_MW_acc = por_MW - pila(1e6 / q)[0] + 1e3 / cand['ACC_sistema']['kWkg'][0] / 1e3
    calor_t = (1000.0 / L_ef(p.phi_objetivo)) * q

    def pila_1(dm, m_t=1.0):
        return [dm * q1 / 1e3 / kw / 1e3, 0.0, r['conversion_t_por_kgps'][0] * dm, fijas(m_t), m_t]
    por_MW_1 = sum(pila_1(1e6 / q1)[:3])
    calor_t_1 = (1000.0 / L_ef(phi1)) * q1
    naves = []
    for c in fe.ORDEN:
        nom = next(n for n, d in p.PLATAFORMAS.items() if d['corto'] == c)
        masa = masas_kg.get(c, 1000.0)
        m_t = masa / 1000.0
        dm = mm.resolver(nom, masa, p.h_max_precision, p.D_lim_oper)['dm_star']
        carga = p.PLATAFORMAS[nom]['carga_max'] / 1e3
        libre = carga - fijas(m_t) - m_t                                            # lo que queda para lo que escala
        P_max, P_acc = libre / por_MW * 1e6, libre / por_MW_acc * 1e6               # [W]
        Q = calor_t * m_t                                                           # calor de la carga liberada [J]
        naves.append({'corto': c, 'nombre': fe.CORTO[c], 'carga_t': carga, 'masa_t': m_t,
                      'techo_t': suelta_unica_max(carga, bp, 0.0, cu), 'compresor_n': unidades_compresor(m_t, cu),
                      'sin': {'mdot_kgps': dm, 'P_MW': dm * q / 1e6, 't_nat_s': masa / (mm.L_eff * dm),
                              't': pila(dm, m_t)},
                      'con': {'mdot_kgps': p.P_ref_masa / q, 'P_MW': p.P_ref_masa / 1e6,
                              't_rec_min': Q / p.P_ref_masa / 60.0, 't': pila(p.P_ref_masa / q, m_t)},
                      'max': None if P_max <= 0 else {'mdot_kgps': P_max / q, 'P_MW': P_max / 1e6,
                                                     't_rec_min': Q / P_max / 60.0, 't': pila(P_max / q, m_t)},
                      'max_acc': None if P_acc <= 0 else {'P_MW': P_acc / 1e6, 't_rec_min': Q / P_acc / 60.0},
                      'max_una_atm': None if libre <= 0 else {
                          'mdot_kgps': libre / por_MW_1 * 1e6 / q1, 'P_MW': libre / por_MW_1,
                          't_rec_min': calor_t_1 * m_t / (libre / por_MW_1 * 1e6) / 60.0,
                          't': pila_1(libre / por_MW_1 * 1e6 / q1, m_t)}})
    return {'componentes': [n for n, _ in APILADO], 't_por_MW': por_MW, 't_por_MW_acc': por_MW_acc,
            'una_atm': {'phi': phi1, 't_por_MW': por_MW_1, 'calor_t_J': calor_t_1},
            'calor_t_J': calor_t, 'fijas_t': r['fijas_t'], 'fijas_sin_comp_t': bp, 'agua_t': w, 'compresor': cu,
            'naves': naves}


def _pila_h(ax, y, masas, x_max, alto, fs):
    """Barra horizontal apilada: lo fijo y lo chico a la izquierda (autor, 2026-10-11), es decir APILADO al revés.
    Rótulo adentro de cada segmento si tiene ancho para él; devuelve el total."""
    base = 0.0
    for (nombre, col), m in zip(APILADO[::-1], masas[::-1]):
        ax.barh(y, m, left=base, height=alto, color=col, edgecolor='white', lw=0.8)
        if m > 0.07 * x_max:
            ax.text(base + m / 2.0, y, _fmt(m), ha='center', va='center', fontsize=fs, color=_tinta(col))
        base += m
    return base


def figura_apilado(d):
    """Por nave, dos barras horizontales apiladas, una debajo de otra (sin freno arriba, con freno abajo), con la
    carga útil; la de con freno, además, en un zoom (a la escala del panel es una franja)."""
    plt = _plt()
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(len(d['naves']), 1, figsize=(7.4, 9.8))
    for ax, nv in zip(axes, d['naves']):
        x_max = 1.16 * max(sum(nv['sin']['t']), nv['carga_t'])
        tot = _pila_h(ax, 1.0, nv['sin']['t'], x_max, 0.5, 8)
        ax.annotate(f'{_fmt(tot)} t', xy=(tot, 1.0), xytext=(5, 0), textcoords='offset points', va='center',
                    ha='left', fontsize=9, color='#0b0b0b')
        _pila_h(ax, 0.0, nv['con']['t'], x_max, 0.5, 8)          # su total se lee en el zoom
        ax.axvline(nv['carga_t'], color='#52514e', lw=1.0, ls=(0, (4, 3)), zorder=3)   # la carga útil (leyenda)
        zx = ax.inset_axes([0.40, 0.17, 0.56, 0.31])     # sus marcas quedan dentro del panel
        z_max = 1.25 * max(sum(nv['con']['t']), nv['carga_t'])
        tz = _pila_h(zx, 0.0, nv['con']['t'], z_max, 0.6, 7.5)
        zx.annotate(f'{_fmt(tz)} t', xy=(tz, 0.0), xytext=(4, 0), textcoords='offset points', va='center', ha='left',
                    fontsize=8, color='#0b0b0b')
        zx.axvline(nv['carga_t'], color='#52514e', lw=1.0, ls=(0, (4, 3)))
        zx.set_xlim(0.0, z_max)
        zx.set_ylim(-0.6, 0.6)
        zx.set_yticks([])
        zx.tick_params(labelsize=7.5, top=False, right=False)
        zx.grid(False)
        zx.set_title('with hold-down, zoom', fontsize=7.5, pad=2, loc='left')
        ax.set_yticks([1.0, 0.0])
        ax.set_yticklabels(['Without\nhold-down', f"With hold-down\n({nv['con']['P_MW']:.0f} MW)"], fontsize=8.5)
        ax.tick_params(axis='y', which='both', left=False, right=False)
        ax.tick_params(top=False)
        ax.set_ylim(-0.45, 1.45)
        ax.set_xlim(0.0, x_max)
        ax.grid(axis='y', visible=False)
        ax.set_title(f"{nv['nombre']}: payload {_fmt(nv['carga_t'])} t ({nv['sin']['mdot_kgps']:.2f} kg/s without)",
                     fontsize=9.5, loc='left')
    axes[-1].set_xlabel('Mass, lightest class (t)')
    fig.legend(handles=[Patch(color=c, label=n) for n, c in APILADO[::-1]] +
                       [Line2D([], [], color='#52514e', lw=1.0, ls=(0, (4, 3)), label='Payload of the airship')],
               loc='upper center', ncol=3, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(rutas.figura('presupuesto_masa_apilado'), dpi=300, bbox_inches='tight')
    plt.close(fig)


# FIGURA DEL PAPER (autor, 2026-10-11): los colores del diagrama del ciclo (figura 1): el rechazo y el motor en la familia
# durazno del diagrama (dos tonos, separados en la pila por el compresor), las fijas (tanques y compresor de H2) en su
# azul; el compresor del escape y la carga no están en el diagrama: grises. Validados en pares vecinos (separación
# normal y daltónica pasan; los grises son neutros a propósito y el durazno es claro como en el diagrama: rótulos)
ESQUEMA = ('#eda04f', '#5f5b55', '#d9711f', '#5b8bd4', '#d6d2cb')     # en el orden de APILADO
# los colores del diagrama del ciclo (autor, 2026-10-12; draw.io: relleno y borde): celeste, compresores y tanques de H2;
# naranja, conversión y rechazo de calor; verde, regeneración; la carga, gris. El segundo de un mismo color, rayado
# la paleta del autor, variante A (2026-10-12): azul, celeste, verde, mostaza y amarillo, en el orden de APILADO
ESTILO_PAPER = {n: (c, 'white', None) for (n, _), c in zip(APILADO, ('#4169E1', '#1EA2E1', '#3CB371', '#E6BC30',
                                                                    '#FFEE00'))}
_CEL, _NAR, _VER, _GRI = ('#dae8fc', '#6c8ebf'), ('#ffe6cc', '#d79b00'), ('#d5e8d4', '#82b366'), ('#f5f5f5', '#666666')
ESTILO_DIAGRAMA = {'Heat rejection (ACC bundles)': (_NAR[0], _NAR[1], None),
                   'Exhaust compressor (14 atm)': (_CEL[0], _CEL[1], '////'),
                   'Conversion (H$_2$ engine)': (_NAR[0], _NAR[1], '////'),
                   'Fixed boxes (tanks, compressor, PEM)': (_CEL[0], _CEL[1], None),
                   'Cargo released (1 t)': (_GRI[0], _GRI[1], None)}


def _tiempo(s):
    if s < 120.0:
        return f'{s:.0f} s'
    return f'{_tmin(s / 60.0)} min' if s < 7200.0 else f'{s / 3600.0:.1f} h'


def figura_apilado_log(d, estilo=None, salida=None):
    """La figura del paper (autor, 2026-10-11): por nave, dos barras horizontales apiladas en escala log, NHD (sin
    freno) y HDM (freno con el rechazo más grande que cabe: el punto de 1 t de la carta de diseño de EXP-18, que llena
    justo la carga útil); a la izquierda qué se instala y cuánto tarda, a la derecha la masa total; colores del
    diagrama; cada corte entre segmentos cae en la masa acumulada real (solo el borde izquierdo es arbitrario). La
    barra de 10 MW salió (quedan los 10 MW en el texto y en la tabla del freno)."""
    plt = _plt()
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    X0 = 0.85 * 0.7 * min(nv['masa_t'] for nv in d['naves'])    # el eje empieza justo antes de la barra más baja
    col = dict(zip([n for n, _ in APILADO], ESQUEMA))
    # estilo opcional {componente: (relleno, borde, rayado)} (autor, 2026-10-12: los colores del diagrama del ciclo)
    est = estilo or {n: (col[n], 'white', None) for n, _ in APILADO}
    todas_1t = all(abs(nv['masa_t'] - 1.0) < 1e-9 for nv in d['naves'])
    x_max = 2.6 * max(sum(nv['sin']['t']) for nv in d['naves'])     # lugar para el rótulo de la barra más larga
    fig, axes = plt.subplots(len(d['naves']), 1, figsize=(7.4, 5.4), sharex=True)
    rotulos = []                       # (eje, texto, desde, hasta, y): se quitan los que no caben en su segmento
    for ax, nv in zip(axes, d['naves']):
        # NHD: sin freno; HDM: freno con el rechazo más grande que cabe; el tiempo al lado del rótulo y la masa total
        # a la derecha de la barra
        # (autor, 2026-10-12) tercera barra: HDM a 1 atm, sin el compresor del escape, dentro de la misma carga útil
        def rot(nombre, x):
            return nombre if x is None else (f"{nombre}   {x['P_MW']:.{1 if x['P_MW'] >= 1 else 2}f} MW   "
                                             f"{_tiempo(x['t_rec_min'] * 60.0)}")
        filas = ((2.0, nv['sin'], f"NHD   {_tiempo(nv['sin']['t_nat_s'])}"),
                 (1.0, nv['max'], rot('HDM 14 atm', nv['max'])),
                 (0.0, nv['max_una_atm'], rot('HDM 1 atm', nv['max_una_atm'])))
        for y, x, _ in filas:
            if x is None:
                ax.text(X0 * 1.15, y, f"{_fmt(nv['masa_t'])} t does not fit: at most "
                        f"{_fmt(nv['techo_t'])} t with a long recovery" if y > 0.5 else
                        'does not fit either: the ceiling does not depend on the pressure', va='center',
                        ha='left', fontsize=7.5, color='#0b0b0b', zorder=4,
                        bbox=dict(facecolor='white', edgecolor='none', pad=1.5))    # la carga útil no lo cruza
                continue
            # (autor, 2026-10-12) todas las barras empiezan en el mismo punto, el de la suelta más chica (la del NT):
            # en las naves de 1 t el segmento de la carga se alarga hacia la izquierda
            acum, x0 = 0.0, 0.7 * min(n_['masa_t'] for n_ in d['naves'])
            for (nombre, _), m in zip(APILADO[::-1], x['t'][::-1]):         # lo fijo y lo chico a la izquierda
                if m <= 0.0:                                                 # sin compresor del escape (1 atm)
                    continue
                ini = max(acum, x0)
                acum += m
                f_, e_, h_ = est[nombre]
                ax.barh(y, acum - ini, left=ini, height=0.72, color=f_, edgecolor=e_, lw=0.7, hatch=h_)
                if acum / ini > 1.25:          # (autor, 2026-10-13) también la carga: 1 t (0.2 t en el NT)
                    tx = ax.text((ini * acum) ** 0.5, y, f'{m:g}' if nombre.startswith('Cargo') else _fmt(m),
                                 ha='center', va='center', fontsize=7.5,
                                 color=_tinta(f_), zorder=5,
                                 bbox=None if not h_ else dict(facecolor=f_, edgecolor='none', pad=0.6))
                    rotulos.append((ax, tx, ini, acum, y))
            ax.annotate(f'{_fmt(round(acum, 6))} t', xy=(acum, y), xytext=(5, 0), textcoords='offset points',
                        va='center', ha='left', fontsize=8, color='#0b0b0b', zorder=4,
                        bbox=dict(facecolor='white', edgecolor='none', pad=1.5))     # la carga útil lo cruza
        # (autor, 2026-10-12) sin la línea vertical de la carga útil: va en el título del panel
        ax.set_xscale('log')
        ax.set_xlim(X0, x_max)
        ax.set_ylim(-0.8, 2.8)            # margen entre las barras y el recuadro (autor)
        ax.set_yticks([2.0, 1.0, 0.0])
        ax.set_yticklabels([f[2] for f in filas], fontsize=8)
        ax.tick_params(axis='y', which='both', left=False, right=False)
        ax.tick_params(top=False)
        ax.grid(axis='y', visible=False)
        # (autor, 2026-10-12) sin negrita y sin "released" (la suelta del NT va en la leyenda de la figura)
        ax.set_title(f"{nv['nombre']}: payload {_fmt(nv['carga_t'])} t", fontsize=9, loc='left', pad=4)
    axes[-1].xaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1000, 10000]))   # solo los que caen en el eje
    axes[-1].xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    axes[-1].xaxis.set_minor_formatter(NullFormatter())
    axes[-1].set_xlabel('Mass, lightest class (t), logarithmic scale', fontsize=9)
    fig.legend(handles=[Patch(facecolor=est[n][0], edgecolor=est[n][1], hatch=est[n][2], lw=0.7,
                              label=n if todas_1t else n.replace(' (1 t)', ''))
                        for n, _ in APILADO[::-1]],
               loc='lower center', ncol=3, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 0.0))   # abajo
    import plot_style
    plot_style.normalizar_paper(fig)
    fig.tight_layout(rect=(0, 0.075, 1, 1), h_pad=1.3)    # aire entre el título y el panel de arriba
    # (2026-10-12) a 7.5 pt, el número de un segmento solo va si cabe en él
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    for ax_, tx, xa, xb, y_ in rotulos:
        ancho = ax_.transData.transform((xb, y_))[0] - ax_.transData.transform((xa, y_))[0]
        if tx.get_window_extent(ren).width > ancho - 3.0:
            tx.remove()
    if salida:
        fig.savefig(salida, dpi=300, bbox_inches='tight')
    else:
        plot_style.guardar_paper(fig, rutas.figura('presupuesto_masa_apilado_log'))
    plt.close(fig)


MASAS_FIG = {'ZNT': 200.0}                 # [kg] suelta del Zeppelin NT en la figura del paper (autor, 2026-10-12)


def main():
    r = datos()
    naves = datos_naves()
    fr = datos_freno()
    ap = datos_apilado()
    figura_apilado(ap)
    # la figura del paper (autor, 2026-10-12): el Zeppelin NT suelta 0.2 t, que sí cabe aun con el compresor de H2 de
    # catálogo fijo (EXP-20: a lo más 0.24 t); es el 20 % de su carga útil, como 1 t en la ZR-3. Las macros Pila* y los
    # chequeos siguen con 1 t en todas las naves ('apilado'); la figura usa 'apilado_fig'
    ap_fig = datos_apilado(MASAS_FIG)
    figura_apilado_log(ap_fig, estilo=ESTILO_PAPER)
    for nv in ap['naves']:
        print(f"{nv['nombre']:<12} carga {nv['carga_t']:5.1f} t | sin freno: " +
              ' + '.join(f'{x:.1f}' for x in nv['sin']['t']) + f" = {sum(nv['sin']['t']):.0f} t | con freno: " +
              ' + '.join(f'{x:.2f}' for x in nv['con']['t']) + f" = {sum(nv['con']['t']):.1f} t")
    rutas.guardar_datos('fig_presupuesto_masa', dict(r, naves=naves, freno=fr, apilado=ap, apilado_fig=ap_fig,
                                                     masas_fig_kg=MASAS_FIG))
    figura_pct(r)
    figura_log(r)
    figura_naves(naves)
    figura_freno_naves(fr)
    print(f"freno: {fr['P_MW']:.0f} MW, {fr['mdot_kgps']:.4f} kg/s, {fr['t_rec_min']:.1f} min; total "
          f"{fr['total_t'][0]:.1f}-{fr['total_t'][1]:.1f} t; " +
          ', '.join(f"{n}: {t[0]:.2f}-{t[1]:.2f} t" for n, t, _ in fr['cajas']))
    for nv in naves:
        print(f"{nv['nombre']:<12} carga {nv['carga_t']:5.1f} t, {nv['mdot_kgps']:.3f} kg/s, al condensador le queda "
              f"{nv['queda_t']:7.1f} t")
    for c in r['cajas']:
        print(f"{c['nombre']:<24} {c['t'][0]:8.2f} - {c['t'][1]:8.2f} t  "
              f"({100 * c['t'][0] / r['carga_t']:.2f} - {100 * c['t'][1] / r['carga_t']:.2f} % de la carga útil)")


if __name__ == '__main__':
    main()
