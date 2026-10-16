# -*- coding: utf-8 -*-
# [EN] Closed form of the burn centred on the release (Section 4, Appendix A).
# =============================================================================
# exp02_quema_centrada.py — EXP-02 (análisis, SIN simular): quema constante
# centrada en la liberación, sin reinyección, de reposo a reposo.
#
# Modelo (ver bitacora/EXP02_quema_centrada.tex para la derivación):
#   (1+k)·M·ÿ = −L·g·B(t) + m·g·H(t − T/2)          B(t) = ṁ·t (masa quemada)
#   - parte en reposo en y0 con flotación(y0) = M·g        (condición del autor)
#   - k: masa añadida vertical (inercia, no peso)
#   - L(φ) = ρ_aire/ρ_H2 − 1 + 9φ: φ = fracción del agua retenida al instante
#   - sin reinyección (ṁ >= 0); total quemado = m/L (vuelve al equilibrio de y0)
#   - arrastre vertical ½ρC_D A v²: corrección de PRIMER ORDEN (válida si
#     kappa = c·H/((1+k)M) <~ 0.1; si es mayor, solo orientativa; 'nan' si no aplica)
# Resultado cerrado: con la liberación en T/2 la nave termina en REPOSO (v(T) = 0)
# y solo baja: H = m·g·T² / (24·(1+k)·M)
#   => T = √(24 (1+k) M H / (m g)),  ṁ = m / (L T),  v_max = m g T / (8 (1+k) M),
#      a_max = m g / (2 (1+k) M).
# Se evalúan fórmulas (instantáneo): no es una simulación.
# =============================================================================

import io
import sys
import csv
import math
import parametros as p
import rutas
import minimax as mm

H_LISTA = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0]      # bajada total permitida [m]
RHO_RATIO = mm.rho_aire / mm.rho_H2                      # ρ_aire/ρ_H2


def L_ef(phi):
    """Fuerza neta por kg de H2 quemado (en kg): pérdida de flotación + O2 retenido - H2 que sale."""
    return RHO_RATIO - 1.0 + (1.0 + p.ratio_O2_H2) * phi


def quema_centrada(M, m, H, k, phi, c=0.0):
    """c = ½·ρ·C_D·A [kg/m]: arrastre vertical como corrección de PRIMER ORDEN.
    Con arrastre la bajada a igual T es H0·(1 − 0.9·c·H0/Me) (H0 = bajada ideal), así que
    para una bajada permitida H se resuelve H0 = [1 − √(1 − 3.6 c H/Me)] / (1.8 c/Me) y
    T = √(24 Me H0/(m g)). kappa = c·H/Me mide si la corrección es chica (<~0.1)."""
    Me = (1.0 + k) * M
    kappa = c * H / Me
    if c > 0:
        disc = 1.0 - 3.6 * kappa
        H0 = (1.0 - math.sqrt(disc)) / (1.8 * c / Me) if disc > 0 else float('nan')
    else:
        H0 = H
    T = math.sqrt(24.0 * Me * H0 / (m * p.g)) if H0 == H0 else float('nan')
    L = L_ef(phi)
    return {'T_s': T, 'dm_kgps': m / (L * T), 'v_max_ms': m * p.g * T / (8.0 * Me),
            'a_max_ms2': m * p.g / (2.0 * Me), 'H2_total_kg': m / L, 'L': L,
            'kappa': kappa, 'alpha_reposo': 0.5 + 0.075 * kappa}


def cargas_eval(nombre):
    cap = min(p.PLATAFORMAS[nombre]['carga_max'] / 1000.0, 30.0)
    return sorted({1.0, cap} | ({20.0} if cap >= 20 else set()))


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    k0, phi0 = p.k_masa_anadida, p.phi_agua_retenida
    filas = []
    for n, d in p.PLATAFORMAS.items():
        kn = d['k']                                        # masa añadida de esta nave (su esbeltez)
        variantes = {'ideal':             (0.0, 1.0, False),   # sin masa añadida (el artículo hasta EXP-13)
                     'con_masa_anadida':  (kn, 1.0, False),
                     'con_arrastre':      (kn, 1.0, True),    # + arrastre (1er orden) = supuestos acordados
                     'escape_sale':       (kn, 0.0, False)}   # agua NO retenida al instante
        for m_t in cargas_eval(n):
            m = m_t * 1000.0
            M = p.rho_aire_std * d['V_casco']          # aire del casco: la inercia es (1 + k)·M (inercia.py)
            for H in H_LISTA:
                fila = {'plataforma': n, 'corto': d['corto'], 'carga_t': m_t, 'H_m': H, 'M_total_kg': M}
                c = 0.5 * p.rho_aire_std * p.Cd * d['Area']
                for nom_v, (k, phi, con_arrastre) in variantes.items():
                    r = quema_centrada(M, m, H, k, phi, c if con_arrastre else 0.0)
                    for kk, vv in r.items():
                        fila[f'{nom_v}_{kk}'] = vv
                filas.append(fila)
    with open(rutas.dato('exp02_quema_centrada', 'csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        for r in filas:
            w.writerow({kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in r.items()})

    def buscar(n, m_t, H):
        return next(f for f in filas if f['plataforma'] == n and f['carga_t'] == m_t and f['H_m'] == H)

    # Comparación con el punto del artículo (caso base, 1 t, excursión 1 m): minimax ṁ* vs quema centrada
    dl = rutas.leer_datos('dlim_flujo')
    pf = buscar(p.nombre_dirigible, 1.0, 1.0)
    fw = buscar('Flying Whales LCA60T', 20.0, 30.0)
    fw5 = buscar('Flying Whales LCA60T', 20.0, 5.0)
    rutas.guardar_datos('exp02_quema_centrada', {
        'k_masa_anadida': k0, 'L_phi1': L_ef(1.0), 'L_phi0': L_ef(0.0),
        'factor_masa_anadida': 1.0 / math.sqrt(1.0 + k0), 'factor_escape': L_ef(1.0) / L_ef(0.0),
        'base_1t_1m': {'minimax_articulo_kgps': dl['plataformas'][p.nombre_dirigible]['dm_star_kgps'],
                     **{kk: pf[kk] for kk in pf if kk.endswith(('dm_kgps', 'T_s'))}},
        'FW_20t_30m': {kk: fw[kk] for kk in fw if kk.endswith(('dm_kgps', 'T_s', 'v_max_ms', 'kappa', 'alpha_reposo'))},
        'FW_20t_5m': {kk: fw5[kk] for kk in fw5 if kk.endswith(('dm_kgps',))},
    })

    print(f'L(φ=1) = {L_ef(1):.2f}   L(φ=0) = {L_ef(0):.2f}   k = {k0}  '
          f'(masa añadida: flujo x{1/math.sqrt(1+k0):.3f}; escape: flujo x{L_ef(1)/L_ef(0):.2f})\n')
    print(f'{"plataforma":<22}{"carga":>6}{"H m":>6}{"T s":>7}{"ṁ ideal":>9}{"ṁ +M_a":>9}{"+arrastre":>10}{"kappa":>7}{"ṁ escape":>10}{"v_max":>7}')
    for f in filas:
        if f['carga_t'] == cargas_eval(f['plataforma'])[-1] or f['carga_t'] == 20.0:
            print(f'{f["plataforma"]:<22}{f["carga_t"]:>6g}{f["H_m"]:>6g}{f["con_masa_anadida_T_s"]:>7.1f}'
                  f'{f["ideal_dm_kgps"]:>9.1f}{f["con_masa_anadida_dm_kgps"]:>9.1f}'
                  f'{f["con_arrastre_dm_kgps"]:>10.1f}{f["con_arrastre_kappa"]:>7.3f}'
                  f'{f["escape_sale_dm_kgps"]:>10.1f}{f["con_masa_anadida_v_max_ms"]:>7.2f}')
    print(f'\n{p.nombre_dirigible}, 1 t, 1 m: minimax del artículo {dl["plataformas"][p.nombre_dirigible]["dm_star_kgps"]:.2f} kg/s'
          f' vs quema centrada ideal {pf["ideal_dm_kgps"]:.2f} kg/s (con masa añadida {pf["con_masa_anadida_dm_kgps"]:.2f})')
    graficar(filas)


def graficar(filas):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import plot_style  # noqa
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.8), sharey=True)
    for n, d in p.PLATAFORMAS.items():
        cap = cargas_eval(n)[-1]
        fs = [f for f in filas if f['plataforma'] == n and f['carga_t'] == cap]
        a1.plot([f['H_m'] for f in fs], [f['con_arrastre_dm_kgps'] for f in fs], color=d['color'],
                lw=2.2, marker='o', ms=4, label=f'{n} ({cap:g} t)')
    a1.set_title('(a) Each airship at its payload capacity (added mass + drag)', fontsize=10)
    fs = [f for f in filas if f['plataforma'] == 'Flying Whales LCA60T' and f['carga_t'] == 20.0]
    col = p.PLATAFORMAS['Flying Whales LCA60T']['color']
    for clave, ls, lab in (('ideal', ':', 'As in the article (no added mass, water kept)'),
                           ('con_masa_anadida', '-.', f'+ added mass (k = {p.k_masa_anadida:g})'),
                           ('con_arrastre', '-', '+ added mass + drag (first order)'),
                           ('escape_sale', '--', '+ added mass, exhaust leaves (water not kept)')):
        a2.plot([f['H_m'] for f in fs], [f[f'{clave}_dm_kgps'] for f in fs], color=col, ls=ls, lw=2.2, label=lab)
    a2.set_title('(b) Flying Whales LCA60T, 20 t: model assumptions', fontsize=10)
    for ax in (a1, a2):
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xticks(H_LISTA); ax.set_xticklabels([f'{h:g}' for h in H_LISTA])
        ax.set_xlabel(r'Allowed total descent  $H$  [m]')
        ax.grid(True, which='both', alpha=0.2)
        ax.legend(fontsize=8, loc='upper right')
    a1.set_ylabel(r'Constant H$_2$ burn rate  $\dot{m}$  [kg/s]')
    fig.tight_layout()
    fig.savefig(rutas.figura('exp02_flujo_vs_bajada'), dpi=200, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
