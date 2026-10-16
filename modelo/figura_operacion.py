# -*- coding: utf-8 -*-
# [EN] Operation of the maneuver, base case (Figure S3).
# =============================================================================
# figura_operacion.py — figura de la operación de la maniobra para el paper 1:
# el caso base (parametros.nombre_dirigible, 1 t, umbral de 1 m) con condensación instantánea,
# con las fases numeradas (el pie de figura las nombra):
#   1 la unidad de rechazo de calor arranca (n_tau_arranque x tau antes de la quema)
#   2 quema previa (precondicionamiento) a la tasa instalada
#   3 quema reactiva con el PD después de soltar
#   4 ajuste fino hasta el nuevo equilibrio
# =============================================================================

import numpy as np
import parametros as p
import rutas
import drop_flujo_optimo_perdrop as sim


def main():
    sim.usar_nave(p.nombre_dirigible, 1000.0)
    r = sim.correr(n_entregas=1, toneladas_por_drop=1.0, graficar=False, verbose=False, devolver_serie=True)
    s = r['serie']
    t = np.asarray(s['t'])
    t_rel = s['t_drops'][0]
    t_ant = r['t_anticipo']
    quema = -np.gradient(np.asarray(s['m_H2']), t)
    t_arr = t_rel - t_ant - p.n_tau_arranque * p.tau_condensador
    tras = (t > t_rel) & (np.abs(quema) > 0.05 * np.abs(quema).max())
    t_fin = float(t[tras][-1]) if tras.any() else t_rel          # fin de la quema reactiva
    t_max = t_rel + 90.0
    rutas.guardar_datos('figura_operacion', {'t_release_s': t_rel, 't_anticip_s': t_ant, 't_arranque_s': t_arr,
                                             't_fin_reactiva_s': t_fin, 'dm_opt': r['dm_opt']})
    import matplotlib.pyplot as plt
    import plot_style  # noqa
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.6), sharex=True)
    fases = [(t_arr, t_rel - t_ant, '#E69F00', '1'), (t_rel - t_ant, t_rel, '#0072B2', '2'),
             (t_rel, t_fin, '#009E73', '3'), (t_fin, t_max, '#999999', '4')]
    for ax in axes:
        for a, b, col, _ in fases:
            ax.axvspan(a, b, color=col, alpha=0.14, lw=0)
        ax.axvline(t_rel, color='black', lw=1, ls='--')
        ax.grid(True, alpha=0.2)
    axes[0].plot(t, np.asarray(s['h']) - p.h0, color='black', lw=1.8)
    axes[0].set_ylabel('Altitude change [m]')
    y0, y1 = axes[0].get_ylim()
    for a, b, col, txt in fases:
        axes[0].text(0.5 * (a + b), y1, txt, ha='center', va='bottom', fontsize=11, fontweight='bold', color=col)
    axes[0].text(t_rel, 0.5 * y0, ' release', fontsize=9, va='center')
    axes[1].plot(t, quema, color='#D55E00', lw=1.8)
    axes[1].axhline(0.0, color='gray', lw=0.8)
    axes[1].set_ylabel('H$_2$ burn rate [kg/s]')
    axes[1].set_xlabel('Time [s]')
    axes[1].set_xlim(t_arr - 5.0, t_max)
    fig.tight_layout()
    plot_style.guardar_paper(fig, rutas.figura('operacion_maniobra'), dpi=150)
    plt.close(fig)
    print(f'liberación {t_rel:.1f} s, anticipación {t_ant:.2f} s, arranque {t_arr:.1f} s, '
          f'fin de la quema reactiva {t_fin:.1f} s')


if __name__ == '__main__':
    main()
