# -*- coding: utf-8 -*-
# [EN] Summary table of verification and validation (supplementary material).
# =============================================================================
# vv_resumen.py — V6 del plan de V&V: TABLA ÚNICA de verificación y validación del
# paper 1, al estilo de ASME V&V 20 (comparación = modelo S contra referencia D,
# diferencia E y la incertidumbre de la referencia u_D, con su veredicto). Junta lo
# que VAL-01, VAL-02, VAL-03 y VER-02 dejan en sus JSON; no calcula nada nuevo.
# Salidas: resultados/tablas/vv_resumen.tex y resultados/vv_resumen.txt
# =============================================================================

import os
import rutas
import parametros as p

AQUI = os.path.dirname(os.path.abspath(__file__))


def filas():
    v1 = rutas.leer_datos('val01_componentes')
    v2 = rutas.leer_datos('val02_condensacion')
    v3 = rutas.leer_datos('val03_cotas')
    r2 = rutas.leer_datos('ver02_numerica')
    ra = rutas.leer_datos('req_aire_casco')
    ma = v1['masa_anadida']['medida_sobre_analitica']
    # las mismas definiciones de las macros de VAL-02 (exportar_numeros.py): contra su apéndice A recalculado,
    # experimentos 1-5 (su tabla tiene los errores registrados en fuentes.ZAK22), y lo recogido / la cota
    dz = [abs(e['esperada_g'] / e['esperada_zak_recalc_g'] - 1.0) for e in v2['experimentos'][:5]]
    km = [e['medido_sobre_modelo'] for e in v2['experimentos']]
    t = v3['tau']
    ex = {a['config']: a for a in r2['exacta']}
    un = max(max(x['u_rel'].values()) for fs in r2['convergencia'].values() for x in fs
             if x['config'].startswith('adoptada'))
    cd = max(abs(f['d_excursion_rel']) for fs in v3['cd'].values() for f in fs)
    de = v3['termica']['deriva']
    dq = (de['5K_minimax']['dm_kg'][1], de['20K_minimax']['dm_kg'][1])
    dp = (de['5K_pasos_10MW_1t']['dm_kg'][1] / 1e3, de['10K_pasos_10MW_1t']['dm_kg'][1] / 1e3)
    return [
        # (tipo, bloque PIRT, comparación, modelo S, referencia D, diferencia E, u_D o criterio, veredicto)
        ('validación', 'P1-P2', f"{v1['n_filas']} constantes contra su fuente primaria", 'modelo', 'ISA, CODATA, IUPAC, NIST, IF97, IEA',
         f"{v1['n_dentro']}/{v1['n_filas']} dentro del criterio", 'criterio por dato', 'validado'),
        # sus curvas y su agua esperada salen del MISMO balance (Tetens): consistencia, no validación (hallazgo M2)
        ('consistencia', 'P3', 'cota de condensación vs su curva teórica recalculada', 'Ec. 5 (IF97)',
         'Zakhvatkin 2022, fig. 4', f"{100 * v2['dif_max_curvas']:.2f} puntos", 'Tetens vs IF97', 'coherente'),
        ('consistencia', 'P3', 'agua esperada con flujos medidos', 'Ec. 5', 'Zakhvatkin 2022, apéndice A',
         f"{100 * min(dz):.1f}-{100 * max(dz):.1f} %", 'transcripción + Tetens', 'coherente'),
        ('validación', 'P3', 'agua recogida / cota', 'cota superior', 'Zakhvatkin 2022, tabla 1',
         f"{100 * min(km):.0f}-{100 * max(km):.0f} % de la cota", 'dispersión entre ensayos', 'cota confirmada; propagada'),
        ('validación', 'P4', 'masa vertical añadida', 'analítica (k)', 'Yamasaki & Goto 2003, tabla 4',
         f"medido/analítico {ma['restringido']:.2f}-{ma['vuelo_libre']:.2f}", 'dos ensayos', 'rango; propagado'),
        ('cota con datos', 'P8', 'retardo del condensador tau', f"{t['tau_supuesto_s']:g} s",
         'c_Al/k (FlyZero, JANAF); Ghatos 2023',
         f"{min(t['tau_FlyZero_s']):.1f}-{max(t['tau_FlyZero_s']):.1f} s; lab. {t['tau_Ghatos_forzada_s'][0]:.0f}-"
         f"{t['tau_Ghatos_forzada_s'][1]:.0f} s", 'orden de magnitud', 'respaldado'),
        ('sensibilidad', 'P5', 'arrastre vertical Cd 0-1', f'Cd {p.Cd:g}', 'JPL (0.8), Robbins & Martone (~0.5)',
         f"<= {100 * cd:.1f} % en la excursión", '--', 'irrelevante'),
        ('sensibilidad', 'P6', 'gas isotérmico (deriva térmica)', 'isotérmico', 'Churchill & Chu + NIST',
         f"{dq[0]:.0f}-{dq[1]:.0f} kg en la quema minimax; {dp[0]:.1f}-{dp[1]:.1f} t en 1 t por escalones",
         'dT_eq 5-20 K (quema), 5-10 K (escalones)', 'limitación declarada'),
        ('requisito', 'P7', 'admisión de aire al casco', 'instantánea (supuesto)', 'sin datos publicados',
         f"{min(x['aire_1m_m3s'] for x in ra['filas'] if x['corto'] != 'ZNT'):.0f}-"
         f"{max(x['aire_1m_m3s'] for x in ra['filas'] if x['corto'] != 'ZNT'):.0f} m3/s cota ideal; "
         f"{ra['filas'][0]['aire_escalon_m3s']:.2f} m3/s en escalones", 'simulación = ṁ/rho_H2',
         'requisito declarado'),
        ('verificación', 'P13', 'código vs solución exacta', 'simulación', 'v = v_t tanh(t/t_v)',
         f"{100 * ex['paso /2, rtol 1e-8']['err_rel_z']:.3f} %", 'tolerancia estricta', 'verificado'),
        ('verificación', 'P13', 'incertidumbre numérica (convergencia)', 'rtol 1e-6', 'rtol 1e-10, paso /4',
         f"<= {100 * un:.2f} %", '--', 'verificado'),
    ]


def main():
    F = filas()
    L = [r'\begin{tabular}{l l p{4.2cm} p{3.4cm} p{3.4cm} l}', r'\toprule',
         r'Tipo & PIRT & Comparación & Referencia & Diferencia & Veredicto \\', r'\midrule']
    for tipo, pb, comp, S, D, E, u, ver in F:
        esc = lambda s: s.replace('%', r'\%').replace('&', r'\&').replace('_', r'\_')
        L.append(f'{tipo} & {pb} & {esc(comp)} & {esc(D)} & {esc(E)} & {esc(ver)} \\\\')
    L += [r'\bottomrule', r'\end{tabular}']
    with open(rutas.tabla('vv_resumen'), 'w', encoding='utf-8') as fh:
        fh.write('% GENERADO por modelo/vv_resumen.py. NO EDITAR A MANO.\n' + '\n'.join(L) + '\n')
    T = ['RESUMEN DE VERIFICACION Y VALIDACION DEL PAPER 1 (generado por modelo/vv_resumen.py)',
         'Estilo ASME V&V 20: comparacion, referencia D, diferencia E, incertidumbre/criterio, veredicto.', '']
    for tipo, pb, comp, S, D, E, u, ver in F:
        T.append(f'[{tipo} | {pb}] {comp}\n    modelo: {S} | referencia: {D}\n    diferencia: {E} | {u} | {ver}')
    texto = '\n'.join(T) + '\n'
    with open(os.path.join(AQUI, 'resultados', 'vv_resumen.txt'), 'w', encoding='utf-8') as fh:
        fh.write(texto)
    print(texto)


if __name__ == '__main__':
    main()
