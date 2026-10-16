# -*- coding: utf-8 -*-
# [EN] Runs the whole pipeline in order, each script in its own process; outputs go to resultados/.
# =============================================================================
# correr_todo.py — regenera TODOS los resultados del modelo, en orden.
#
# Uso (desde la carpeta modelo/):
#     python correr_todo.py            todo (~6 min)
#     python correr_todo.py --rapido   sin los dos pasos lentos (~3 min)
#
# Cada script corre en su propio proceso y deja:
#     resultados/figuras, resultados/datos, resultados/tablas  (sus salidas)
#     resultados/logs/<script>.txt                            (su consola)
# Al final exportar_numeros.py arma resultados/numeros.tex.
#
# ¿Qué cambió respecto a la última versión guardada?   git diff -- resultados/
# =============================================================================

import os
import sys
import io
import time
import subprocess
import rutas

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ORDEN = [
    ('autores.py',                     'autores del manuscrito (no es parte del modelo)'),
    ('apendice_dlim_flujo.py',         'minimax vs D_lim (4 plataformas)'),
    ('condensacion.py',                'potencia mínima de condensación'),
    ('irradiancia_casco.py',           'irradiancia sobre los paneles del casco (SENAMHI + factor de forma)'),
    ('comparacion_O2_corrected.py',    'fuente de O2 + regeneración solar'),
    ('t2_gain_scheduled.py',           'esquema estándar (1 y 4 entregas)'),
    ('drop_flujo_optimo.py',           'minimax, caso base (parametros.nombre_dirigible)'),
    ('drop_flujo_optimo_perdrop.py',   'minimax, 4 plataformas + multi-drop'),
    ('reactivo_vs_limpia.py',          'control reactivo puro vs quema limpia, mismo L (sec. 4, tabla 3)'),
    ('t3_barrido_reactivo.py',         'barrido reactivo (Tabla 3)'),
    ('apendice_E_sweep.py',            'barrido de carga 1-30 t (Ap. E)'),
    ('apendice_E_flujo_vs_payload.py', 'figura del Ap. E'),
    ('apendice_F_gen.py',              'tablas del Ap. F'),
    ('mapa_requisitos.py',             'mapa de requisitos verificado (Paso 1)'),
    ('exp02_quema_centrada.py',        'EXP-02: quema centrada sin reinyección (fórmulas)'),
    ('exp05_condensacion_tiempo_real.py', 'EXP-05: condensación en tiempo real (fórmulas)'),
    ('exp06_requisitos_minimos.py',    'EXP-06: requisitos mínimos, condensar todo (fórmulas)'),
    ('exp09_condensacion_dinamica.py', 'EXP-09: simulaciones del artículo con dinámica del condensador'),
    ('exp10_escalones_condensador.py', 'EXP-10: liberación por escalones con retardo del condensador'),
    ('exp14_calor_por_temperatura.py', 'EXP-14: calor por nivel de temperatura (banda caliente y fría)'),
    ('exp15_altura.py',                'EXP-15: operación en altura (atmósfera estándar desde el día caluroso)'),
    ('exp11_balance_energia.py',       'EXP-11: cajas negras, balance de energía y regeneración real'),
    ('exp12_helio.py',                 'EXP-12: comparación con helio (cajas negras)'),
    ('val01_componentes.py',           'VAL-01: validación por componente contra datos publicados'),
    ('val02_condensacion.py',          'VAL-02: Ec. 5 contra experimentos H2-aire publicados'),
    ('req_aire_casco.py',              'requisito de admisión de aire al casco (paper 1)'),
    ('tablas_articulo.py',             'tablas del artículo (inglés), desde los resultados'),
    ('diagrama_ciclo.py',              'figura 1 del paper 1 (diagrama del ciclo)'),
    ('figura_operacion.py',            'figura de la operación de la maniobra (paper 1)'),
    ('fig_esferoides.py',              'figura de los cascos como esferoides y k de Lamb (paper 1)'),
    ('balance_masa.py',                'balance de masa por caja negra (paper 1)'),
    ('req_masa_rechazo.py',            'requisito de masa del rechazo y brecha con lo que existe (paper 1)'),
    ('helices_vectoriales.py',         'hélices que la nave ya tiene, fuera de la premisa (discusión del paper 1)'),
    ('exp16_balance_neto.py',          'EXP-16: balance neto por tonelada, materia y energía (con compresión)'),
    ('freno.py',                       'el freno como 4.ª caja negra: reparto, recuperación, energía (paper 1)'),
    ('fig_presupuesto_masa.py',        'presupuesto de masa sin y con freno, las cuatro naves (figura del paper 1)'),
    ('exp18_por_nave.py',              'EXP-18: cascos en 3D y alturas de las cuatro naves (figuras del paper 1), tablas por nave'),
    ('exp20_reparto_sueltas.py',       'EXP-20: repartir la carga en sueltas con el freno (figura 2x2 del paper 1)'),
    ('exp21_balance_energia.py',       'EXP-21: balance de masa y energía por tonelada (lo lee la verificación)'),
    ('exp22_energia_potencia.py',      'EXP-22: energía y potencia por etapa; la entrega de referencia de la sección 5.3'),
    ('derivaciones.py',                'pasos intermedios de las ecuaciones del paper 1'),
    ('ver02_numerica.py',              'VER-02: verificación numérica (solución exacta, convergencia)'),
    ('val03_cotas.py',                 'VAL-03: cotas y sensibilidades de la PIRT (tau, arrastre, gas, incertidumbre)'),
    ('credibilidad.py',                'CRED-01: evaluación de credibilidad (NASA-STD-7009B), sensibilidades'),
    ('vv_resumen.py',                  'V&V: tabla única de verificación y validación del paper 1'),
    ('fuentes.py',                     'VAL-00: registro y auditoría de las fuentes de validación'),
    ('verificacion.py',                'VER-01: verificación interna del subsistema'),
    ('exportar_numeros.py',            'numeros.tex'),
]
LENTOS = {'t3_barrido_reactivo.py', 'mapa_requisitos.py'}
# (2026-10-16) no vienen en la copia pública (ninguno escribe algo que lea otro script): si faltan, se saltan
OPCIONALES = {'autores.py', 'diagrama_ciclo.py', 'apendice_E_sweep.py', 'apendice_E_flujo_vs_payload.py',
              'apendice_F_gen.py'}


def main():
    rapido = '--rapido' in sys.argv
    env = dict(os.environ, PIPELINE='1', MPLBACKEND='Agg', PYTHONIOENCODING='utf-8')
    fallos = []
    t_inicio = time.time()
    print(f'Corriendo el modelo en {rutas.MODELO}\n')
    for script, que in ORDEN:
        if script in OPCIONALES and not os.path.exists(os.path.join(rutas.MODELO, script)):
            print(f'  --     {script:<32}          (no está en esta copia)')
            continue
        if rapido and script in LENTOS:
            print(f'  --     {script:<32}          (omitido: --rapido)')
            continue
        t0 = time.time()
        r = subprocess.run([sys.executable, script], cwd=rutas.MODELO, env=env,
                           capture_output=True, text=True, encoding='utf-8', errors='replace')
        log = r.stdout + (f'\n--- STDERR ---\n{r.stderr}' if r.stderr.strip() else '')
        with open(os.path.join(rutas.LOGS, script.replace('.py', '.txt')), 'w', encoding='utf-8') as f:
            f.write(log)
        estado = 'ok' if r.returncode == 0 else 'ERROR'
        print(f'  {estado:<6} {script:<32} {time.time() - t0:6.1f} s   {que}')
        if r.returncode != 0:
            fallos.append(script)
            print(f'         -> ver resultados/logs/{script.replace(".py", ".txt")}')
            if script == 'exportar_numeros.py':
                break
    print(f'\nTotal: {time.time() - t_inicio:.0f} s')
    if fallos:
        print(f'FALLARON: {", ".join(fallos)}')
        sys.exit(1)
    print('Todo en orden. Para ver qué cambió:  git diff -- resultados/')


if __name__ == '__main__':
    main()
