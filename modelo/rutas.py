# [EN] Paths of the outputs in resultados/.
# =============================================================================
# rutas.py — dónde guarda cada script sus resultados (pipeline, Paso 0)
#
#   resultados/figuras/<nombre>.png   figuras (nombre fijo: se sobreescriben;
#                                     git guarda el historial)
#   resultados/datos/<nombre>.json    números que calcula cada script (los lee
#                                     exportar_numeros.py para armar numeros.tex)
#   resultados/tablas/<nombre>.tex    tablas LaTeX generadas
#   resultados/logs/<script>.txt      salida de consola de cada corrida
#
# Regla: ningún script guarda nada fuera de resultados/.
# =============================================================================

import os
import json

MODELO     = os.path.dirname(os.path.abspath(__file__))
RESULTADOS = os.path.join(MODELO, 'resultados')
FIGURAS    = os.path.join(RESULTADOS, 'figuras')
DATOS      = os.path.join(RESULTADOS, 'datos')
TABLAS     = os.path.join(RESULTADOS, 'tablas')
LOGS       = os.path.join(RESULTADOS, 'logs')

for _d in (FIGURAS, DATOS, TABLAS, LOGS):
    os.makedirs(_d, exist_ok=True)


def figura(nombre):
    """Ruta fija de una figura: resultados/figuras/<nombre>.png"""
    return os.path.join(FIGURAS, f'{nombre}.png')


def tabla(nombre):
    """Ruta fija de una tabla LaTeX: resultados/tablas/<nombre>.tex"""
    return os.path.join(TABLAS, f'{nombre}.tex')


def dato(nombre, ext='json'):
    """Ruta fija de un archivo de datos: resultados/datos/<nombre>.<ext>"""
    return os.path.join(DATOS, f'{nombre}.{ext}')


def _nativo(x):
    """Convierte tipos de numpy a tipos nativos para poder escribir JSON."""
    if isinstance(x, dict):
        return {str(k): _nativo(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_nativo(v) for v in x]
    if hasattr(x, 'item'):          # numpy escalar
        return x.item()
    return x


def guardar_datos(nombre, datos):
    """Guarda los números de un script en resultados/datos/<nombre>.json"""
    with open(dato(nombre), 'w', encoding='utf-8') as f:
        json.dump(_nativo(datos), f, indent=2, ensure_ascii=False)
        f.write('\n')


def sello_datos(excluir=('verificacion.json',)):
    """Huella (SHA-256) de todos los archivos de resultados/datos/, en orden, salvo los excluidos. verificacion.py
    la guarda con sus chequeos y exportar_numeros.py la recalcula: si no coincide, los datos cambiaron después de
    verificar (o la verificación se cayó antes de escribir) y no se exporta (M-6, auditoría 3)."""
    import hashlib
    h = hashlib.sha256()
    for nom in sorted(os.listdir(DATOS)):
        if nom in excluir:
            continue
        h.update(nom.encode('utf-8'))
        with open(os.path.join(DATOS, nom), 'rb') as fh:
            h.update(fh.read())
    return h.hexdigest()


def leer_datos(nombre):
    with open(dato(nombre), encoding='utf-8') as f:
        return json.load(f)


def mostrar():
    """plt.show() solo al correr un script a mano; correr_todo.py lo desactiva."""
    if os.environ.get('PIPELINE') != '1':
        import matplotlib.pyplot as plt
        plt.show()
