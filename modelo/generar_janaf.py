# -*- coding: utf-8 -*-
# [EN] Builds janaf_datos.py from the NIST-JANAF tables in datos_fuente/janaf/.
# =============================================================================
# generar_janaf.py — escribe janaf_datos.py copiando las filas de las tablas NIST-JANAF
# guardadas en datos_fuente/janaf/<id>.txt (texto plano de
# https://janaf.nist.gov/tables/<id>.txt, consulta 2026-09-24). No está en correr_todo.py:
# se corre solo si cambian las tablas. Las funciones están en janaf.py.
# =============================================================================

import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(AQUI, 'datos_fuente', 'janaf')
TABLAS = (('N2', 'N-023'), ('O2', 'O-029'), ('Ar', 'Ar-001'), ('H2O', 'H-064'), ('CO2', 'C-095'), ('H2', 'H-050'))
M = {'N2': '28.0134e-3', 'O2': '31.9988e-3', 'Ar': '39.948e-3', 'H2O': '18.01528e-3', 'CO2': '44.0095e-3',
     'H2': '2.01588e-3'}                                   # [kg/mol] IUPAC
T_MIN, T_MAX = 200.0, 3000.0


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    out = ['# -*- coding: utf-8 -*-',
           '# GENERADO por generar_janaf.py desde datos_fuente/janaf/*.txt. NO EDITAR A MANO.',
           '# Tablas NIST-JANAF (Chase, M. W., NIST-JANAF Thermochemical Tables, 4th ed.,',
           '# J. Phys. Chem. Ref. Data Monograph 9, 1998), consulta 2026-09-24.',
           '# especie: (id JANAF, M [kg/mol] (IUPAC), dfH°(298.15) [kJ/mol], filas), con filas',
           '#   (T [K], Cp [J/mol K], S° [J/mol K] (p° = 0.1 MPa), H - H(298.15) [kJ/mol]).',
           '', 'DATOS = {']
    for esp, tid in TABLAS:
        filas, dfh = [], None
        for linea in open(os.path.join(D, f'{tid}.txt'), encoding='utf-8'):
            c = linea.strip().split('\t')
            try:
                T = float(c[0])
            except ValueError:
                continue
            if c[0] == '298.15':
                dfh = float(c[5])
            if T_MIN <= T <= T_MAX:
                filas.append((T, float(c[1]), float(c[2]), float(c[4])))
        out.append(f"    '{esp}': ('{tid}', {M[esp]}, {dfh!r}, (")
        out += [f'        ({f[0]!r}, {f[1]!r}, {f[2]!r}, {f[3]!r}),' for f in filas]
        out.append('    )),')
        print(f'{esp:4s} {tid}: {len(filas)} filas {filas[0][0]:g}-{filas[-1][0]:g} K, dfH {dfh} kJ/mol')
    out.append('}')
    with open(os.path.join(AQUI, 'janaf_datos.py'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(out) + '\n')
    print('janaf_datos.py escrito')


if __name__ == '__main__':
    main()
