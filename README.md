# A regenerative hydrogen cycle for post-release buoyancy compensation in heavy-lift airships: real-time condensation limits and hold-down operation

Antony Elthon Espinoza Leon<sup>a</sup>, Jorge Enrique Ortiz Porras<sup>b</sup>, Flor de María Laura Ortiz Pinedo<sup>b</sup>, Rosiand Manrique Alva<sup>b</sup>, Francisco Javier Wong Cabanillas<sup>b</sup>, Alcides Guillermo Joo Aguayo<sup>b</sup>, Gina Coral Tejada Estrada<sup>b</sup>, Ysabel Zevallos Parave<sup>b</sup>, Fidel Tadeo Soria Cuellar<sup>c</sup>

<sup>a</sup> National University of Engineering (UNI), Lima, Peru
<sup>b</sup> National University of San Marcos (UNMSM), Lima, Peru
<sup>c</sup> National University of Education Enrique Guzmán y Valle (UNE), Lima, Peru

Corresponding author: Antony Elthon Espinoza Leon (UNI).

**This repository contains only the model and the verification of this article.** Other studies of the authors are not included.

## Status

Manuscript under submission (Elsevier). The PDFs in this repository are the submitted version and have not been peer reviewed.

## Contents

| Path | Description |
|---|---|
| `manuscript.pdf` | Submitted manuscript |
| `supplementary_material.pdf` | Supplementary material: credibility assessment with NASA-STD-7009B and details of the results |
| `modelo/` | Python model that generates every number, table and figure of the article (`modelo` is Spanish for *model*) |
| `modelo/resultados/` | The outputs used by the article: `numeros.tex`, the generated tables and figures, and the verification log |
| `modelo/datos_fuente/` | Reference data used by the model (excerpts of NIST-JANAF and the NIST Chemistry WebBook) |

## How to run

Tested with Python 3.14 on Windows (numpy 2.4.2, scipy 1.17.0, matplotlib 3.10.8). LaTeX is not needed.

```
pip install -r modelo/requirements.txt
cd modelo
python correr_todo.py            # whole pipeline, about 6 min
python correr_todo.py --rapido   # without the two slow sweeps, about 3 min
```

Each script runs in its own process and writes to `modelo/resultados/` (`datos/`, `tablas/`, `figuras/`, and its console output in `logs/`). The last two steps are:

- `verificacion.py`: about 1100 automated checks (conservation of mass and energy in the simulations, closed forms against simulation, an exact solution and numerical convergence, consistency between scripts that compute the same quantity, every derivation step against an independent computation, and the reference data). It stores a SHA-256 seal of the data it verified.
- `exportar_numeros.py`: writes `resultados/numeros.tex`. It refuses to write if any check fails or if the data changed after verification.

## Where each number comes from

Every number of the model in the article is a LaTeX macro (`\val...`) defined in `modelo/resultados/numeros.tex`; the manuscript and the supplementary material contain no number of the model written by hand. The article uses about 325 of these macros; the others are intermediate values kept for transparency and for the checks.

### Tables

| Article | File in `resultados/tablas/` | Script |
|---|---|---|
| Table 2 | `art_plataformas_rev.tex` | `tablas_articulo.py` |
| Table 3 | `art_cota_ideal_rev2.tex` | `tablas_articulo.py` |
| Table 4 | `art_ambiente.tex` | `tablas_articulo.py` |
| Table 6 | `art_helio_rev2.tex` | `tablas_articulo.py` |
| Table S3 | `art_validacion.tex` | `tablas_articulo.py` |
| Table S4 | `art_sensibilidad.tex` | `credibilidad.py` |
| Table S6 | `art_reactivo.tex` | `tablas_articulo.py` |
| Table S7 | `art_masas.tex` | `req_masa_rechazo.py` |
| Table S8 | `art_brecha_rechazo.tex` | `req_masa_rechazo.py` |
| Table S9 | `art_potencia_escalon.tex` | `tablas_articulo.py` |
| Table S10 | `art_freno_rev.tex` | `freno.py` |

Tables 1 and 5 of the article and Tables S1, S2 and S5 of the supplementary material are written in the manuscript, with their numbers taken from `numeros.tex`.

### Figures

| Article | File in `resultados/figuras/` | Script |
|---|---|---|
| Figure 1 | — | cycle diagram drawn by hand (not code) |
| Figure 2 | `exp18_cascos_3d.png` | `exp18_por_nave.py` |
| Figure 3 | `exp18_alturas_por_nave.png` | `exp18_por_nave.py` |
| Figure 4 | `exp18_diseno_carga_tiempo.png` | `exp18_por_nave.py` |
| Figure 5 | `presupuesto_masa_apilado_log.png` | `fig_presupuesto_masa.py` |
| Figure 6 | `exp22_carga_naves.png` | `exp22_energia_potencia.py` |
| Figure 7 | `exp22_potencia_suelta_1atm.png` | `exp22_energia_potencia.py` |
| Figure 8 | `exp22_regen_nomograma.png` | `exp22_energia_potencia.py` |
| Figure S1 | `val02_condensacion.png` | `val02_condensacion.py` |
| Figure S2 | `credibilidad.png` | `credibilidad.py` |
| Figure S3 | `operacion_maniobra.png` | `figura_operacion.py` |
| Figure S4 | `limpia_vs_reactivo.png` | `reactivo_vs_limpia.py` |
| Figure S5 | `exp20_reparto_naves.png` | `exp20_reparto_sueltas.py` |
| Figure S6 | `exp20_reparto_universal.png` | `exp20_reparto_sueltas.py` |

### Sections of the article and main scripts

| Section | Scripts |
|---|---|
| 3 Modeling | `parametros.py` (all constants), `inercia.py`, `minimax.py`, `condensador.py`, `drop_flujo_optimo_perdrop.py`, `exp05_condensacion_tiempo_real.py` |
| 4 Ideal bound | `minimax.py`, `exp02_quema_centrada.py`, `reactivo_vs_limpia.py`, `t3_barrido_reactivo.py`, `apendice_dlim_flujo.py` |
| 5.1 Without a hold-down device | `exp06_requisitos_minimos.py`, `balance_masa.py`, `req_masa_rechazo.py`, `req_aire_casco.py`, `exp09_condensacion_dinamica.py`, `exp18_por_nave.py` |
| 5.2 With a hold-down device | `freno.py`, `helices_vectoriales.py`, `fig_presupuesto_masa.py`, `exp18_por_nave.py` |
| 5.3 Delivering in releases | `exp20_reparto_sueltas.py`, `exp22_energia_potencia.py` |
| 6 Energy and regeneration | `exp21_balance_energia.py`, `exp22_energia_potencia.py`, `exp11_balance_energia.py`, `exp16_balance_neto.py`, `irradiancia_casco.py`, `comparacion_O2_corrected.py` |
| 7 Comparison with other methods | `exp12_helio.py`, `tablas_articulo.py` |
| 8 Verification, validation and sensitivities | `verificacion.py`, `ver02_numerica.py`, `val01_componentes.py`, `val02_condensacion.py`, `val03_cotas.py`, `credibilidad.py`, `vv_resumen.py`, `fuentes.py` |
| 9 Discussion | `exp14_calor_por_temperatura.py`, `exp15_altura.py`, `val03_cotas.py`, `fig_esferoides.py` |
| Appendices A–C | `derivaciones.py` |

## Sources

`fuentes.py` keeps the register of the published sources used by the model, each checked with a six-step protocol and identified by the SHA-256 fingerprint of its document. The documents themselves are not distributed (copyright); without them, the identity check of a source is reported as "not checked" and does not fail.

## Language

The code comments are in Spanish. This README and the header of each script are in English.

## How to cite

Please cite the article and the Zenodo record: https://doi.org/10.5281/zenodo.23396619

## License

- Code (`modelo/*.py`): MIT License (see `LICENSE`).
- Manuscript and supplementary material (PDFs): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- `modelo/datos_fuente/`: excerpts of NIST Standard Reference Data (NIST-JANAF Thermochemical Tables; NIST Chemistry WebBook), redistributed with attribution for reproducibility; they are not covered by the MIT License.
