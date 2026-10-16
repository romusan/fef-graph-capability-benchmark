# Índice de evidencia para revisores

El suplemento se extrae conservando la carpeta raíz `FEF-Graph/`. Las rutas siguientes parten de `FEF-Graph/envio_Applied_Sciences/revision_2026/`.

| Archivo o carpeta | Qué permite comprobar |
|---|---|
| `11_validacion_completa/topology_audit.json` | Grafos e identidad de los modelos |
| `11_validacion_completa/kernel_verification.json` | Equivalencia del núcleo geométrico |
| `11_validacion_completa/archived_validation.json` | Revalidación de diseños históricos |
| `11_validacion_completa/self_recovery_v2/` | Generadores, objetivos y búsquedas de recuperación propia |
| `13_cierre_editorial/results/analysis_v2.json` | Conteos, incertidumbre y contrastes entre perfiles |
| `13_cierre_editorial/results/extended_analysis_v2.json` | CMA-ES, ablación, componentes QD, rejillas y lotes |
| `13_cierre_editorial/results/per_seed_v2.csv` | Indicadores por semilla y checkpoint |
| `13_cierre_editorial/results/resolution_v2.json` | Vector fijo, resolución y comprobación física |
| `13_cierre_editorial/results/resolution_aggregate.csv` | Conteos de sensibilidad a la resolución |
| `13_cierre_editorial/results/tolerance_aggregate.csv` | Las 16 combinaciones de tolerancia |
| `07_statistics/reanalisis_estadistico.json` | Estadística histórica reparada |
| `12_campanas_corregidas/results/` | 2.640 búsquedas formales por modelo |
| `12_campanas_corregidas/engineering_results/` | 240 búsquedas del contorno de referencia NACA0012 |
| `12_campanas_corregidas/qd_results/` | 240 corridas QD formales más cinco pilotos separados por semilla |
| `15_auditoria_reenvio/model_diagnostics.json` | 240 vectores RC-G/B-O2, errores recalculados, forma y validación; prueba entre lotes |
| `15_auditoria_reenvio/example_*.csv` | Coordenadas de ejemplos cercanos a cada mediana |
| `15_auditoria_reenvio/envelope_diagnostics.json` | 960 selecciones factoriales y envolventes muestreadas |
| `13_cierre_editorial/results/revised_line_locations.json` | Páginas y líneas reales del manuscrito limpio |
| `15_auditoria_reenvio/highlight_report.json` | Base documental y método de resaltado |
| `13_cierre_editorial/REPRODUCIR_V2.md` | Instrucciones para regenerar análisis |

El manifiesto SHA-256 del archivo ZIP identifica los bytes exactos de todas las entradas. El suplemento debe adjuntarse al reenvío o publicarse en un depósito accesible: la presencia local de estos archivos no equivale a que ya estén disponibles para los revisores. No se declara un DOI hasta que el depósito exista.
