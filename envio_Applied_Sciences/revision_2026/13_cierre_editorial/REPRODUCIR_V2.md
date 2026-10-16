# Reproducir la revisión corregida

Extraer el suplemento editorial conservando la estructura `FEF-Graph/`. Los módulos QD están en `FEF-Graph/qd_experiments/`, los datos históricos en `FEF-Graph/repo_codigo_nuevo/datos/` y la revisión en `FEF-Graph/envio_Applied_Sciences/revision_2026/`. El archivo raíz `MANIFEST_SHA256.json` identifica los bytes de cada entrada. No se distribuye el entorno virtual de la máquina autora.

La rama pública `revision/applsci-4642776` de `romusan/fef-graph-capability-benchmark` conserva esa disposición para motor y revisión; los datos históricos están en `datos/` y los lectores admiten ambas ubicaciones. Su `PUBLICATION_MANIFEST_SHA256.json` identifica los archivos de la ampliación y `verify_publication.py` comprueba la evidencia científica. Las evaluaciones originales, la carta y los respaldos editoriales no forman parte de ese depósito público. Para regenerar allí los resultados, usar los pasos numéricos siguientes y compilar directamente el `.tex` canónico; el constructor editorial de cartas requiere el respaldo privado y solo se usa en el paquete editorial completo.

Desde `revision_2026/`, crear un entorno Python 3.12 e instalar `12_campanas_corregidas/requirements.lock.txt`. El entorno numérico usado y las versiones por corrida se conservan. Se usa Numba de doble precisión sin fastmath; no hay dependencia de GPU. PDF requiere MiKTeX/LaTeX y Poppler, además de pypdf/pdfplumber/Pillow en el entorno de documentos.

## Regenerar análisis sin nuevas búsquedas

Ejecutar con el Python del entorno numérico, en este orden:

```text
python 13_cierre_editorial/summarize_v2.py
python 13_cierre_editorial/extended_analysis.py
python 13_cierre_editorial/validate_resolution_v2.py
python 13_cierre_editorial/render_tables_figures.py
python 13_cierre_editorial/finish_diagnostics.py
python 13_cierre_editorial/build_documents.py
```

Las tablas rechazan campañas incompletas y requieren semillas 1–30 por celda. `summarize_v2.py --allow-partial` escribe exclusivamente en `partial/`; no produce tablas finales. La validación de resolución no ejecuta optimización ni selecciona vectores nuevos. Los CSV de tolerancia y resolución y las coordenadas de objetivos quedan en `13_cierre_editorial/results/`.

El constructor de documentos usa la copia anterior conservada en `13_cierre_editorial/antes_correcciones_codex.zip` para datos de autores, referencias y comentarios. **Al regenerar documentos se debe volver a compilar y revisar los PDF**, actualizar el registro de revisión visual y cerrar la matriz; los hashes del control anterior ya no prueban la nueva salida. `pdf_qa.py` registra compilación y renderizado, pero no se autoaprueba visualmente. Sus rutas de herramientas son las de Windows utilizadas en esta entrega y deben adaptarse en otra instalación. `verify_release.py` comprueba los PDF locales mediante los hashes del cierre registrado.

## Repetir búsquedas

Los cuatro archivos de motor congelados no se deben editar para reutilizar resultados: la caché compara el hash y rechaza versiones incompatibles. Las ejecuciones existentes se reutilizan sin repetir cálculo; para una réplica nueva utilizar **otra copia extraída**, conservar los resultados suministrados como referencia y mover las carpetas de resultados de esa copia a una ubicación de respaldo antes de ejecutar:

```text
python 12_campanas_corregidas/run_campaigns.py factorial
python 12_campanas_corregidas/run_campaigns.py easy
python 12_campanas_corregidas/run_campaigns.py independent
python 12_campanas_corregidas/run_campaigns.py new_targets
python 12_campanas_corregidas/run_campaigns.py singularity
python 12_campanas_corregidas/run_campaigns.py repeatability
python 12_campanas_corregidas/run_engineering_campaign.py
python 12_campanas_corregidas/source/self_recovery_v2.py
python 12_campanas_corregidas/run_qd_campaigns.py
```

La recuperación escribe en `11_validacion_completa/self_recovery_v2/`; respaldar también esa carpeta al repetir. Los controladores usan como máximo tres trabajadores. Pilotos 871 y pruebas de determinismo no se incorporan a las tasas de evaluación. No mezclar resultados históricos, corregidos y nuevos de una réplica.

El objetivo NACA0012 usa el archivo de coordenadas suministrado, sin necesidad de descargar de nuevo la web. Conserva su procedencia y hash; no contiene medidas experimentales de un mecanismo. Para repetir análisis históricos, `07_statistics/reanalisis_estadistico.py` usa los datos adyacentes incluidos y `11_validacion_completa/validate_archived.py` reevalúa vectores sin corregirlos silenciosamente.

## Indicadores y denominadores

`any_strict` se refiere a cualquier candidato evaluado con RMS≤1 mm y MAX≤2 mm dentro del presupuesto. Movimiento y calidad en las tablas principales corresponden al mejor candidato estricto retenido por el criterio de coste, verificado posteriormente; no se busca retrospectivamente entre todos los candidatos un ganador físico. El análisis de resolución usa el candidato mínimo-RMS a 20k, que puede ser distinto. QD conserva candidatos estrictos por modelo/celda de referencia y cuenta éxito por semilla, nunca por número de élites como si fueran réplicas independientes. Los checkpoints son observaciones dependientes de la misma búsqueda.
