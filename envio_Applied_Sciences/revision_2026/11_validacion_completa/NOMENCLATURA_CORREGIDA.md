# Nomenclatura derivada de los grafos implementados

La inspección durante la ejecución encontró una limitación adicional al informe inicial: no hay ocho familias clásicas distintas representadas. Los nombres del código se conservan como identificadores históricos, pero no deben usarse como certificación taxonómica.

En el grafo de cuerpos, Watt tiene sus dos cuerpos ternarios adyacentes; Stephenson los tiene separados por cuerpos binarios. La distinción se documenta en la fuente académica de STU Bratislava y en la publicación de síntesis de Stephenson III. El script `topology_audit.py` reconstruye conexiones a partir de las funciones del simulador; no deduce la topología del nombre del archivo.

| ID histórico | Etiqueta nueva | Segundo bucle | Clasificación del grafo implementado |
|---|---|---|---|
| fourbar | F4 | Ninguno | Cuatro barras |
| slidercrank | SC | Guía prismática | Slider-crank, geometría rígida corregida |
| wattI | CD-G | C sobre AB, C-D-G | Stephenson, ground ternario, no “Watt I” |
| wattII | RC-G | C sobre rocker B-O4, C-D-G | Watt con ground ternario |
| stephensonI | A-O4 | A-D-O4 | Watt con ground y crank ternarios; articulaciones coincidentes |
| stephensonII | B-O2 | B-D-O2 | Stephenson con articulaciones coincidentes |
| stephensonIII | C-O4 | C sobre AB, C-D-O4 | Stephenson con ground coincidente |
| gearedFiveBar | G5 | Accionamientos relacionados | Cinco barras con relación discreta compatible con un ciclo |

C-O4 es una especialización de CD-G fijando G=O4. A-O4 traza mediante el subbucle O2-A-D-O4, mientras que el bucle con B restringe validez pero no determina adicionalmente el trazador; no representa un espacio de curvas independiente de todo cuatro barras. Las comparaciones son por modelo, cajas de parámetros y localización del trazador, no por todas las inversiones o capacidades intrínsecas de una familia.

Consecuencia editorial: usar “eight kinematic models” en título, objetivo, tablas y conclusiones; mostrar etiquetas nuevas y equivalencia con IDs heredados. Explicar las especializaciones y parámetros redundantes. No conservar el reclamo de comparar las cinco inversiones clásicas Watt/Stephenson. Las campañas mantienen sus IDs y ecuaciones para no destruir trazabilidad; el cambio de nomenclatura no cambia resultados numéricos.

Fuentes:

- STU Bratislava, material académico sobre síntesis de mecanismos: https://www.sjf.stuba.sk/buxus/docs/pracoviska/ATC/PDF/TM_III/12_L_Synthesis.pdf
- Plecnik y McCarthy, Kinematic synthesis of Stephenson III six-bar function generators: https://www.sciencedirect.com/science/article/pii/S0094114X15002384
