# Protocolo de QD corregido

Especificado el 16 de octubre de 2026 antes de las semillas formales 1–30. El piloto usa solo semilla 871. Motor y reparación geométrica: `PROTOCOLO_V2.md`.

Cada brazo usa exactamente 40 000 simulaciones totales sobre medium/compact y ocho familias. Se guarda el vector de cada élite estricto retenido por familia y celda de referencia, incluidos los del baseline; se comprueban a posteriori movimiento y transmisión con 1920 intervalos. Se separan precisión y movimiento; la calidad se refiere a los élites retenidos, no a todos los vectores transitorios visitados.

Se mantienen descriptores, límites congelados del piloto histórico, mutación gaussiana de amplitud 0.12 y surrogate de ridge local ponderado de 30 vecinos del código original. No se ajustan contra las nuevas semillas. Todas las variantes usan el motor corregido, y DE está desactivado. El archivo entre familias es compartido en todas las variantes del factorial: por tanto no es una componente que varíe en esa comparación.

Factorial de componentes: MAP-Elites presente/ausente × surrogate presente/ausente, 30 semillas por combinación, rejilla nativa 24×24. Cuando MAP-Elites está ausente se conserva un conjunto elitista de hasta 576 vectores por RMS, con selección uniforme de padres; no se mantiene un mapa de nichos. Cuando está presente se elige uniformemente entre nichos ocupados. Surrogate activo guía la mitad de las mutaciones tras reunir 120 observaciones válidas por familia, con los parámetros históricos. En ambos casos hay 4000 evaluaciones de inicialización, inmigración 5 %, fase de búsqueda de 24 000 simulaciones y fase de refinamiento de 16 000 (2000 por familia). El refinador es común y ninguna llamada de diferencias finitas queda fuera del presupuesto. Se conserva lo evaluado al agotar cada cuota, sin coste oculto.

Baseline adicional: motor FEF sin DE, 120 muestras, 16 plantillas, vecindarios 6×4 y top-4, refinamiento Powell/LS con configuración histórica de QD. Para imponer un total exacto se asignan 5000 simulaciones a cada una de las ocho familias. Se reinicia si el procedimiento termina antes. Esta es una versión explícitamente acotada del baseline, no una ejecución idéntica al baseline histórico de presupuesto natural.

Sensibilidad: MAP-Elites+surrogate con 12, 18, 24 y 36 celdas por lado, 30 semillas por resolución. El brazo 24 se reutiliza del factorial sin contarlo otra vez. Cobertura y QD-score se evalúan siempre en la misma malla de referencia 24×24; no comparar porcentajes de ocupación de denominadores distintos como si midieran igual diversidad. Las llamadas de refinamiento y las inválidas cuentan en todas las variantes.

Total formal: (4 combinaciones + baseline +3 rejillas adicionales) ×30 =240 ejecuciones. Contrastes de componentes condicionados a este diseño, sin atribuir diferencias del baseline a una sola componente. Reportar éxito por semilla con Wilson 95 %, diferencias por bloque de semilla, incertidumbre y ajuste de Holm sobre los contrastes definidos para cada bloque. No interpretar no significación como equivalencia ni forzar una conclusión favorable a QD.

Archivar resultados, configuración, hashes de fuentes originales y nuevas, vectores, tiempo y verificaciones. No mezclar las 15 semillas históricas con las nuevas para fabricar n=45: los protocolos difieren.
