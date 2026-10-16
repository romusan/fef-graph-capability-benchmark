# Extensión con contorno de referencia de ingeniería

16 de octubre de 2026, fijado antes de ejecutar sus semillas. Responde a la solicitud de añadir un objetivo procedente de ingeniería además de las curvas sintéticas. No cambia las curvas ni los resultados del factorial en ejecución.

Se usa el contorno alterado NACA0012 de la Turbulence Modeling Resource enlazada oficialmente por NASA. La página especifica su fórmula y cierre de borde de salida; el archivo de coordenadas se conserva íntegro y con SHA-256. Es una geometría de referencia aerodinámica, no una trayectoria medida ni evidencia de que un mecanismo se haya usado industrialmente.

Fuente y coordenadas:

- https://www.nasa.gov/nasa-turbulence-modeling-resource/
- https://tmbwg.github.io/turbmodels/naca0012_val.html
- https://tmbwg.github.io/turbmodels/NACA0012_grids/n0012points_superbig_clust_fix.dat

Preprocesamiento: descartar z=0, preservar el orden del contorno cerrado, remuestrear por longitud de arco, centrar y normalizar a 90 mm de ancho. No suavizar ni ajustar mediante ningún mecanismo. El adaptador de datos queda versionado con hash propio; los cuatro archivos congelados del motor no cambian. Cada salida incluye hash de datos, adaptador y núcleo.

Diseño: ocho modelos ×30 semillas (1–30), compact, motor híbrido, exactamente 20 000 llamadas al simulador por corrida, mismos indicadores y comprobación física que los demás objetivos. Son 240 corridas adicionales en `engineering_results`, separadas de las 2640 corridas ya especificadas. Total de búsquedas a nivel de modelo para la entrega: 2880, más 60 repeticiones, 8 recuperaciones y 240 ejecuciones QD. Los checkpoints no aumentan el n independiente.

El estudio pasa de seis a siete contornos. La referencia de ingeniería amplía el origen geométrico, pero no convierte los datos en mediciones industriales ni justifica afirmaciones de fabricación, precisión de máquina o comportamiento aerodinámico. No se usan las mediciones de fuerzas de la página de validación como si fueran coordenadas de trayectoria.
