# Protocolo de las campañas corregidas

Fijado el 16 de octubre de 2026, después de comprobar la infraestructura con semillas 871 y 877 y antes de ejecutar las semillas de evaluación 1–30. Es una especificación prospectiva de esta ampliación, no una preinscripción del artículo original.

## Motor y validez

Motor histórico congelado en `source/historical_engine.py`. Los kernels de `fast_geometry.py` reproducen las trayectorias y métricas históricas en 800 pruebas aleatorias, excepto cuando se habilitan las reparaciones explícitas. No se usa fastmath. La versión corregida desplaza también el pivote y la manivela del slider cuando cambia x0, y limita su biela a 180/260 mm. En geared five-bar se admiten relaciones −3, −2, −1, 1, 2, 3 para que ambos accionamientos completen un número entero no nulo de vueltas durante un ciclo de entrada. Se proyecta la variable categórica al entero admisible más cercano y se archiva el vector realmente simulado. No se interpreta como optimización diferenciable de una relación continua.

Cada vector seleccionado se revalida con 1920 intervalos y extremo 2π incluido. Se comprueba cierre de estado completo, residuales de eslabones y constancia de signo de áreas orientadas de todas las díadas. La transmisión se mide con ángulo agudo para detectar ambas alineaciones; para el slider se usa el complemento de la oblicuidad de la biela respecto a la guía. El umbral de 30 grados es un filtro de calidad añadido, no una definición universal de factibilidad. Los controles muestreados no son una certificación analítica entre muestras.

## Presupuesto e indicadores

La unidad de coste es una invocación del simulador, válida o inválida, incluidas las diferencias finitas de LS. Se registra el mejor vector encontrado también durante refinamiento. Al agotar el presupuesto se interrumpe la búsqueda conservando lo efectivamente evaluado, sin completar un refinamiento con coste oculto.

Primario: existencia de al menos un candidato con RMS ≤1 mm y MAX ≤2 mm encontrado dentro del presupuesto, acompañado de su verificación de movimiento. Secundarios: RMS del mejor candidato por RMS, resultado del mejor candidato por coste compuesto y éxito más transmisión ≥30 grados. Siempre se conserva el denominador completo de 30 semillas. Los éxitos puramente geométricos permanecen separados de los físicamente válidos.

Factorial: 8 familias ×2 objetivos (medium, medium-hard) ×2 perfiles ×30 semillas = 960 búsquedas. Cada búsqueda tiene 40 000 simulaciones y checkpoints acumulados a 10 000, 20 000 y 40 000, correspondientes a 0.5, 1 y 2 veces B=20 000. Son curvas acumuladas de una misma búsqueda, no 2880 ejecuciones reiniciadas ni observaciones independientes entre presupuestos. La elección B se fija aquí y no depende de qué familia ganó el piloto. Se conserva la configuración del motor medium-hard y se reinicia al terminar naturalmente hasta alcanzar el límite.

Optimizador independiente: implementación oficial pycma, mismos objetivos, simulador, cotas y coste compuesto, semillas 1–30, perfil compact, ambas curvas medium y medium-hard; 20 000 simulaciones por caso. Sin plantillas de FEF. Es una verificación independiente de política, no una prueba de imposibilidad de las familias que fallen.

Objetivos nuevos: Bernoulli con cruce, superelipse de exponente 4 y curva de Fourier con concavidad asimétrica. Fórmulas/procedencia definidas en `campaign.target_points` antes de evaluar; 8 familias ×30 semillas, compact y 20 000 simulaciones. Bernoulli es una referencia matemática publicada; ninguna trayectoria se presenta como medición industrial. Todas se normalizan a ancho 90 mm. Los objetivos originales conservan la discretización histórica en su comparación a N=120.

Singularidad: control con/sin penalización para medium, compact, 8 familias ×30 semillas, 20 000 simulaciones; el control sin penalización elimina también el residual correspondiente de LS. Las filas con penalización se reutilizan solo desde el checkpoint compatible del factorial, sin añadir réplicas ficticias.

## Reproducibilidad y análisis

Semillas de evaluación 1–30. Pilotos 871 y pruebas de determinismo 877 no entran en las tasas. Tres familias se repitieron tres veces con semilla 877 y dieron resultados numéricos idénticos. Se guarda SHA-256 de cuatro fuentes del motor, configuración efectiva por variable y versiones; resultados atómicos y reanudación rechazan fuentes incompatibles. Tres procesos como máximo en cuatro núcleos lógicos.

Las tasas llevan Wilson 95 %. Comparaciones de perfiles se bloquean por semilla dentro del protocolo y reportan diferencias e intervalos; no se agrupan como equivalentes las campañas históricas con configuraciones distintas. Se distinguen análisis descriptivos de contrastes y se declara la familia de multiplicidad para cada bloque. No se equipara no significación con equivalencia ni fracaso de búsqueda con incapacidad cinemática.

QD tendrá especificación separada antes de sus semillas formales, con presupuesto total exacto y todas las componentes experimentales controladas. Los datos históricos se mantienen íntegros y se etiquetan como protocolo anterior.
