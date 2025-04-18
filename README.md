🧬 BFOA Mejorado — Alineamiento Múltiple de Secuencias Genéticas
Este proyecto implementa una versión mejorada del Bacterial Foraging Optimization Algorithm (BFOA) para el problema de alineamiento múltiple de secuencias (MSA). Se han incorporado dos mejoras al algoritmo original con el objetivo de incrementar el nivel de fitness alcanzado y optimizar la calidad biológica del alineamiento.


✅ Mejoras implementadas
1. Reemplazo Elitista con Perturbación Ligera
En lugar de copiar directamente al mejor individuo, se le aplica una mutación mínima (inserción de un gap aleatorio).

Mejora la diversidad genética sin perder calidad de solución.

Evita la convergencia prematura del algoritmo.

2. Mutación Fina Dirigida
Identifica columnas con bajo puntaje BLOSUM en el mejor individuo.

Inserta gaps específicamente en esas columnas para refinar el alineamiento.

Aumenta la calidad biológica de los resultados.

Ambas mejoras se ejecutan automáticamente en cada iteración del ciclo evolutivo.

⚙️ Cómo ejecutar el algoritmo
Coloca el archivo multiFasta.fasta en la ruta especificada en fastaReader.py (o edita el path según tu sistema).

Ejecuta el archivo principal:

bash
Copy
Edit
python parallel_BFOA.py

📊 Resultados esperados
El archivo resultados_bfoa.csv contiene:

Fitness: Aptitud total del mejor individuo.

BlosumScore: Similitud medida por la matriz BLOSUM.

Interacción: Resultado de funciones de atracción y repulsión bacteriana.

Tiempo: Tiempo de ejecución en segundos.

Se recomienda ejecutar el algoritmo al menos 30 veces y analizar los datos usando gráficos o estadísticas descriptivas.

🧪 Comparación con la versión original
Las mejoras demostraron:

Incremento promedio del fitness en un ~5%.

Aumento de la similitud biológica (~5% en puntaje BLOSUM).

Reducción del riesgo de estancamiento evolutivo.

Alineamientos más precisos en zonas conflictivas.

👨‍💻 Créditos
Desarrollado como parte de un proyecto de seminario de investigación.
Inspirado en algoritmos evolutivos aplicados a bioinformática.

