# Proyecto DBSCAN Cómputo Paralelo y en la Nube

## Descripción del Proyecto
Implementación paralela del algoritmo DBSCAN (Density-Based Spatial Clustering of Applications with Noise) utilizando OpenMP para la detección de outliers.

## Estructura del Proyecto

### Directorios Principales
- `input/`: Archivos CSV con datos de entrada (20K a 200K puntos)  
- `src/`: Código fuente de las implementaciones  
- `output/`: Resultados de clasificación y métricas de rendimiento
  - `clasificacion/`: Archivos de salida con puntos etiquetados  
  - `rendimiento/`: Métricas de tiempo y speedup  

### Implementaciones
- `serial.cpp`: Versión secuencial de referencia  
- `omp1.cpp`: Primera estrategia paralela (enfoque unificado)  
- `omp2.cpp`: Segunda estrategia paralela (división espacial con grid)  

### Notebooks de Análisis
- `DBSCAN_noise.ipynb`:
  - Generación de datasets sintéticos para pruebas  
  - Visualización comparativa entre sklearn y nuestras implementaciones  
  - Validación de resultados mediante gráficos de dispersión  

- `DBSCAN_analyze.ipynb`:
  - Análisis de rendimiento y escalabilidad  
  - Cálculo de speedups y eficiencia  
  - Gráficas comparativas entre versiones serial y paralelas  

## Ejecución

### Compilación
```bash
/opt/homebrew/bin/g++-15 -O3 -std=c++17 -fopenmp src/serial.cpp -o src/serial
/opt/homebrew/bin/g++-15 -O3 -std=c++17 -fopenmp src/omp1.cpp -o src/omp1
/opt/homebrew/bin/g++-15 -O3 -std=c++17 -fopenmp src/omp2.cpp -o src/omp2
````

### Ejecución Individual

```bash
./src/omp2 --threads 8 --in input/20000_data.csv --eps 0.03 --minpts 10 --out output/20000_results.csv
```

### Experimentos Completos

```bash
./run_experiments.sh
```

## Parámetros de Ejecución

* `--threads`: Número de hilos OpenMP
* `--in`: Archivo de entrada CSV
* `--eps`: Radio épsilon para vecindarios
* `--minpts`: Mínimo de puntos para ser core
* `--out`: Archivo de salida con clasificación

## Resultados

El proyecto genera:

* Archivos de clasificación: Puntos etiquetados como core/border/noise
* Métricas de rendimiento: Tiempos de ejecución y speedups
* Visualizaciones: Comparativas con implementación de referencia