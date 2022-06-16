# TRABAJO FIN DE MÁSTER -- Gonzalo Soria
Repositorio conteniendo los **scripts** y algunos de los **ficheros generados** en el análisis llevado a cabo durante el TFM realizado por Gonzalo Soria.

## PASO 2: Preprocesado y normalización de los datos
En este segundo paso, se lleva a cabo el preprocesado de los datos por medio del software de Bollito. Se crea un entorno de conda que permite su ejecución.
Dada la gran diferencia en cuanto a número de genes y moléculas detectadas entre los distintos tipos de elementos a analizar (AM vs. tumor vs. PIC), el filtrado se realiza en tres pasos por separado y luego se unifican los resultados.

***Input:***
* Los ficheros de UMIs (AB1017\[1-8\]\_noTcells.txt) y metadata (AB1017\[1-8\]\_no_Tcells_metadata.txt) generados en el **PASO 1**.
* Los ficheros samples.tsv: como vamos a hacerlo por separado, se generan los ficheros *samples_AM.tsv, samples_tumor.tsv y samples_PIC.tsv*. Este fichero únicamente contiene dos columnaes indicando (1) el nombre de la muestra y (2) la condición (AM, tumor ó PIC, respectivamente).
* Los ficheros units.tsv: como vamos a hacerlo por separado, se generan los ficheros *units_AM.tsv, units_tumor.tsv y units_PIC.tsv*. Este fichero contiene 3 columnas indicando (1) el nombre de la muestra, tal y como se indica en su respectivo fichero samples_\*.tsv; (2) el path al fichero donde se encuentra la matriz de expresión correspondiente a cada muestra y (3) el path al fichero donde se encuentran los metadatos correspondientes a cada muestra.
* Los ficheros config.yaml: como vamos a hacerlo por separado, se generan los ficheros *units_AM.tsv, units_tumor.tsv y units_PIC.tsv*. Este fichero contiene toda la información necesaria para ejecutar el software de bollito. Por ejemplo, se especifica para cada uno de los 3 ficheros, qué ficheros de units y samples va a emplear como input, qué directorio va a generar como output, qué tipo de datos vamos a analizar, etc. Además, incluye también los parámetros específicos (*tabla 1 de la memoria*) empleados para filtrar células y genes poco representativos, y que varía entre los 3 tipos de células. La mediana se extrae con el código que se muestra al final del script 1.2, del objeto seurat generado previamente al procesamiento.


Hay que tener en cuenta que, cada vez que queramos ejecutar bollito, hay que renombrar el config_\*.yaml que queramos emplear a config.yaml, para que se ejecute correctamente.

Además, tal y como se indica en el script 2.1, el análisis se ejecuta en tres pasos: 
* Primero, el filtrado de los ficheros de los 3 tipos celulares.
* Segundo, se unifican los ficheros .rds generados en cada directorio out\_\*/ desde R con Seurat (v.3.2.0, que es la misma que emplea bollito) tal y como se especifica en el script 2.2.
* Tercero, se continua con bollito para ejecutar la normalización. 

***Output:*** 
* Directorios con los resultados del preprocesado para los 3 tipos celulares (directorios out\_\*/ especificados en cada config.yaml)
* Objeto seurat **seurat_find-clusters.rds** ya procesado que se empleará en el análisis de aquí en adelante, en la ubicación: "./out_AMs/seurat/merged/3_clustering/". 

**FIGURAS:** 
* *Figura 5B* generada automáticamente por bollito (al especificarle la condición a emplear para el clustering en el config_AM.yaml): "./out_AMs/seurat/merged/3_clustering/5_lisi_score_by_group"
* *Figura 7* (script 2.3)
