# MASTER_THESIS
Scripts, input &amp; output files y todos los ficheros supplementarios empleados en el análisis de la MSc Thesis se han subido a este repositorio de github.

Analysis starts with 8 UMI tab delimited files (AB10171-8.txt) and some metadata files containing per-cell annotations.

## PASO 1
En el primer paso del análisis, se ejecuta el script 1.1 desde el software de R (v. 4.2.0) para eliminar las células T de mi análisis, y generar una visión global de los datos iniciales proporcionados. Además, se divide el fichero de metadatos proporcionado inicialmente generando un fichero de metadatos por cada uno de los ficheros de datos iniciales. Esta división es necesaria para ejecutar Bollito en el **PASO 2**
El script 1.2 permite obtener por Seurat (v.4.0) un UMAP plot de los datos antes de su procesado.

## PASO 2
En este segundo paso, se lleva a cabo el preprocesado de los datos por medio del software de Bollito. Se crea un entorno de conda que permite su ejecución.
Dada la gran diferencia en cuanto a número de genes y moléculas detectadas entre los distintos tipos de elementos a analizar (AM vs. tumor vs. PIC), el filtrado se realiza en tres pasos por separado y luego se unifican los resultados.

Para ello, se necesitan:
* Los ficheros de UMIs (AB1017\[1-8\]\_noTcells.txt) y metadata (AB1017\[1-8\]\_no_Tcells_metadata.txt) generados en el **PASO 1**.
* Los ficheros samples.tsv: como vamos a hacerlo por separado, se generan los ficheros *samples_AM.tsv, samples_tumor.tsv y samples_PIC.tsv*. Este fichero únicamente contiene dos columnaes indicando (1) el nombre de la muestra y (2) la condición (AM, tumor ó PIC, respectivamente).
* Los ficheros units.tsv: como vamos a hacerlo por separado, se generan los ficheros *units_AM.tsv, units_tumor.tsv y units_PIC.tsv*. Este fichero contiene 3 columnas indicando (1) el nombre de la muestra, tal y como se indica en su respectivo fichero samples_\*.tsv; (2) el path al fichero donde se encuentra la matriz de expresión correspondiente a cada muestra y (3) el path al fichero donde se encuentran los metadatos correspondientes a cada muestra.
* Los ficheros config.yaml: como vamos a hacerlo por separado, se generan los ficheros *units_AM.tsv, units_tumor.tsv y units_PIC.tsv*. Este fichero contiene toda la información necesaria para ejecutar el software de bollito. Por ejemplo, se especifica para cada uno de los 3 ficheros, qué ficheros de units y samples va a emplear como input, que directorio va a generar como output, qué tipo de datos vamos a analizar, etc. Además, incluye también los parámetros específicos (*tabla 1 de la memoria*) empleados para filtrar células y genes poco representativos, y que varía entre los 3 tipos de células.


Hay que tener en cuenta que, cada vez que queramos ejecutar bollito, hay que renombrar el config_\*.yaml que queramos emplear a config.yaml, para que se ejecute correctamente.

Además, tal y como se indica en el script 2.1, el análisis se ejecuta en tres pasos: 
* Primero, el filtrado de los ficheros de los 3 tipos celulares.
* Segundo, se unifican los ficheros .rds generados en cada directorio out\_\*/ desde R con Seurat (v.3.2.0, que es la misma que emplea bollito) tal y como se especifica en el script 2.2.
* Tercero, se continua con bollito para ejecutar la normalización. 

***Output:*** Objeto seurat ya procesado que se empleará en el análisis de aquí en adelante, en la ubicación: "./out_AMs/seurat/merged/3_clustering/seurat_find-clusters.rds". Se muestra en el script 2.3 cómo obtener la figura 7 como overview de los datos ya procesados.

## PASO 3
Una vez normalizados los datos, se procede al análisis transcriptómico *per se*.
Se emplea Seurat (V.4.0.0) de aquí en adelante para mí análisis.
Siguiendo el script 3.1, se generan (entre otros) los ficheros:
* *DEGs\_EMTPICs\_EMT.tsv* conteniendo los 76 DEGs tumorales iniciales
* *DEGs\_EMTPICs\_EMT\_final.tsv* conteniendo los 29 genes tumorales finales
* *DEGs_PICs_AM.tsv* conteniendo los 749 genes asociados a AM iniciales
* *DEGs_PICs_AM_final.tsv* conteniendo los 362 genes tumorales finales
* los ficheros Degs
