# MASTER_THESIS
Repositorio conteniendo los **scripts** y algunos de los **ficheros generados** en el análisis llevado a cabo durante el TFM realizado por Gonzalo Soria.

El análisis
Analysis starts with 8 UMI tab delimited files (AB10171-8.txt) and some metadata files containing per-cell annotations.

## PASO 1
En el primer paso del análisis, se ejecuta el script 1.1 desde el software de R (v. 4.2.0) para eliminar las células T de mi análisis, y generar una visión global de los datos iniciales proporcionados. Además, se divide el fichero de metadatos proporcionado inicialmente generando un fichero de metadatos por cada uno de los ficheros de datos iniciales. Esta división es necesaria para ejecutar Bollito en el **PASO 2**

***Input:***
* Ficheros de cuentas: *AB1017\[1-8\].txt)*
* Ficheros de metadatos: *final_metadata.tsv*
* Ficheros de metacell (PIC-seq annotation): *PICs_metacell.txt ; singlets_metacell.txt*

***Output:***
* *AB1017\[1-8\]\_noTcells.txt*
* *AB1017\[1-8\]\_no_Tcells_metadata.txt*

**FIGURAS:** 
* *Figura 5A* (script 1.2) 
* *Figura 14* (script 1.1)
* *Figura Suplementaria 1* (script 1.1)

## PASO 2
En este segundo paso, se lleva a cabo el preprocesado de los datos por medio del software de Bollito. Se crea un entorno de conda que permite su ejecución.
Dada la gran diferencia en cuanto a número de genes y moléculas detectadas entre los distintos tipos de elementos a analizar (AM vs. tumor vs. PIC), el filtrado se realiza en tres pasos por separado y luego se unifican los resultados.

***Input:***
* Los ficheros de UMIs (AB1017\[1-8\]\_noTcells.txt) y metadata (AB1017\[1-8\]\_no_Tcells_metadata.txt) generados en el **PASO 1**.
* Los ficheros samples.tsv: como vamos a hacerlo por separado, se generan los ficheros *samples_AM.tsv, samples_tumor.tsv y samples_PIC.tsv*. Este fichero únicamente contiene dos columnaes indicando (1) el nombre de la muestra y (2) la condición (AM, tumor ó PIC, respectivamente).
* Los ficheros units.tsv: como vamos a hacerlo por separado, se generan los ficheros *units_AM.tsv, units_tumor.tsv y units_PIC.tsv*. Este fichero contiene 3 columnas indicando (1) el nombre de la muestra, tal y como se indica en su respectivo fichero samples_\*.tsv; (2) el path al fichero donde se encuentra la matriz de expresión correspondiente a cada muestra y (3) el path al fichero donde se encuentran los metadatos correspondientes a cada muestra.
* Los ficheros config.yaml: como vamos a hacerlo por separado, se generan los ficheros *units_AM.tsv, units_tumor.tsv y units_PIC.tsv*. Este fichero contiene toda la información necesaria para ejecutar el software de bollito. Por ejemplo, se especifica para cada uno de los 3 ficheros, qué ficheros de units y samples va a emplear como input, que directorio va a generar como output, qué tipo de datos vamos a analizar, etc. Además, incluye también los parámetros específicos (*tabla 1 de la memoria*) empleados para filtrar células y genes poco representativos, y que varía entre los 3 tipos de células.


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


## PASO 3
Una vez normalizados los datos, se procede al análisis transcriptómico *per se*.
Se emplea Seurat (V.4.0.0) de aquí en adelante para mí análisis.
En el script 3.1, se extraen los *DEGs* empleados en mi análisis.
En el script 3.2 se realiza el análisis de vías de enriquecimiento sobre estos *DEGs*.

***Input:***
* seurat data object: seurat_find-clusters.rds"
* Firmas de D. Cook *et al.*, 2022: *EMP_conserved_signature.tsv* y *Malignant_cell_signature.tsv*

***Output:*** Se generan, entre otros, los siguientes ficheros:
* *DEGs\_EMTPICs\_EMT.tsv* conteniendo los 76 DEGs tumorales iniciales (script 3.1)
* *DEGs\_EMTPICs\_EMT\_final.tsv* conteniendo los 29 genes tumorales finales (script 3.1)
* *DEGs\_PICs\_AM.tsv* conteniendo los 749 genes asociados a AM iniciales (script 3.1)
* *DEGs\_PICs\_AM\_final.tsv* conteniendo los 362 genes tumorales finales (script 3.1)
* El fichero *combined_res.tsv* y *combined_DEGs_names_res.tsv* con los 388 DEGs empleados en el resto del análisis (script 3.1)
* Los ficheros *DEGs\_\*\_noMith.tsv* para el análsis de vías enriquecidas (script 3.2)


**FIGURAS:**
* *Figura 6* (script 3.1)
* *Figura 8A-D* (script 3.2)
* *Figura 9* (script 3.2)

## PASO 4
### NicheNetR Analysis: scripts 4.1 y 4.2
En el script 4.1 se especifica la metodología a seguir para realizar el análisis empleando *NicheNetR* (v.1.0.0).
En el script 4.2 se realiza un filtrado de los ligandos obtenidos en el análisis durante 4.1, que se filtran de acuerdo a su expresión en datos de *AM* generados en **Maria Casanova *et al.,* 2021**. Luego, el fichero aquí generado se carga de nuevo en el script 4.1 para continuar con el análisis.

***Input***:
* seurat data object: seurat_find-clusters.rds"
* Fichero *combined_res.tsv*
* Fichero de datos de expresión de genes en *AMs* de Maria Casanova *et al.,* 2021: "./TRMs_cancer-Maria-data/SUPPLEMENTARY_TABLE_3-ALL_DATA.xlsx"
* (Fichero *best_upstream_ligands_filtered.tsv* generado en el script 4.2 y que se requiere para continuar el script 4.1)

***Output***:
* (Fichero *best_upstream_ligands.tsv* generado en el script 4.2 y que se requiere para continuar el script 4.1)
* Fichero *prioritized_tbl_oi_filtered_top.tsv* conteniendo la información de las 66 interacciones de interés predichas para el nicho de *EMT*.
* Fichero *prioritized_tbl_oi-EMT_low.tsv* conteniendo la información de las 95 interacciones predichas para el nicho pobre en *EMT*.
* Fichero *prioritized_tbl_oi-Singlets.tsv* conteniendo la información de las 92 interacciones predichas para los singletes en lugar de dobletes.

**FIGURAS:**
* *Figura 10A-D* (script 4.1)
* *Figura 11B* (script 4.1)

### CellphoneDB Analysis:
En el script 4.3a se generan los ficheros necesarios para ejecutar cellphoneDB desde R.
En el script 4.3b, se ejecuta cellphone desde el terminal de linux

***Input 4.3a***: 
* seurat data object: seurat_find-clusters.rds"
* Ficheros de genes no mitocondriales: *DEGs_PICs_AM_final_noMith.tsv* y *DEGs_EMTPICs_EMT_final_noMith.tsv*

***Output 4.3a e Input 4.3b***:
* Directorio *counts_mtx/* conteniendo *matrix.mtx*, *features.tsv* y *barcodes.tsv*
* Fichero con los genes tal y como pide cellphoneDB: *human-cellphone_genes-EMTp_x_AM.txt*

***Output 4.3b:***
Se obtienen dos directorios de resultados, pero empleamos:
* *out_subsampling*: Conteniendo entre otros *relevant_interactions.txt* con las 14 interacciones significativas.

## PASO 5
### PRIORIZACIÓN DE INTERACCIONES
En el script 5, se analizan los resultados obtenidos y se lleva a cabo una priorización final de las interacciones obtenidas.

*** Input:***
* Resultados de NicheNetR: *prioritized_tbl_oi_filtered_top.tsv*, *prioritized_tbl_oi-EMT_low.tsv*, *prioritized_tbl_oi-Singlets.tsv*
* Firmas de D. Cook *et al.*, 2022: *EMP_conserved_signature.tsv* y *Malignant_cell_signature.tsv*
* Resultados de CellphoneDB: *relevant_interactions.txt*

***Output:***
* Fichero conteniendo el ranking final ordenado de interacciones: *final_interaction_ranking.tsv*

**FIGURAS:**
* *Figura 11A-B*
* *Figura 12A-B*
* *Figura 13*
* 


