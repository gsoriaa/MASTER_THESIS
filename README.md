# TRABAJO FIN DE MÁSTER -- Gonzalo Soria
Repositorio conteniendo los **scripts** y algunos de los **ficheros generados** en el análisis llevado a cabo durante el TFM realizado por Gonzalo Soria.

## PASO 4: Análisis de interacción célula-célula
### NicheNetR Analysis: scripts 4.1 y 4.2
En el script 4.1 se especifica la metodología a seguir para realizar el análisis empleando *NicheNetR* (v.1.0.0).
En el script 4.2 se realiza un filtrado de las interacciones obteniads en el análisis durante 4.1, que se filtran de acuerdo a la expresión de los ligandos identificados en datos de *AM* generados en **Maria Casanova *et al.,* 2021** para comprobar que estos ligandos tienen un mínimo de expresión en *AMs*. 
A continuación, el fichero aquí generado se carga de nuevo en el script 4.1 para continuar con el análisis.

***Input***:
* seurat data object: *seurat_find-clusters.rds*
* Fichero *combined_res.tsv*
* Fichero de datos de expresión de genes en *AMs* de Maria Casanova *et al.,* 2021: "./TRMs_cancer-Maria-data/SUPPLEMENTARY_TABLE_3-ALL_DATA.xlsx"
* (Fichero *best_upstream_ligands.tsv* generado en el script 4.1 y que se requiere para el script 4.2)
* (Fichero *best_upstream_ligands_filtered.tsv* generado en el script 4.2 y que se requiere para continuar el script 4.1)

***Output***:
* (Fichero *best_upstream_ligands.tsv* generado en el script 4.1 y que se requiere para el script 4.2)
* (Fichero *best_upstream_ligands_filtered.tsv* generado en el script 4.2 y que se requiere para continuar el script 4.1)
* Fichero *prioritized_tbl_oi_filtered_top.tsv* conteniendo la información de las 66 interacciones de interés predichas para el nicho de *EMT*.
* Fichero *prioritized_tbl_oi-EMT_low.tsv* conteniendo la información de las 95 interacciones predichas para el nicho pobre en *EMT*.
* Fichero *prioritized_tbl_oi-Singlets.tsv* conteniendo la información de las 92 interacciones predichas para los singletes en lugar de dobletes.

**FIGURAS:**
* *Figura 10A-D* (script 4.1)
* *Figura 11C* (script 4.1)

### CellphoneDB Analysis:
En el script 4.3a se generan los ficheros necesarios para ejecutar cellphoneDB desde R.
En el script 4.3b, se ejecuta cellphone desde el terminal de linux

***Input 4.3a***: 
* seurat data object: *seurat_find-clusters.rds*
* Ficheros de genes no mitocondriales: *DEGs_PICs_AM_final_noMith.tsv* y *DEGs_EMTPICs_EMT_final_noMith.tsv*

***Output 4.3a e Input 4.3b***:
* Directorio *counts_mtx/* conteniendo *matrix.mtx*, *features.tsv* y *barcodes.tsv*
* Fichero con los genes tal y como pide cellphoneDB: *human-cellphone_genes-EMTp_x_AM.txt*

***Output 4.3b:***
* *out_subsampling*: Conteniendo entre otros *relevant_interactions.txt* con las 14 interacciones significativas.
