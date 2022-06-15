# TRABAJO FIN DE MÁSTER -- Gonzalo Soria
Repositorio conteniendo los **scripts** y algunos de los **ficheros generados** en el análisis llevado a cabo durante el TFM realizado por Gonzalo Soria.

## PASO 3: Análisis del perfil transcripcional causado por la interacción
Una vez normalizados los datos, se procede al análisis transcriptómico *per se*.
Se emplea Seurat (V.4.0.0) de aquí en adelante para mí análisis.

En el script 3.1, se extraen los *DEGs* empleados en mi análisis.

En el script 3.2 se realiza el análisis de vías enriquecidas en estos *DEGs*.

***Input:***
* seurat data object: *seurat_find-clusters.rds*
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
