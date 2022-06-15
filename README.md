## PASO 1: Generación de los ficheros iniciales
En el primer paso del análisis, se ejecuta el script 1.1 desde el software de R (v. 4.2.0) para eliminar las células T de mi análisis, y generar una visión global de los datos iniciales proporcionados. Además, se divide el fichero de metadatos proporcionado inicialmente generando un fichero de metadatos por cada uno de los ficheros de datos iniciales. Esta división es necesaria para ejecutar Bollito en el **PASO 2**

***Input:***
* Ficheros de cuentas: *AB1017\[1-8\].txt)*
* Ficheros de metadatos: *final_metadata.tsv*
* Ficheros de metacell (PIC-seq annotation): *PICs_metacell.txt ; singlets_metacell.txt*

***Output:***
* *AB1017\[1-8\]\_no\_Tcells.txt*
* *AB1017\[1-8\]\_no\_Tcells_metadata.txt*

**FIGURAS:** 
* *Figura 5A* (script 1.2) 
* *Figura 14* (script 1.1)
* *Figura Suplementaria 1* (script 1.1)
