## Hago un merge de los tres ficheros post-qc.rds 
## (cada uno en su directorio out_*/seurat/merged/3_clustering/ ),
## sobreescribiendo el último post-qc.rds obtenido, para continuar 
## con el proceso de análisis de bollito
## 
## ## Desde R:
## Seurat(v.3.2.3) *WARNING: LA v.4.0 no es compatible con snakemake, por lo que tenemos que instalar la v3 para crear este objeto seurat
library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyverse)
combined<-merge(x = readRDS("./out_AMs/seurat/merged/1_preprocessing/seurat_post-qc.rds"), 
                y = c(readRDS("./out_tumor/seurat/merged/1_preprocessing/seurat_post-qc.rds"), 
                      readRDS("./out_PIC/seurat/merged/1_preprocessing/seurat_post-qc.rds")))

##AM.rds = 496 samples; TUM.rds = 971 samples; PIC.rds = 551
##merged.rds = 2018 samples

saveRDS(combined, "./out_AMs/seurat/merged/1_preprocessing/seurat_post-qc.rds") 
# en mi caso, el AM es el último en ser analizado, por eso se sobreescribe el rds de este directorio