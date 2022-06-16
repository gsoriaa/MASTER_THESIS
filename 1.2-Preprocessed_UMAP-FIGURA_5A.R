################################################################
################################################################
#####
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(peRReo) # devtools::install_github("jbgb13/peRReo")
##
##
## Antes de realizar el preprocesado, se quiere tener una visión previa de los datos sin filtrar.
## Para ello, se genera un objeto a partir de todos los datos de expresión previos al empleo de bollito:
for (file in c("AB10171","AB10172","AB10173","AB10174","AB10175","AB10176","AB10177","AB10178")){ # para todos los ficheros
  count_df = read.table(paste0("./trimmed/counts/", file, "_no_Tcells.txt"))
  meta_df = read.table(paste0("./trimmed/metadata/",file, "_no_Tcells_metadata.txt"), header = T, sep = "\t")
  seurat = CreateSeuratObject(counts = count_df, meta.data = meta_df)
  saveRDS(seurat, paste0("./trimmed/",file, ".Rds"))
}
rds_full <- lapply(list.files(path = "./trimmed/", pattern = ".Rds"), function(x) readRDS(x))
preprocessed_assay <- do.call(merge, args = c(x = rds_full[1], y = list(rds_full[2:8])))

dim(preprocessed_assay@meta.data %>% 
      filter(Batch %in% 
               c("AB10171_no_Tcells","AB10174_no_Tcells" ))) #748
dim(preprocessed_assay@meta.data %>% 
      filter(Batch %in% 
               c('AB10172_no_Tcells','AB10175_no_Tcells', 'AB10176_no_Tcells'))) #1147
dim(preprocessed_assay@meta.data %>% 
      filter(Batch %in% 
               c('AB10173_no_Tcells','AB10177_no_Tcells', 'AB10178_no_Tcells'))) #877

# Visualización por UMAP
preprocessed_assay <- preprocessed_assay %>%
  FindVariableFeatures() %>% # %>% #2000 variable features
  ScaleData() %>% #requires scaled data
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE, umap.method = 'umap-learn')
pdf("preprocessed_assay_UMAP.pdf")
DimPlot(preprocessed_assay, label = TRUE, group.by= 'group', reduction='umap') + ggplot2::theme(legend.position = "none") #sin leyenda pero poniendo el nombre en el custer
DimPlot(preprocessed_assay, label = FALSE, group.by= 'group', reduction='umap') # con leyenda
dev.off()

## Para ver la mediana de Counts y features empleadas en bollito para cada uno de los 3 conjuntos de datos,
## y especificar los parámetros de los config.yaml mostrados en la tabla 1, podríamos emplear el siguiente código:
### Por ejemplo: en el caso de los AMs
                   
#a. N_features (genes)
median(as.vector(as.matrix(
  preprocessed_assay@meta.data %>% 
    filter(Batch %in% 
             c('AB10173_no_Tcells','AB10177_no_Tcells', 'AB10178_no_Tcells')) %>% 
    select(nFeature_RNA))))
#287 es la mediana
287/2.5 #114.8 como limite inferior
287*2.5 #717.5 como límite superior

#b. N_counts (moléculas)
median(as.vector(as.matrix(
  preprocessed_assay@meta.data %>% 
    filter(Batch %in% 
             c('AB10173_no_Tcells','AB10177_no_Tcells', 'AB10178_no_Tcells')) %>% 
    select(nCount_RNA))))

#591
591/3 #197 como limite inferior
591*3 #1773 como límite superior

# y esto mismo aplicado a los tres tipos de fichero editando los nombres del fichero de origen.
