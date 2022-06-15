########
######## Realizar análisis por cellphone:
## DESDE R 
library(Seurat)
library(SeuratObject)
library(Matrix)
library(nichenetr) # únicamente para emplear la función de traducción de genes
library(dplyr)
#######################################################################################
#######################################################################################
### 1- crear una carpeta 'counts_mtx'
dir.create('./counts_mtx')
assay_seurat <- readRDS("./out_AMs/seurat/merged/3_clustering/seurat_find-clusters.rds")
Idents(assay_seurat) = assay_seurat$group
writeMM(assay_seurat@assays$SCT@counts, 
        file = 'counts_mtx/matrix.mtx') # guardamos la matriz
# NULL
write(x = rownames(assay_seurat@assays$SCT@counts), 
      file = "counts_mtx/features.tsv") #guardamos los feature names (genes)
write(x = colnames(assay_seurat@assays$SCT@counts), 
      file = "counts_mtx/barcodes.tsv") # guardamos los barcodes (células)
# generamos metadata
# table(assay_seurat@meta.data$group)
assay_seurat@meta.data$Cell = rownames(assay_seurat@meta.data)
df = assay_seurat@meta.data[, c('Cell', 'group')]
write.table(df, file ='meta.tsv', sep = '\t', quote = F, row.names = F)
# generamos el fichero de DEGs necesario para el análisis,
# especificando el cluster de mayor expresión: EMT_PICs ó AMs
# 
# se emplean como input los ficheros de DEGs 
# libres de genes mitocondriales
DEGs_PICs_AM_final_noMith = read.table("./DEGs_PICs_AM_final_noMith.tsv")
DEGs_EMTPICs_EMT_final_noMith = read.table("./DEGs_EMTPICs_EMT_final_noMith.tsv")

cellphone_genes_TPIC1<-data.frame(
  cluster = "EMT_PICs", 
  gene = rownames(
    DEGs_EMTPICs_EMT_final_noMith[
      DEGs_EMTPICs_EMT_final_noMith$avg_log2FC > 0,]
  )
) #los genes UP de EMT_PICs_vs_EMT para EMT_P
cellphone_genes_AM<-data.frame(
  cluster = "AMs", 
  gene = rownames(
    DEGs_PICs_AM_final_noMith[
      DEGs_PICs_AM_final_noMith$avg_log2FC < 0,]
  )
) # los genes down de PICs vs AM para AMs
cellphone_genes_TPIC2<-data.frame(
  cluster = "EMT_PICs", 
  gene = rownames(
    DEGs_PICs_AM_final_noMith[
      DEGs_PICs_AM_final_noMith$avg_log2FC > 0,]
  )
) # los genes UP de PICs vs AM a EMT_P
cellphone_genes = rbind(cellphone_genes_AM,
                        cellphone_genes_TPIC1, 
                        cellphone_genes_TPIC2) # 386 genes

# una vez creado el fichero de genes, es necesario traducir 
# los genes a humano, ya que es la referencia empleada por cellphone

cellphone_human_genes = cellphone_genes[!is.na(
  cellphone_genes$gene %>% 
           convert_mouse_to_human_symbols()),] # 346 genes
rownames(cellphone_human_genes) = 
  toupper(cellphone_human_genes$gene)
cellphone_human_genes$gene = 
  toupper(cellphone_human_genes$gene)

write.table(cellphone_human_genes, 
            "human-cellphone_genes-EMTp_x_AM.txt", 
            quote=FALSE, row.names = FALSE, sep = "\t")
# Y SE EJECUTA CELLPHONE DESDE EL ENTORNO DE CONDA EN EL TERMINAL 
# TAL Y COMO SE INDICA EN EL SCRIPT 4.3b