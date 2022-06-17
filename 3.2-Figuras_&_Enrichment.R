library("ggrepel")
library("stringr")
library("ggvenn")
library(readxl)
library(ggplot2)
#
library(clusterProfiler); packageDescription ("clusterProfiler", fields = "Version") 
library(topGO)
#library(org.Hs.eg.db) #homo sapiens
library(org.Mm.eg.db) #mus musculus
library(stringr)

# FUNCTION TO TRANSLATE FROM HUMAN TO MOUSE GENES (EMT_SIGNATURES)
convertHumantoMouse <- function(x){
  
  require("biomaRt")
  # human = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  # mouse = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl") # these are giving mirror errors
  human = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  mousex <- unique(genesV2[, 2]) # unique the 2nd column values
  
  human_genes_number <- length(x)
  mouse_genes_number <- length(mousex)
  
  if(human_genes_number != mouse_genes_number){
    genes_not_trans <- setdiff(x, genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"),sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  
  return(mousex)
}
###########################1-
############################# CARGAR LOS DEGs de EMT_PICs vs. EMT_singlets

DEGs_EMTPICs_EMT<-read.table(file = "./DEGs_EMTPICs_EMT.tsv", sep = "\t", header = TRUE)
dim(DEGs_EMTPICs_EMT) #76

DEGs_EMTPICs_EMT_final<-read.table(file = "./DEGs_EMTPICs_EMT_res.tsv", sep = "\t", header = TRUE)
dim(DEGs_EMTPICs_EMT_final) #29

## ELIMINAMOS LOS GENES MITOCONDRIALES
DEGs_EMTPICs_EMT_final_noMith<-DEGs_EMTPICs_EMT_final[
  str_detect(rownames(DEGs_EMTPICs_EMT_final), 
             regex("^mt-", ignore_case = TRUE)) == FALSE,]

dim(DEGs_EMTPICs_EMT_final_noMith) #27

write.table(rownames(DEGs_EMTPICs_EMT_final_noMith), "genes_TUM_nomith.txt",
            quote=FALSE, row.names = FALSE)
write.table(DEGs_EMTPICs_EMT_final_noMith, "DEGs_EMTPICs_EMT_final_noMith.tsv", 
            quote=FALSE, row.names = TRUE, sep = "\t")

###VOLCANOPLOT DE ESTE DATASET (FIGURA 8A)
DEGs_EMTPICs_EMT$Type<-factor(DEGs_EMTPICs_EMT$Type, 
                              levels = c("Up-regulated", "Down-regulated"))
pdf("./FIGURA_8A-volcanoplot_DEGs_EMT.pdf")
p<- ggplot(data=DEGs_EMTPICs_EMT, 
           aes(x=avg_log2FC, y=-log10(p_val_adj), col=Type)) + 
  geom_point() + theme_minimal()
# Añadimos las líneas
p2 <- p + geom_vline(xintercept = c(-0.25,0.25),col="red", linetype = "dashed") + 
  geom_vline(xintercept = 0,col="red") + 
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggtitle("VolcanoPlot DEGs EMT_singlets vs EMT_PICs")
# Destacamos los genes filtrados finales
p3 <- p2 +   geom_point(data=DEGs_EMTPICs_EMT[rownames(DEGs_EMTPICs_EMT_final_noMith),], 
                        aes(x=avg_log2FC,y=-log10(p_val_adj)), 
                        color='gold',
                        pch=6,
                        size=4)

p3
#PLOT con nombres
####
ggplot(DEGs_EMTPICs_EMT_final_noMith,
       aes(x=avg_log2FC, y = -log(p_val_adj)),
       ylim())+
  geom_point() + # Show dots
  geom_text_repel(aes(label = rownames(DEGs_EMTPICs_EMT_final_noMith)),
                  min.segment.length = 0,
                  box.padding = 0.2,
                  max.overlaps = 15
  ) +
  ggtitle("EMT- PIC vs. singlets DEGs")
dev.off()
###########################2-
############################# CARGAMOS LOS DEGs de AMs vs. EMT_PICs

DEGs_PICs_AM<-read.table(file = "DEGs_PICs_AM.tsv", sep = "\t", header = TRUE)
dim(DEGs_PICs_AM) #749

DEGs_PICs_AM_final<-read.table(file = "DEGs_PICs_AM_res.tsv", sep = "\t", header = TRUE)
dim(DEGs_PICs_AM_final) #362

## ELIMINAMOS GENES MITOCONDRIALES
DEGs_PICs_AM_final_noMith<-DEGs_PICs_AM_final[
  str_detect(rownames(DEGs_PICs_AM_final),
             regex("^mt-", ignore_case = TRUE)) == FALSE,]

dim(DEGs_PICs_AM_final_noMith) #360

write.table(rownames(DEGs_PICs_AM_final_noMith), "genes_AM_nomith.txt",
            quote=FALSE, row.names = FALSE)
write.table(DEGs_PICs_AM_final_noMith, "DEGs_PICs_AM_final_noMith.tsv", 
            quote=FALSE, row.names = TRUE, sep = "\t")

###VOLCANOPLOT DE ESTE DATASET (FIGURA 8B)

pdf("./FIGURA_8B-volcanoplot_DEGs_AM.pdf")
DEGs_PICs_AM$Type<-factor(DEGs_PICs_AM$Type, levels = c("Up-regulated", "Down-regulated"))
p<- ggplot(data=DEGs_PICs_AM, aes(x=avg_log2FC, y=-log10(p_val_adj), col=Type)) +
  geom_point() + theme_minimal()

p2 <- p + geom_vline(xintercept = c(-0.25,0.25),col="red", linetype = "dashed") +
  geom_vline(xintercept = 0,col="red") + geom_hline(yintercept=-log10(0.05), col="red") +
  ggtitle("VolcanoPlot DEGs AM_singlets vs. EMT_PICs")

p3 <- p2 +   geom_point(data=DEGs_PICs_AM_final_noMith, 
                        aes(x=avg_log2FC,y=-log10(p_val_adj)), 
                        color='gold',
                        pch=6,
                        size=3)

p3

ggplot(DEGs_PICs_AM_final_noMith,
       aes(x=avg_log2FC, y = -log(p_val_adj)),
       ylim())+
  geom_point() + # Show dots
  geom_text_repel(aes(label = rownames(DEGs_PICs_AM_final_noMith)),
                  min.segment.length = 0,
                  box.padding = 0.2,
                  max.overlaps=20
  ) +
  ggtitle("AM_singlets VS. EMT_PIC DEGs")

dev.off()

###############################
###############################

###ANÁLISIS DE RUTAS ENRIQUECIDAS EN MI SET DE GENES DIFERENCIALMENTE EXPRESADOS: 
############### 
############### CLUSTERPROFILER
############### 

genes_full_combined = c(rownames(DEGs_EMTPICs_EMT_final_noMith), 
                        rownames(DEGs_PICs_AM_final_noMith))
genes_full_combined #387 genes

ego_full <- enrichGO(gene          = genes_full_combined, #genes
                     OrgDb         = org.Mm.eg.db, #organismo de referencia (ratón, mus musculus)
                     keyType       = 'SYMBOL', #código de los genes
                     ont           = "BP", # Ontología de referencia: Proceso biológico
                     pAdjustMethod = "BH", # Corrección por Benjamini - Hochberg
                     pvalueCutoff  = 0.05, # límite de p-valor 0.05
                     qvalueCutoff  = 0.05) # límite de q-value 0.05
#* En ont: BP, MF (Molecular function), CC (Cellular component) o ALL (para emplear todas las 
#* ontologías como referencia).

# Cuántos procesos se asocian significativamente
dim(ego_full) #636 if q-val = 0.05

# FIGURA 8C
pdf("FIGURA_8C-enrichment_dotplot_TOP25.pdf")
dotplot(ego_full, showCategory=nrow(ego_full[0:25])) 
dev.off()
##### Visualización de rutas relevante para el proceso de EMT
key_terms<- c("epith","epiderm","mesench","adhesion",
              "migration","angiogen","growth","hypoxia",
              "beta","p53","nfkb", "epidermal",
              "endothelial","extracel", "matrix", "kinase B",
              "tgf", "move", "chemotax", "remodel","grow","blood") 

# Me quedo con las categorías que contengan palabras clave
enriched_cathegories<-c()
for (term in key_terms){
  for (cathegory in ego_full@result$Description[ego_full@result$p.adjust < 0.05]){
    if (is.na(cathegory) == FALSE){
      if (str_detect(cathegory, term) == TRUE){
        enriched_cathegories<-c(enriched_cathegories, cathegory)
      }
    }
  }
}
ego_full_key<-dplyr::filter(ego_full,
                            ego_full@result$Description %in% enriched_cathegories)

# 72
pdf("FIGURA_8D-enrichment_dotplot_key.pdf")
dotplot(ego_full_key, showCategory=nrow(ego_full_key[0:30]))
dev.off()

################
################ FIRMAS DE EMP (Cook & Vanderhyden, 2022)
################

EMP_conserved_signature <- as.vector(as.matrix(
  read.table("EMP_conserved_signature.tsv", header = F)))
malignant_conserved_signature <- as.vector(as.matrix(
  read.table("Malignant_cell_signature.tsv", header = F)))

length(EMP_conserved_signature) #328 genes en la firma general conservada de EMP
length(malignant_conserved_signature) #128 genes en la firma conservada específica de célula tumoral
sum(EMP_conserved_signature %in% genes_full_combined)/328*100 #17 genes
sum(malignant_conserved_signature %in% genes_full_combined)/128*100 #5 genes, 4.13%

EMP_conserved_signature <- convertHumantoMouse(toupper(EMP_conserved_signature))
malignant_conserved_signature <- convertHumantoMouse(toupper(malignant_conserved_signature))

length(EMP_conserved_signature) #346 genes en la firma general conservada de EMP
length(malignant_conserved_signature) #121 genes en la firma conservada específica de célula tumoral
sum(EMP_conserved_signature %in% genes_full_combined)/346*100 #19 genes, 5.49%
sum(malignant_conserved_signature %in% genes_full_combined)/121*100 #5 genes, 4.13%
#  
#  
pdf("FIGURA_9-DCOOK_EMP_SIGNATURES-venn.pdf") #FIGURA 9
x <- list(
  "EMP_conserved_signature" = EMP_conserved_signature, 
  "DEGs_of_interest" = genes_full_combined)
ggvenn(x,fill_color=c("cyan","burlywood2"), show_percentage = FALSE, text_size = 5)

y <- list(
  "EMP_malignant" = malignant_conserved_signature, 
  "DEGs_of_interest" = genes_full_combined)
ggvenn(y,fill_color=c("cyan","burlywood2"), show_percentage = FALSE, text_size = 5)
dev.off()

# NEXT: INTERACTION ANALYSIS
