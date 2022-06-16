library(Seurat)
library(dplyr)
library(ggrepel)
library(SeuratObject)
##
### Leemos el objeto creado como resultado del filtrado y normalización de mis datos en Bollito
dir.create("./degs_analysis")
assay_seurat <- readRDS("./out_AMs/seurat/merged/3_clustering/seurat_find-clusters.rds")

#
## Para llevar a cabo el análisis de expresión diferencial (DEA) 
## primero debemos ver qué tipo de datos tenemos. Hay que comprobar
## qué distribución siguen mis datos para escoger el test adecuado.
#

##SHAPIRO-WILK TEST PARA VER LA DISTRIBUCIÓN:
shapiro.test(assay_seurat@meta.data$nFeature_SCT) # W = 0.95477, p-value < 2.2e-16
shapiro.test(assay_seurat@meta.data$nCount_SCT) # W = 0.93633, p-value < 2.2e-16
#p<0.05 por tanto es distribucion no normal

###Al no seguir una distribución normal, puedo aplicar el wilcoxon rank sum test como DEA
#*si p>0.05, la H0 de que mis datos siguen una distribución normal 
#no se puede descartar y por tanto tendría que emplear otro
#tipo de test estadístico (como un t-test) como método de DEA

## DEFINIR EL assay_seurat DE SCT
DefaultAssay(assay_seurat) <- "SCT"
# DEFINIR LA COLUMNA DE METADATA QUE SE EMPLEARÁ POR DEFECTO PARA DIFERENCIAR MIS CÉLULAS
Idents(assay_seurat) = assay_seurat@meta.data$group

## IDENTIFICAR LOS GENES MARCADORES QUE DEFINEN CADA UNO DE LOS CLUSTERS DE PICs
## FIGURA 6
MARKER_PICs <- FindAllMarkers(subset(assay_seurat, group %in% c("EMT_PICs", "PROLIF_PICs", "STRESS_PIC")),
                              test.use = "wilcox", group.by = "group")
dim(MARKER_PICs) #462
MARKER_PICs = MARKER_PICs %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC) > 0.45)
length(unique(MARKER_PICs$gene)) # 62 genes únicos
# Están los tres genes definitorios de cada cluster en mi set?
# 
c("Ccnd1", "Clu", "Hmga2") %in% MARKER_PICs$gene # TRUE TRUE TRUE

# Ordenamos los genes para representarlos de forma más visual en un heatmap
# 
stress_genes <- c()
prolif_genes <- c()
emt_genes <- c()
for (gene_id in unique(MARKER_PICs$gene)){
top_C = MARKER_PICs %>% filter(gene == gene_id) %>% arrange(desc(avg_log2FC)) %>% top_n(1, avg_log2FC)
#ordenar segun la expresion del gen y coger el top 1
if (top_C$cluster == "STRESS_PIC"){ #asignarlo al cluster que me interese
stress_genes <- c(stress_genes, top_C$gene)
}
if (top_C$cluster == "PROLIF_PICs"){ #asignarlo al cluster que me interese
prolif_genes <- c(prolif_genes, top_C$gene)
}
if (top_C$cluster == "EMT_PICs"){ #asignarlo al cluster que me interese
emt_genes <- c(emt_genes, top_C$gene)
}
}
marker_features <- c(emt_genes, prolif_genes, 
                     rev(stress_genes)) #rev() para cambiar el cluster de estrés y 
                                        # agrupar Ccnd1 con su cluster de expresión (proliferación)

pdf("FIGURA_6-PICs_marker_genes.pdf") 
DoHeatmap(subset(assay_seurat, group %in% c("EMT_PICs", "PROLIF_PICs", "STRESS_PIC")), 
          marker_features, group.by= "group")
dev.off()

#######################################################
########## IDENTIFICAR EMT_PIC DEGs ###################
#######################################################
####
####Aquí se procede a realizar el flujo esquematizado en la FIGURA 16
####
### PASO 1: DEFINIR LA PRIMERA APROXIMACIÓN DE DEGs
### MIN_LFC = 0.25
## A) TUMOR DEGs: EMT_S x EMT_P
DEGs_EMTPICs_EMT<- FindMarkers(assay_seurat, ident.1 = "EMT_PICs", ident.2 = "EMT_singlet", group.by = "group")
dim(DEGs_EMTPICs_EMT) #196 DEGs (sin filtrar)
Type<-c()
for (gene in rownames(DEGs_EMTPICs_EMT)) {
  if (DEGs_EMTPICs_EMT[gene, "avg_log2FC"] > 0){Type<-c(Type,"Up-regulated")} 
  else {Type<-c(Type,"Down-regulated")}
}
DEGs_EMTPICs_EMT$Type<-factor(Type, levels = c("Up-regulated", "Down-regulated"))
'
head(DEGs_EMTPICs_EMT)
#         p_val        avg_log2FC pct.1 pct.2    p_val_adj         Type
# Lyz2   4.535575e-42   3.825358 0.965 0.066 5.824132e-38 Up-regulated
# Chil3  3.784814e-36   2.789321 0.843 0.000 4.860080e-32 Up-regulated
# Ccl6   6.601848e-29   1.601593 0.713 0.000 8.477433e-25 Up-regulated
# Ctss   5.188125e-28   1.507657 0.696 0.000 6.662071e-24 Up-regulated
# Ear2   5.607122e-28   1.507657 0.696 0.000 7.200106e-24 Up-regulated
# Tmsb4x 1.481568e-24   1.727377 0.948 0.678 1.902482e-20 Up-regulated
'
write.table(DEGs_EMTPICs_EMT, "DEGs_EMTPICs_EMT_nofilter.tsv", quote=FALSE, row.names = TRUE, sep = "\t")
DEGs_EMTPICs_EMT<-dplyr::filter(DEGs_EMTPICs_EMT,
                                 p_val_adj < 0.05 )

write.table(DEGs_EMTPICs_EMT, "DEGs_EMTPICs_EMT.tsv", quote=FALSE, row.names = TRUE, sep = "\t")
dim(DEGs_EMTPICs_EMT) #76 DEGs significativos asociados a tumor (log2FC > 0.25)

##B) AM DEGs: AMs x EMT_P
DEGs_PICs_AM<-FindMarkers(assay_seurat, ident.1 = "EMT_PICs", ident.2 = "AMs", group.by = "group")
dim(DEGs_PICs_AM) #847
Type<-c()
for (gene in rownames(DEGs_PICs_AM)) {
  if(DEGs_PICs_AM[gene, "avg_log2FC"] > 0){Type<-c(Type,"Up-regulated")} 
  else {Type<-c(Type,"Down-regulated")}
}
DEGs_PICs_AM$Type<-factor(Type, levels = c("Up-regulated", "Down-regulated"))
write.table(DEGs_PICs_AM, "DEGs_PICs_AM_nofilter.tsv", quote=FALSE, row.names = TRUE, sep = "\t")
DEGs_PICs_AM<-dplyr::filter(DEGs_PICs_AM,
                                 p_val_adj < 0.05 )
write.table(DEGs_PICs_AM, "DEGs_PICs_AM.tsv", quote=FALSE, row.names = TRUE, sep = "\t")
dim(DEGs_PICs_AM) #749 DEGs significativos asociados a AM (log2FC > 0.25)

### STEP 2: FILTRAR AMBOS SETS DE GENES
## DEFINIMOS GENES BASALES QUE SE EXPRESAN MÁS EN UNO U
## OTRO TIPO DE CÉLULAS:  singletes de AM vs. singletes de tumor
DEGs_singlets <- FindMarkers(assay_seurat, ident.1 = c("EMT_singlet","PROLIF_singlet","STRESS_singlet"), ident.2 = "AMs", group.by = "group")
DEGs_singlets<-dplyr::filter(DEGs_singlets,
                                p_val_adj < 0.05 )
dim(DEGs_singlets) #[1] 829 5
DEGs_singlets_AM<-dplyr::filter(DEGs_singlets,
                                 avg_log2FC<0)
dim(DEGs_singlets_AM) #194
#
DEGs_singlets_TUM<-dplyr::filter(DEGs_singlets,
                                 avg_log2FC>0)
dim(DEGs_singlets_TUM) #635

#####PROBLEMA: los genes pueden aparecer como consecuencia del ruido de fondo
##### aportado por la segunda célula en interacción (y que no es de interés)

## Para empezar, hay 2 opciones: 
## 1- genes COMPARTIDOS entre ambos sets
## 2- genes ÚNICOS para uno de los sets

#1- Si los genes son COMUNES, no aparecen como consecuencia del ruido de fondo (al menos en uno de los casos)
## 4 OPCIONES:
## A.Si el gen se expresa más en ambos cluster de singlete que en el PIC:
#	| T_sing		   AM_sing
#	| 	      PIC			          => No podemos saber el origen del DEG, podría verse afectado 
# |______________________
## B.Si el gen se expresa menos en ambos cluster de singlete que en el PIC:
#|          PIC
#|	T_sing	  	 AM_sing			=> No podemos saber el origen del DEG, podría verse afectado 
#|_______________________
## C.Si el gen se expresa más en tumor pero menos en AM
#|	T_sing	  
#|		      PIC	
#|			         AM_sing			=> El gen se puede asociar a tumor
#|_______________________
## D.Si el gen se expresa menos en tumor pero más en AM
#|		  	       AM_sing
#|		      PIC	
#|	T_sing					          => El gen se puede asociar a AM
#|_______________________

############## Número de genes en cada opción:
up_DEGs_EMTPICs_EMT = DEGs_EMTPICs_EMT %>% filter(Type == 'Up-regulated') %>% rownames() # 71
down_DEGs_EMTPICs_EMT = DEGs_EMTPICs_EMT %>% filter(Type == 'Down-regulated') %>% rownames() #5 
up_DEGs_PICs_AM = DEGs_PICs_AM %>% filter(Type == 'Up-regulated') %>% rownames() #711
down_DEGs_PICs_AM = DEGs_PICs_AM %>% filter(Type == 'Down-regulated') %>% rownames() #38

sum(up_DEGs_EMTPICs_EMT %in% up_DEGs_PICs_AM) #0
sum(down_DEGs_EMTPICs_EMT %in% down_DEGs_PICs_AM) #0
sum(up_DEGs_EMTPICs_EMT %in% down_DEGs_PICs_AM) #26, 100% tumor genes
sum(down_DEGs_EMTPICs_EMT %in% up_DEGs_PICs_AM) #4, 100% AM genes


####SI LOS GENES NO SON COMUNES
# Se comprueba la expresión en singletes para filtrar estos genes 
# de acuerdo al tipo celular al que mayormente se puede asociar 
# su expresión  

#### EN RESUMEN, simplemente estamos eliminando genes
#### no comunes SOBREEXPRESADOS en los PICs que son:
#1- A-DEGs %in% AM_singlet genes (BG-noise-affected)
#2- B-DEGs %in% TUM_singlet genes (BG-noise-affected)

#1-
#Definimos el set de genes únicos asociados a AM sobre-expresados
# en el PIC:
#UP-GENES (71 genes)
#unique up-genes (45 genes)
unique_up_TUM = up_DEGs_EMTPICs_EMT[up_DEGs_EMTPICs_EMT %in% 
                                      rownames(DEGs_PICs_AM) == 
                                      FALSE]

#Definimos el set de genes eliminados que se asocian a AM erróneamente
trimmed_genes_TUM = unique_up_TUM[unique_up_TUM %in% 
                                    rownames(DEGs_singlets_AM)]

#filtramos los DEGs:
TUM_DEGs = rownames(DEGs_EMTPICs_EMT[
  rownames(DEGs_EMTPICs_EMT) %in% trimmed_genes_TUM == FALSE,])
length(TUM_DEGs) #33 genes

#2-
#Definimos el set de genes únicos asociados a tumor y sobre-expresados en el PIC:
#UP-GENES (711 genes)
#unique up-genes (707 genes)
unique_up_AM = up_DEGs_PICs_AM[up_DEGs_PICs_AM %in% rownames(DEGs_EMTPICs_EMT) == FALSE]

#Definimos el set de genes eliminados que se asocian erróneamente a tumor
trimmed_genes_AM = unique_up_AM[unique_up_AM %in% rownames(DEGs_singlets_TUM)]

#filtramos los DEGs:
AM_DEGs = rownames(DEGs_PICs_AM[rownames(DEGs_PICs_AM) %in% trimmed_genes_AM == FALSE,])
length(AM_DEGs) #388 genes

sum(TUM_DEGs %in% AM_DEGs) # 30 de los genes son comunes 
TUM_DEGs = TUM_DEGs[TUM_DEGs %in% up_DEGs_PICs_AM ==FALSE] #elimino los 4 genes que se asocian a AM
AM_DEGs = AM_DEGs[AM_DEGs %in% up_DEGs_EMTPICs_EMT ==FALSE] #elimino los 26 genes comunes que se asocian a tumor

length(TUM_DEGs) # 29
length(AM_DEGs) # 362

DEGs_EMTPICs_EMT_res = DEGs_EMTPICs_EMT[TUM_DEGs,] #DF
write.table(DEGs_EMTPICs_EMT_res, "DEGs_EMTPICs_EMT_res.tsv", quote=FALSE, row.names = TRUE, sep = "\t")
DEGs_PICs_AM_res = DEGs_PICs_AM[AM_DEGs,] #DF
write.table(DEGs_PICs_AM_res, "DEGs_PICs_AM_res.tsv", quote=FALSE, row.names = TRUE, sep = "\t")

## Combino los genes
combined_DEGs_names = unique(c(TUM_DEGs,AM_DEGs))
combined_DEGs <- rbind(DEGs_EMTPICs_EMT_res, DEGs_PICs_AM_res)


'
head(combined_DEGs, 10)
#              p_val avg_log2FC pct.1 pct.2    p_val_adj         Type
#Lyz2   4.535575e-42  3.8253583 0.965 0.066 5.824132e-38 Up-regulated
#Chil3  3.784814e-36  2.7893210 0.843 0.000 4.860080e-32 Up-regulated
#Ccl6   6.601848e-29  1.6015932 0.713 0.000 8.477433e-25 Up-regulated
#Ctss   5.188125e-28  1.5076568 0.696 0.000 6.662071e-24 Up-regulated
#Ear2   5.607122e-28  1.5076568 0.696 0.000 7.200106e-24 Up-regulated
#Tmsb4x 1.481568e-24  1.7273773 0.948 0.678 1.902482e-20 Up-regulated
#Laptm5 1.639314e-24  1.1934289 0.626 0.000 2.105043e-20 Up-regulated
#Lpl    4.076515e-22  1.0310269 0.574 0.000 5.234653e-18 Up-regulated
#Tyrobp 1.582895e-18  0.8687555 0.496 0.000 2.032596e-14 Up-regulated
#Clec7a 1.583651e-18  0.8824304 0.496 0.000 2.033566e-14 Up-regulated
'

write.table(combined_DEGs_names, "combined_DEGs_names_res.tsv", quote=FALSE, row.names = F, sep = "\t")
write.table(combined_DEGs, "combined_res.tsv", quote=FALSE, row.names = T, sep = "\t")
##
