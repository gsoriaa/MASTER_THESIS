## PRIORIZACIÓN DE INTERACCIONES
## 
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(stringr)
library(openxlsx)
library(ggvenn)
library(pheatmap)
library(nichenetr)
# enrichment
library(clusterProfiler)
library(topGO)
library(org.Mm.eg.db) #mus musculus


##########################
# Función para escalar los resultados (scores) de la FIG.13: Rango de 0 a 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
##########################
# Función para traducir los genes de humano a ratón (EMT conserved signatures)
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
##########################
## Cargamos el fichero de interacciones obtenido
prioritized_tbl_filter <- read.csv(
  "prioritized_tbl_oi_filtered_top.tsv", 
  header = TRUE, sep = "\t") #66
receptor <- unique(prioritized_tbl_filter$receptor) # 34 unique receptor
write.table(receptor, "receptors_filtered_ligand.tsv", quote = FALSE, row.names = FALSE)
ligands<-unique(prioritized_tbl_filter$ligand) #42 unique ligands

# Análisis de promiscuidad de ligandos y receptores

AM_ligands_table <- as.data.frame(table(prioritized_tbl_filter$ligand)) %>% arrange(desc(Freq))
AM_ligands_table$Freq = as.numeric(AM_ligands_table$Freq)
TUM_receptors_table <- as.data.frame(table(prioritized_tbl_filter$receptor)) %>% arrange(desc(Freq))
TUM_receptors_table$Freq = as.numeric(TUM_receptors_table$Freq)

pdf("FIGURA_11A-Promiscuidad_elementos.pdf")
par(mfrow=c(2,1))
coul <- brewer.pal(8, "Set2") 
barplot(height = AM_ligands_table$Freq[1:10], names = AM_ligands_table$Var1[1:10],
        col = coul,
        horiz = T, las = 1, xlim = c(0,12), xlab = 'Number of interactions',
        main = "Occurences of each ligand")

barplot(height = TUM_receptors_table$Freq[1:10], names = TUM_receptors_table$Var1[1:10],
        col = coul,
        horiz = T, las = 1, xlim = c(0,12), xlab = 'Number of interactions',
        main = "Occurences of each receptor")
dev.off()
######################
####ENRICHMENT LIGAND-RECEPTORS??
####
ego_LR <- enrichGO(gene          = unique(c(receptor, ligands)),
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

# how many significant Biological Processes are there?
dim(ego_LR) # 1287 if p-val = 0.05

pdf("FIGURA_11B-LR_enrichment_T25.pdf")
dotplot(ego_LR, showCategory=nrow(ego_LR[0:25]))
dev.off()

################################
##########L-R IN EMTp SIGNATURE??
EMP_conserved_signature <- as.vector(as.matrix(
  read.table("EMP_conserved_signature.tsv", header = F)))
malignant_conserved_signature <- as.vector(as.matrix(
  read.table("Malignant_cell_signature.tsv", header = F)))

length(EMP_conserved_signature) #328 genes en la firma general conservada de EMP
length(malignant_conserved_signature) #128 genes en la firma conservada específica de célula tumoral

EMP_conserved_signature <- convertHumantoMouse(toupper(EMP_conserved_signature))
malignant_conserved_signature <- convertHumantoMouse(toupper(malignant_conserved_signature))

length(EMP_conserved_signature) #346 genes en la firma general conservada de EMP (30 no traducidos)
length(malignant_conserved_signature) #121 genes en la firma conservada específica de célula tumoral (11 no traducidos)

length(unique(c(receptor,ligands))) #74 unique elements
sum(c(receptor,ligands) %in% EMP_conserved_signature)/length(unique(c(receptor,ligands))) * 100 #8 common genes out of 88 L-R, 10.81 % of my elements
sum(c(receptor,ligands) %in% malignant_conserved_signature)/length(unique(c(receptor,ligands)))*100 #3 genes, 4.05%

c(receptor,ligands)[c(receptor,ligands) %in% EMP_conserved_signature]
##interactions involving: "F3"     "Itgb1" "Cdh1" "Grn"    "Lgals3" "Rtn4"   "Nampt"  "Efna1"  MIGHT BE DETERMINANT
c(receptor,ligands)[c(receptor,ligands) %in% malignant_conserved_signature]
#F3, "Cdh1, rtn4

#################################################
##COMPARE EMT-RICH VS LOW NICHENET INTERACTIONS##
#################################################
dim(prioritized_tbl_filter) #66 interactions in EMT_rich
prioritized_tbl_EMTlow <- read.csv("./prioritized_tbl_oi-EMT_low.tsv", header = TRUE, sep = "\t")
dim(prioritized_tbl_EMTlow) #95 interactions in EMT_low

##filter common interactions in both niches:
##(non EMT specific interactions)
prioritized_tbl_filter[prioritized_tbl_filter$ligand_receptor %in% prioritized_tbl_EMTlow$ligand_receptor,]
## 0 common interactions
pdf("FIGURA_12A-ggvenn_common-rich_low-interact.pdf")
PIC_LR <- list(
  "EMT_rich_interactions" = prioritized_tbl_filter$ligand_receptor, 
  "EMT_low_interactions" = prioritized_tbl_EMTlow$ligand_receptor)
ggvenn(PIC_LR,fill_color=c("cyan","burlywood2"), show_percentage = FALSE)
dev.off()

#################################################
##COMPARE WITH SINGLETS INTERACTIONS#############
#################################################

dim(prioritized_tbl_filter) #66 interactions in EMT_rich
prioritized_tbl_singlets <- read.csv("./prioritized_tbl_oi-Singlets.tsv", header = TRUE, sep = "\t")
dim(prioritized_tbl_singlets) #92 interactions in EMT_low

##filter common interactions in both niches:
##(non EMT specific interactions)
sum(prioritized_tbl_filter$ligand_receptor %in% prioritized_tbl_singlets$ligand_receptor)/66*100
##62/66common interactions!! (93.94%)
pdf("FIGURA_12B-ggvenn-singlet_interactions.pdf")
EMT_LR <- list(
  "EMT_PIC_interactions" = prioritized_tbl_filter$ligand_receptor, 
  "EMT_SING_interactions" = prioritized_tbl_singlets$ligand_receptor)
ggvenn(EMT_LR,fill_color=c("cyan","burlywood2"), show_percentage = FALSE)
dev.off()
###########################################################################################################################
############################PRIORITIZATION SCORE FINAL ####################################################################
###########################################################################################################################
# 
# table to score (copy from prioritized_tbl_filter)
table_prioritized <- prioritized_tbl_filter %>% dplyr::select(c("ligand", "receptor"))
rownames(table_prioritized) <- prioritized_tbl_filter$ligand_receptor
###################
###################
##1- prioritization score from NICHENET:
prior_score = range01(prioritized_tbl_filter$prioritization_score)
##2- EMT_general: EMP_conserved_signature$Nompropio
# sum(c(receptor,ligands_final) %in% EMP_conserved_signature$Nompropio)/73 * 100 #8 common genes out of 88 L-R, 10.96% of my elements
# interactions involving: F3, ITGB1, CDH1, CALM1, ANXA1, IL11, JAG1, DUSP18 MIGHT BE DETERMINANT
EMP_general_score <- c()
for (interaction in rownames(table_prioritized)){
  EMP_general_score = c(EMP_general_score, length(grep(paste0('^',table_prioritized[interaction, 'ligand'],
                                                              '$|^', table_prioritized[interaction, 'receptor'], '$'),
                                                       ignore.case = T, 
                                                       unique(EMP_conserved_signature), value = T)))
} # find genes in general signature and score it
EMP_general_score <- range01(EMP_general_score)
###################
## 3-EMT_specific: malignant_conserved_signature$nompropio
# sum(c(receptor,ligands_final) %in% malignant_conserved_signature$nompropio)/73*100 #3 genes, 4.1%
# F3, CDH1, RTN4
EMP_malignant_score <- c()
for (interaction in rownames(table_prioritized)){
  EMP_malignant_score = c(EMP_malignant_score, length(grep(paste0('^',table_prioritized[interaction, 'ligand'],
                                                                  '$|^', table_prioritized[interaction, 'receptor'], '$'),
                                                           ignore.case = T, 
                                                           unique(malignant_conserved_signature), value = T)))
} # find genes in general signature and score it
EMP_malignant_score <- range01(EMP_malignant_score)
###################
## 4- frequency: AM_ligands_table & TUM_receptors_table contains data of L-R frequency
freq_score <- c()
for (interaction in rownames(table_prioritized)){
  freq_score = c(freq_score, sum(TUM_receptors_table %>% filter(Var1 == table_prioritized[interaction, 'receptor']) %>% dplyr::select(Freq),
                                 AM_ligands_table %>% filter(Var1 == table_prioritized[interaction, 'ligand']) %>% dplyr::select(Freq)))
}
freq_score = range01(freq_score)
###################
## 5- NICHENET + CELLPHONE in EMT: 
cellphone_res =  read.table("./out_subsampling/relevant_interactions.txt", sep = "\t", header = T)
cellphone_res_emt = cellphone_res %>% 
  filter(AMs.EMT_PICs > 0 | AMs.EMT_singlet > 0 | 
           EMT_PICs.EMT_PICs > 0 | EMT_singlet.EMT_PICs > 0 | 
           EMT_PICs.AMs > 0 | EMT_PICs.EMT_singlet) %>% dplyr::select(interacting_pair)
match_score <- c()
for (interaction in rownames(table_prioritized)){
  # print(interaction)
  int = gsub('--', '_', interaction)
  # print(int)
  if (toupper(int) %in% cellphone_res_emt$interacting_pair){
    match_score = c(match_score, 1)} else{match_score = c(match_score, 0)}
}
match_score = range01(match_score)
######
## Y qué interacciones no son identificadas específicamente en EMT?
cellphone_res %>% filter(!interacting_pair %in% cellphone_res_emt$interacting_pair) %>% dplyr::select(interacting_pair)
#NRP1_VEGFA, PVR_NECTIN3.
###################
###################
# FINAL SCORE: weight each score & sum all values:
a = 0.25 # prioritization nichenet weight
b = 0.1 # general EMT = 10% weight
c = 0.15 # malignant EMT = 20% weight
d = 0.2 # PROMISQ = 25% weight
e = 0.3 # nichenet+cellphone = 45% weight
# And get final score
FINAL_SCORE = prior_score * a +
  EMP_general_score * b +
  EMP_malignant_score * c + 
  freq_score * d +
  match_score * e
# heatmap
FINAL_SCORE = range01(FINAL_SCORE)
Breaks<-c(seq(0,1, length=6)) #set breaks

cols <- rev(brewer.pal(6,"RdBu")) #set color for breaks (rev = reverse)
# write full pair ranking
'
manual_prioritization = data.frame("nichenet_prioritization" = prior_score, 
                                   "EMP_general_signature" = EMP_general_score,
                                   "EMP_malignant_signature" = EMP_malignant_score,
                                   "Promisquity" = freq_score,
                                   "cellphone_predict" = match_score,
                                   "final_score" = FINAL_SCORE, 
                                   row.names = prioritized_tbl_filter$ligand_receptor) %>%
  arrange(desc(FINAL_SCORE))  
write.table(manual_prioritization, "final_interaction_ranking.tsv", quote=FALSE, row.names = T, sep = "\t")
'
# plot top 30
manual_prioritization = data.frame("nichenet_prioritization" = prior_score, 
                                   "EMP_general_signature" = EMP_general_score,
                                   "EMP_malignant_signature" = EMP_malignant_score,
                                   "Promisquity" = freq_score,
                                   "cellphone_predict" = match_score,
                                   "final_score" = FINAL_SCORE, 
                                   row.names = prioritized_tbl_filter$ligand_receptor) %>%
  arrange(desc(FINAL_SCORE))  %>% top_n(30, FINAL_SCORE)
pdf("FIGURA_13-FINAL_SCORE.pdf")
pheatmap(t(manual_prioritization), color = cols, breaks = Breaks, 
         cluster_rows = F, cluster_cols = F, show_colnames = T, border_color = "grey",
         cellwidth = 10, cellheight = 10, 
         labels_row = c("Priorización por Nichenetr", "Firma de EMP general",
                        "Firma de EMP en cáncer", "Promiscuidad", "Predicho en cellphone", 
                        "Priorización global"))

dev.off()
'
head(manual_prioritization)
#                   nichenet_prioritization EMP_general_signature EMP_malignant_signature Promisquity cellphone_predict final_score
#Spp1--Cd44                       0.9362196                   0.0                       0  0.45454545                 1 1.0000000
#Lgals3--Itgb1                    0.8228880                   1.0                       0  1.00000000                 0 0.8092018
#Spp1--Itgb1                      0.9492709                   0.5                       0  1.00000000                 0 0.7797533
#Tnfsf12--Tnfrsf12a               0.9458264                   0.5                       1  0.09090909                 0 0.7274634
#Rtn4--Cd44                       0.6347269                   0.5                       1  0.45454545                 0 0.7193867
#F7--F3                           0.9977706                   0.5                       1  0.00000000                 0 0.7191496                 
'
