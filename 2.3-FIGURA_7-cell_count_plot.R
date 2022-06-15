# pre- vs. post- object:
##2772 cells vs. 2018 cells (SCT)
##52634 features vs. 13460 features (SCT)
####################################################################################
## 
## FIGURA 7.  
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(peRReo) # devtools::install_github("jbgb13/peRReo")
##
assay_seurat <- readRDS("./out_AMs/seurat/merged/3_clustering/seurat_find-clusters.rds")
# pre- vs. post- object:
##2772 cells vs. 2018 cells (SCT)
##52634 features vs. 13460 features (SCT)
####################################################################################
table(assay_seurat@meta.data$group)
##	    AMs       EMT_PICs    EMT_singlet       PROLIF_PICs      PROLIF_singlet 
#           392            115            121            181            370 
##    STRESS_PIC STRESS_singlet     Unassigned 
#           103            464            272 
##
####PLOT
cell_counts = data.frame(Cluster = 
                           c("AMs","EMT_singlet","PROLIF_singlet", 
                             "STRESS_singlet", "EMT_PICs", 
                             "PROLIF_PICs", "STRESS_PICs"),
                         Cell_count = c(392, 121, 370, 
                                        464, 115, 
                                        181, 103))

pdf("FIGURA_7-cell_count.pdf")
# coul <- latin_palette("badbunny2", n = 8)
# coul <- latin_palette("badgyal", n = 8)
coul <- rev(latin_palette("karolg", n = 7))
barplot(height = cell_counts$Cell_count, names = cell_counts$Cluster,
        col = coul,
        horiz = F, las = 1, ylab = 'Número de células', ylim = c(0,500),
        main = "Recuento de estados celulares", width = 0.5)

# pie chart PICs
cell_perc = cell_counts$Cell_count[5:7]/sum(cell_counts$Cell_count[5:7]) * 100
pie(cell_perc, labels = paste0(round(cell_perc,2), "%"),
    col = coul[5:7], border = "white", main = "Proporción de dobletes identificados")

legend("bottomleft", legend = cell_counts$Cluster[5:7],
       fill =  coul[5:7])
dev.off()
####################################################################################

