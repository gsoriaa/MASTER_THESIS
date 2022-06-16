#######################
###FILTERING T-CELLS###
#######################
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
setwd("./data")
 
## data/raw con los ficheros originales y luego:
# dir.create("./trimmed")
# dir.create("./trimmed/counts")
# dir.create("./trimmed/metadata")

#### Filtrado de células T
#### Cargamos los ficheros de metadata:
#### Tenemos unos ficheros de metadata iniciales con los datos de las T_cells
singlets_metadata<-read.table("./raw/metadata/singlets_metacell.txt")
PICs_metadata<-read.table("./raw/metadata/pics_metacell.txt")
PICs_metadata$interaction <- paste0(PICs_metadata$a_name, '--', PICs_metadata$b_name)

##### Cargamos los ficheros de cuentas (UMIs)
AB10171<-read.table("raw/counts/AB10171.txt") #PIC
AB10172<-read.table("raw/counts/AB10172.txt")
AB10173<-read.table("raw/counts/AB10173.txt")
AB10174<-read.table("raw/counts/AB10174.txt") #PIC
AB10175<-read.table("raw/counts/AB10175.txt")
AB10176<-read.table("raw/counts/AB10176.txt")
AB10177<-read.table("raw/counts/AB10177.txt")
AB10178<-read.table("raw/counts/AB10178.txt")

Tcells_data<- c(singlets_metadata %>% filter(sin_names == "T_cells") %>% rownames(), # seleccionamos los singletes de células T
                rownames(PICs_metadata[grep("t_cell", PICs_metadata$interaction, ignore.case = T),])) # y los dobletes en los que intervienen

for (file in c("AB10171","AB10172","AB10173","AB10174","AB10175","AB10176","AB10177","AB10178")){ # para todos los ficheros, eliminamos células T
  file_df<-read.table(paste("./raw/counts/", as.character(file),".txt", sep = ""))
  file_trim<-file_df %>% select(colnames(.)[!colnames(.) %in% Tcells_data])
  write.table(file_trim, paste("./trimmed/counts/",as.character(file),"_no_Tcells.txt", sep=""), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
}

#### Leemos los ficheros generados libres de contaminación por células T
AB10171_nT<-read.table("./trimmed/counts/AB10171_no_Tcells.txt") #PIC
AB10172_nT<-read.table("./trimmed/counts/AB10172_no_Tcells.txt")
AB10173_nT<-read.table("./trimmed/counts/AB10173_no_Tcells.txt")
AB10174_nT<-read.table("./trimmed/counts/AB10174_no_Tcells.txt") #PIC
AB10175_nT<-read.table("./trimmed/counts/AB10175_no_Tcells.txt")
AB10176_nT<-read.table("./trimmed/counts/AB10176_no_Tcells.txt")
AB10177_nT<-read.table("./trimmed/counts/AB10177_no_Tcells.txt")
AB10178_nT<-read.table("./trimmed/counts/AB10178_no_Tcells.txt")

# gráfico de barras: recuento de células T ANTES del filtrado (FS. 1)
data<-data.frame(name = c("AB10171","AB10172","AB10173","AB10174","AB10175","AB10176","AB10177","AB10178"),
                 value = c(sum(colnames(AB10171) %in% Tcells_data),
                           sum(colnames(AB10172) %in% Tcells_data),
                           sum(colnames(AB10173) %in% Tcells_data),
                           sum(colnames(AB10174) %in% Tcells_data),
                           sum(colnames(AB10175) %in% Tcells_data),
                           sum(colnames(AB10176) %in% Tcells_data),
                           sum(colnames(AB10177) %in% Tcells_data),
                           sum(colnames(AB10178) %in% Tcells_data)))
pdf("FS1-Recuento_células_T.pdf")
par(mar=c(5.1,5.1,4.1,2.1))
coul <- brewer.pal(8, "Set2") 
barplot(height = data$value, names = data$name,
        col = coul,
        horiz = T, las = 1, xlim = c(0,250), xlab = 'T-cell count',
        main = paste("T-cell composition BEFORE trimming", sep = ""))
####
# gráfico de barras: recuento de células T DESPUÉS del filtrado (FS. 1)
data<-data.frame(name = c("AB10171","AB10172","AB10173","AB10174","AB10175","AB10176","AB10177","AB10178"),
                 value = c(sum(colnames(AB10171_nT) %in% Tcells_data),
                           sum(colnames(AB10172_nT) %in% Tcells_data),
                           sum(colnames(AB10173_nT) %in% Tcells_data),
                           sum(colnames(AB10174_nT) %in% Tcells_data),
                           sum(colnames(AB10175_nT) %in% Tcells_data),
                           sum(colnames(AB10176_nT) %in% Tcells_data),
                           sum(colnames(AB10177_nT) %in% Tcells_data),
                           sum(colnames(AB10178_nT) %in% Tcells_data)))
barplot(height = data$value, names = data$name,
        col = coul,
        horiz = T, las = 1, xlim = c(0,250), xlab = 'T-cell count',
        main = paste("T-cell composition AFTER trimming", sep = ""))
dev.off()
#### VISIÓN GENERAL DE LOS DATOS QUE TENEMOS?
#### 
## Cargamos metadata full
metadata_full<-read.table("./raw/metadata/final_metadata.tsv", sep = "\t", header = TRUE, row.names = 2)
dim(metadata_full) #2772 cells 
'
head(metadata_full)
           BBrowser_barcodes             Batch Total.count Total.expressed.feature     module      group  condition
WMC2607521      1_WMC2607521 AB10178_no_Tcells         525                     261  AM_MHC_hi        AMs AM_singlet
WMC2607522      1_WMC2607522 AB10178_no_Tcells         189                     121 Unassigned Unassigned Unassigned
WMC2607523      1_WMC2607523 AB10178_no_Tcells         333                     183 Unassigned Unassigned Unassigned
WMC2607524      1_WMC2607524 AB10178_no_Tcells        1582                     504         AM        AMs AM_singlet
WMC2607525      1_WMC2607525 AB10178_no_Tcells         507                     258         AM        AMs AM_singlet
WMC2607526      1_WMC2607526 AB10178_no_Tcells         624                     281         AM        AMs AM_singlet
'
for (file in c("AB10171","AB10172","AB10173","AB10174","AB10175","AB10176","AB10177","AB10178")){
  file_df<-read.table(paste("./trimmed/counts/",as.character(file),"_no_Tcells.txt", sep = ""))
  metadata_file<-metadata_full[colnames(file_df),]
  write.table(metadata_file, paste("./trimmed/metadata/",as.character(file),"_no_Tcells_metadata.txt", sep=""), row.names = T, col.names = TRUE, quote = FALSE, sep = "\t")
}

# Gráfico de composición original de cada fichero (fig. 14)
pdf("FIGURA_14-Composicion_por_fichero.pdf")
par(mfrow=c(4,2))
# Se reordenan los ficheros por motivos visuales: ficheros con mayormente AMs a la izquierda,
                                              # ficheros con mayormente células KP a la derecha, 
                                              # finalmente los dos ficheros de PICs
for (file in c("AB10173","AB10172","AB10177","AB10175","AB10178","AB10176","AB10171","AB10174")){ #files
  file_df<-read.table(paste("./trimmed/counts/",as.character(file),"_no_Tcells.txt", sep = ""), header = T)
  file_meta<-read.table(paste0("./trimmed/metadata/",as.character(file),"_no_Tcells_metadata.txt", sep = ""),  header = T)
  data<-data.frame(name = c("Tumor","AMs", "PIC", "Unlabeled"),
                   value = c(sum(colnames(file_df) %in% rownames(file_meta[str_detect(file_meta$condition, 
                                                                                      regex("Tumor", ignore_case = TRUE)),])),
                             sum(colnames(file_df) %in% rownames(file_meta[str_detect(file_meta$condition, 
                                                                                      regex("AM", ignore_case = TRUE)),])),
                             sum(colnames(file_df) %in% rownames(file_meta[str_detect(file_meta$condition, 
                                                                                      regex("PIC", ignore_case = TRUE)),])),
                             sum(colnames(file_df) %in% rownames(file_meta[str_detect(file_meta$condition, 
                                                                                      regex("Unassigned", ignore_case = TRUE)),])))
  )
  barplot(height = data$value, names = data$name,
          col = coul, xlim = c(0,400), xlab = 'Cell count',
          horiz = T, las = 1,
          main = paste(as.character(file)," cell composition", sep = "")
  )
}
dev.off()
####################################
####################################

######## SIGUIENTE PASO: Preprocesado por BOLLITO
