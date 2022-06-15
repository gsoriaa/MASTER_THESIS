###FILTERING LIGANDS NOT EXPRESSED IN NAIVE AM
######################################
library(dplyr)
library(openxlsx)

# Búsqueda manual de información de los ligandos de interés en el browser de immgen: 
# http://rstats.immgen.org/Skyline/skyline.html 
ligands =  read.csv("best_upstream_ligands-immgen.tsv", sep = "\t", header = FALSE)
colnames(ligands) = c("gene", "extra_info")
out_ligands = ligands[ligands$info == "no expression lung",]$gene ##ligands filtered out through immgen browser
length(out_ligands) #18

# 18 ligandos aparecen como no expresados en AMs según immgen,
# se procede a analizar su expresión en un dataset de referencia:
AM_expression_data = read.xlsx("./TRMs_cancer-Maria-data/SUCIO-SUPPLEMENTARY_TABLE_3-ALL_DATA.xlsx")
## AM_Expression_data contiene los datos de expresión de AMs del paper de M. Casanova et al., 2021 (ST.3)

AM_expression_ligands = AM_expression_data[AM_expression_data$X1 %in% out_ligands,] #10 ligands identificados en la matriz

ligands_toKeep <- c() # ligandos que no quiero eliminar de mi dataset
for (ligand in AM_expression_ligands$X1){
  if (sum(AM_expression_ligands[AM_expression_ligands$X1 == ligand,2:9] > 5) >= 7 ){
    ligands_toKeep <- c(ligands_toKeep,ligand)
  }
} 
# Para aquellos genes que habrían sido eliminado por immgen, comprobamos aquellos que presenten
# más de 5 TPMs de expresión en al menos 7 de las 8 muestras .
# 8/18 genes se mantienen de acuerdo a los datos de María

out_ligands[out_ligands %in% ligands_toKeep ==FALSE]

# También se mantienen ADAM12 & DLL1 ya que podrían ser relevantes para el estudio
ligands_toKeep<-c(ligands_toKeep, "Dll1", "Adam12")
ligands_final <- ligands[ligands$gene %in% out_ligands == FALSE | ligands$gene %in% ligands_toKeep ==TRUE,"gene"]
length(ligands_final)
write.table(ligands_final, "best_upstream_ligands_filtered.tsv", quote=FALSE, row.names = FALSE, sep = "\t")
