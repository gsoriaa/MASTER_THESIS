library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #
library(car)
## Adaptado de https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md
## 
###0. PREPARAMOS EL ANÁLISIS
##
assay_seurat <- readRDS("./out_AMs/seurat/merged/3_clustering/seurat_find-clusters.rds")
# setwd("./Nichenet")
## Definimos los nichos
fake_niches <- c()
for (cluster in assay_seurat@meta.data$group){
if (cluster == "PROLIF_PICs" || cluster == "STRESS_PIC"){
fake_niches <- c(fake_niches, "EMT_Low")} # nicho "pobre" en EMT
else{
if (cluster == "EMT_PICs"){
fake_niches <- c(fake_niches, "EMT_High")} # nicho "rico" en EMT
else{
if(cluster == "AMs"){
fake_niches <- c(fake_niches, "AMs")} # AMs 
else{
fake_niches <- c(fake_niches, "none")} # células que no me interesan
}
}
}
assay_seurat@meta.data$fake_niches <- fake_niches


celltype_id <- "fake_niches" # metadata column name of the cell type of interest
assay_seurat <- SetIdent(assay_seurat, value = assay_seurat[[celltype_id]])

## Leemos las matrices de referencia que asocian ligando y receptor (lr_network) y
## ligando y gen diana de interacción (ligand_target_matrix)
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network <- lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network <- lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
  distinct(ligand, receptor, bonafide)

## Traducimos los genes a ratón (la referencia es de humano)
organism = "mouse"
if(organism == "mouse"){
  lr_network = lr_network %>% 
    mutate(ligand = convert_human_to_mouse_symbols(ligand), 
           receptor = convert_human_to_mouse_symbols(receptor)) %>% drop_na()

  colnames(ligand_target_matrix) = ligand_target_matrix %>% 
    colnames() %>% convert_human_to_mouse_symbols()
  rownames(ligand_target_matrix) = ligand_target_matrix %>% 
    rownames() %>% convert_human_to_mouse_symbols()
  ligand_target_matrix = ligand_target_matrix %>% 
    .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
}

### Una vez definidos mis datos y la referencia, empieza el análisis:
##1- Definimos los nichos de interés
niches = list(
  "pEMT_rich_niche" = list(
    "sender" = c("AMs"), # expresa ligando
    "receiver" = c("EMT_High")), # expresa receptor
  "pEMT_poor_niche" = list(
    "sender" = c("AMs"),
    "receiver" = c("EMT_Low")
    )
  )

# definimos qué datos nos interesan
assay_oi = "SCT"

## 2- Se calcula la expresión diferencial entre nichos de interés
## tanto para sender como para receiver: Wilcoxon test
## 
###A. DE_sender: como en ambos casos empleo las mismas células,
### se diseña una dataframe 'mock' donde LFC = 0, y p_v= 1, 
### para que su valor no sea significativo
# DE_sender = calculate_niche_de(seurat_proj assay_seurat %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi)
DE_sender = tibble(data.frame(gene = rownames(assay_seurat %>% subset(features = lr_network$ligand %>% unique())), p_val = 1, avg_log2FC = 0, pct.1 = 0, pct.2 = 0, p_val_adj = 1, sender = "AMs", sender_other_niche = "AMs"))

# Definimos df diferencial para los receptores (EMT_High & EMT_Low)
DE_receiver = calculate_niche_de(seurat_obj = assay_seurat %>% 
                                   subset(features = lr_network$receptor %>% unique()), 
                                 niches = niches, type = "receiver", assay_oi = assay_oi)


expression_pct = 0.10 # determina el % mínimo de células expresando el gen (0.1 = 10%)
DE_sender_processed = process_niche_de(DE_table = DE_sender, 
                                       niches = niches, 
                                       expression_pct = expression_pct, 
                                       type = "sender")
DE_receiver_processed = process_niche_de(DE_table = DE_receiver, 
                                         niches = niches, 
                                         expression_pct = expression_pct, 
                                         type = "receiver")

# Se calcula la expresión diferencial de los ligandos y receptores
# diferenciando según tipo celular (no solo nicho).
# 
# Se procesa la información de ambas variables para definir el ranking de
# ligandos a emplear. Para ello, se puede emplear tanto el mínimo lfc como la media
# aunque se indica en la referencia que el mínimo sería recomendado al ser 
#   un criterio más restrictivo.
specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, 
                                                DE_receiver_processed, 
                                                lr_network, 
                                                specificity_score = specificity_score_LR_pairs)


## 3- Datos de expresión espacial: No disponible en este ensayo
include_spatial_info_sender = FALSE 
include_spatial_info_receiver = FALSE

# como antes, crearemos variables neutras que no influirán en el análisis
# A. Para mis células sender:
if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
}

if(include_spatial_info_sender == TRUE){ # en caso de que tuviesemos datos espaciales
  sender_spatial_DE = calculate_spatial_DE(seurat_proj = assay_seurat %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
  sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)

  # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
  sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)

  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))

} else { # Esto se ejecuta, añadiendo un score neutro para todas mis células.
  sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
  sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  

}

# B. Para mis células receiver:
if(include_spatial_info_receiver == TRUE){
  receiver_spatial_DE = calculate_spatial_DE(seurat_proj = assay_seurat %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)

  # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
  receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)

  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))

} else {
  receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
  receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}

##4- Se calcula la actividad del ligando y se infieren conexiones entre ligando y genes diana
lfc_cutoff = 0.15 # 0.15 es el recomendado, aunque si solo hay dos receiver niches se recomienda 
                  # aumentarlo a 0.25, pero como los cambios en el PIC se esperan más pequeños, se mantiene.
specificity_score_targets = "min_lfc"
expression_pct = 0.10

# Se define el conjunto de genes diana a emplear como los genes DE
# entre ambos nichos de células 'receiver0: Rico y pobre en EMT
DE_receiver_targets = calculate_niche_de_targets(seurat_obj = assay_seurat, 
                                                 niches = niches, 
                                                 lfc_cutoff = lfc_cutoff, 
                                                 expression_pct = expression_pct, 
                                                 assay_oi = assay_oi) 
##[1] "Calculate receiver DE between: rich and poor"
##  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01s  
##[1] "Calculate receiver DE between: poor and rich"

# Del conjunto de DE genes, conservo los que aparecen entre 
# mis 369 DEGs obtenidos en la etapa anterior
# 
DEGs_combined = read.table("./combined_res.tsv")
DE_receiver_targets = DE_receiver_targets %>% 
  dplyr::filter(gene %in% rownames(DEGs_combined)) #my DEGs

DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, 
                                                           niches = niches, expression_pct = expression_pct, 
                                                           specificity_score = specificity_score_targets)
# se define una variable de genes como ruido de fondo que
# comprende TODOS los genes identificados en uno u otro cluster
# de receiver cells
background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
geneset_niche1 = DE_receiver_processed_targets %>% 
  filter(receiver == niches[[1]]$receiver & 
           target_score >= lfc_cutoff & 
           target_significant == 1 & 
           target_present == 1) %>% 
  pull(target) %>% unique()
geneset_niche2 = DE_receiver_processed_targets %>% 
  filter(receiver == niches[[2]]$receiver & 
           target_score >= lfc_cutoff & 
           target_significant == 1 & 
           target_present == 1) %>% 
  pull(target) %>% unique()
  
# Se comprueban los genes que quedan fuera del análisis
# (no presentes en la ligand_target_matrix)
# Si se quedan fuera muchos genes, quizás haya un problema en cuanto
# al nombre de los genes (eg mapeo humano-ratón erróneo)
length(geneset_niche1)
## [1] 58
length(geneset_niche2)
## [1] 10
length(background)
## [1] 391
length(geneset_niche1 %>% setdiff(rownames(ligand_target_matrix)))
## [1] 9
length(geneset_niche2 %>% setdiff(rownames(ligand_target_matrix)))
## [1] 3
# Es aceptable, si me saliesen muchos genes en las dos ultimas 
# lineas, tendria que revisar que ha pasado con mis genes

# El número de genes a analizar es bastante reducido, aunque teniendo en
# cuenta la precisión del estudio es algo esperable.

top_n_target = 250

# se define mi lista a emplear en el análisis
niche_geneset_list = list(
  "pEMT_rich_niche"= list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "pEMT_poor_niche"= list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
  )

# y se calcula la actividad de ligando sobre genes diana  
ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, 
                                                          ligand_target_matrix = ligand_target_matrix, 
                                                          top_n_target = top_n_target)

##5.Se calcula la expresión ESCALADA de los ligandos, receptores y
## genes diana en cada nicho (expresión logarítmica y fracciones de expresión)

features_oi = union(lr_network$ligand, lr_network$receptor) %>% 
  union(ligand_activities_targets$target) %>% setdiff(NA)

# se representa la expresión por dotplot de seurat
dotplot = suppressWarnings(Seurat::DotPlot(
  assay_seurat %>% 
    subset(idents = niches %>% 
             unlist() %>% 
             unique()), 
  features = features_oi, 
  assay = assay_oi))
pdf("dotplot_step_5.pdf")
dotplot # visualización
dev.off()
# definimos df global
exprs_tbl = dotplot$data %>% as_tibble()
exprs_tbl = exprs_tbl %>% rename(
  celltype = id, gene = features.plot,
  expression = avg.exp, 
  expression_scaled = avg.exp.scaled, 
  fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% as_tibble() %>% 
  select(celltype, gene, expression, expression_scaled, fraction) %>%
  distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
# df de ligando
exprs_tbl_ligand = exprs_tbl %>% 
  filter(gene %in% lr_network$ligand) %>% 
  rename(sender = celltype, ligand = gene, 
         ligand_expression = expression, 
         ligand_expression_scaled = expression_scaled, 
         ligand_fraction = fraction) 

# df de receptor
exprs_tbl_receptor = exprs_tbl %>% 
  filter(gene %in% lr_network$receptor) %>% 
  rename(receiver = celltype, receptor = gene, 
         receptor_expression = expression, 
         receptor_expression_scaled = expression_scaled, 
         receptor_fraction = fraction)

# df de genes diana
exprs_tbl_target = exprs_tbl %>% 
  filter(gene %in% ligand_activities_targets$target) %>% 
  rename(receiver = celltype, target = gene, 
         target_expression = expression, 
         target_expression_scaled = expression_scaled, 
         target_fraction = fraction)

exprs_tbl_ligand = exprs_tbl_ligand %>%  
  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>%
  mutate(ligand_fraction_adapted = ligand_fraction) %>% 
  mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% 
  mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = exprs_tbl_receptor %>% 
  mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% 
  mutate(receptor_fraction_adapted = receptor_fraction) %>% 
  mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% 
  mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

##6. Aquí se puntua la interacción ligando-receptor en base a la 
## expresión del receptor. Se le dará un mayor score a los receptores
## más fuertemente expresados en uno u otro nicho.No afecta a la hora de
## priorizar ligandos, pero sí de priorizar receptores por ligando

exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% 
  inner_join(DE_sender_receiver %>% 
               distinct(niche, sender, receiver))
  
ligand_scaled_receptor_expression_fraction_df = 
  exprs_sender_receiver %>% 
  group_by(ligand, receiver) %>%
  mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 

##7. Prioritización de las conexiones ligando-receptor & ligando-gen diana
##

#variable para definir el peso que se da a cada parámetro creado a
#lo largo del análisis
prioritizing_weights = 
  c("scaled_ligand_score" = 0, #5 originariamente, se establece en 0 porque el sender es el mismo en ambos (AMs), mock DE_sender data
    "scaled_ligand_expression_scaled" = 1, # expresión del ligando en AMs (recomendado 1)
    "ligand_fraction" = 1, # prioriza ligandos expresados en mayor % de células (pero sin eliminar los que no superan el cutoff definido en 0.1)
    "scaled_ligand_score_spatial" = 0, #2 originariamente, se establece en 0 porque no tenemos datos espaciales
    "scaled_receptor_score" = 0.5, # igual que ligand pero para los receptores
    "scaled_receptor_expression_scaled" = 0.5,
    "receptor_fraction" = 1, 
    "ligand_scaled_receptor_expression_fraction" = 1,
    "scaled_receptor_score_spatial" = 0,
    "scaled_activity" = 0, # mejor ignorar valores absolutos de expresión
    "scaled_activity_normalized" = 1, #y centrarnos en los normalizados
    "bona_fide" = 1) # se le da mayor peso a las interacciones de bases de datos curadas

#* La red de interacciones ligando-receptor que emplea NichenetR como referencia
# consiste en interacciones documentadas en bases de datos curadas (definidas como 
# 'bona fide interactions'); e interacciones predichas en base a la 
# anotación del gen y  a interacciones proteína-proteína (PPIs)
# 

output = list(DE_sender_receiver = DE_sender_receiver, 
              ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
              sender_spatial_DE_processed = sender_spatial_DE_processed, 
              receiver_spatial_DE_processed = receiver_spatial_DE_processed,
              ligand_activities_targets = ligand_activities_targets,
              DE_receiver_processed_targets = DE_receiver_processed_targets,
              exprs_tbl_ligand = exprs_tbl_ligand,  
              exprs_tbl_receptor = exprs_tbl_receptor, 
              exprs_tbl_target = exprs_tbl_target)

prioritization_tables = get_prioritization_tables(output, prioritizing_weights)

###8. VISUALIZACIÓN DE LOS RESULTADOS
# Se define las interacciones L-R TOP
# 
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
  group_by(ligand) %>% 
  top_n(1, prioritization_score) %>% ungroup() %>% 
  select(ligand, receptor, niche) %>% 
  rename(top_niche = niche)

top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
  group_by(ligand, receptor) %>% 
  top_n(1, prioritization_score) %>% ungroup() %>%
  select(ligand, receptor, niche) %>% 
  rename(top_niche = niche)

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, prioritization_score) %>%
  group_by(ligand, niche) %>% 
  top_n(1, prioritization_score) %>% ungroup() %>% 
  distinct() %>% inner_join(top_ligand_niche_df) %>% 
  filter(niche == top_niche) %>% 
  group_by(niche) %>% 
  top_n(50, prioritization_score) %>%
  ungroup() # get the top50 ligands per niche

# Una vez obtenidos, se comprueban las interaccion L-R
# en las células 'receiver' de rich_emt (EMT_PICs), que son las
# que nos interesan:
receiver_oi = "EMT_High" 

filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% pull(ligand) %>% unique()

write.table(filtered_ligands, "best_upstream_ligands.tsv", quote = FALSE, row.names = FALSE, sep ="\t")

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
  distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% 
  group_by(ligand) %>% filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% ungroup() 

# 92 interacciones únicas asocciadas a mi nicho de interés

write.table(prioritized_tbl_oi, "prioritized_tbl_oi-singlets.tsv", quote=FALSE, row.names = FALSE, sep = "\t")

## A continuación, se lleva a cabo un filtro posterior, según aquellos
## ligandos que presenten mayor expresión en AMs.
## Se consulta la web immgen y el paper de Maria

******INSERTAR FILTRO*******
  
## SE ELIMINAN 8 LIGANDOS POCO RELEVANTES
## remove from original best upstream ligands, 
## the ones that are not relevant

filtered_ligands = as.vector(as.matrix(read.csv("./best_upstream_ligands_filtered.tsv", header = FALSE)))

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
  distinct() %>% 
  inner_join(top_ligand_receptor_niche_df) %>% 
  group_by(ligand) %>% filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% ungroup() 

# 78 interacciones únicas asocciadas a mi nicho de interés

write.table(prioritized_tbl_oi, "prioritized_tbl_oi_filtered.tsv", quote=FALSE, row.names = FALSE, sep = "\t")

## Y me quedo con las interacciones que se asocien a mi top_niche
## añadiendo ese ultimo filter(top_niche == niche)
prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
  distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% 
  group_by(ligand) %>% filter(receiver == receiver_oi) %>% 
  filter(top_niche == niche) %>% 
  top_n(2, prioritization_score) %>% ungroup() 

# 65 interacciones únicas asocciadas a mi nicho de interés

write.table(prioritized_tbl_oi, "prioritized_tbl_oi_filtered_top.tsv", quote=FALSE, row.names = FALSE, sep = "\t")

## PLOTs: 1- se visualiza una comparación de la expresión por medio del
## mínimo lfc observado en uno y otro nicho para cada ligando de interés:

lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, 
                                         prioritization_tables$prioritization_tbl_ligand_receptor, 
                                         plot_legend = FALSE, heights = NULL, widths = NULL)

## Visualization de la expresión de ligandos y 
## actividad de ligandos en relación a los links L-gen diana

######## Se corre el código crudo de las funciones preestablecidas
######## dado que provocaban errores (con algunas modificaciones)
# función editada: make_ligand_activity_target_exprs_plot = function(receiver_oi, prioritized_tbl_oi, prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_target, exprs_tbl_ligand, exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = TRUE, heights = NULL, widths = NULL)
best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()

# expresión de los ligandos
ordered_ligands = 
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% best_upstream_ligands) %>% 
  select(niche, sender, ligand, ligand_score) %>% 
  distinct() %>% group_by(ligand) %>% 
  summarise(ligand_score = max(ligand_score)) %>% 
  inner_join(prioritization_tables$prioritization_tbl_ligand_receptor %>% 
               select(niche, sender, ligand, ligand_score) %>% 
               distinct()) %>%
  arrange(sender, ligand_score)

ordered_ligands = ordered_ligands %>% 
  mutate(ligand_ordered = factor(ligand, ordered = T)) %>% 
  distinct(ligand, ligand_ordered, niche) %>% 
  rename(niche_prior = niche)

plot_data = output$exprs_tbl_ligand %>% 
  inner_join(ordered_ligands) %>% 
  filter(sender %in% (prioritization_tables$prioritization_tbl_ligand_receptor$sender 
                      %>% unique()))
plot_data = plot_data %>% group_by(ligand) %>% 
  mutate(ligand_expression_scaled_sender = nichenetr::scaling_zscore(ligand_expression)) %>% 
  inner_join(prioritization_tables$prioritization_tbl_ligand_receptor %>% 
               distinct(sender, receiver, niche))
p1 = plot_data  %>%
    # ggplot(aes(ligand_ordered, sender , color = ligand_expression_scaled_myeloid, size = ligand_fraction )) +
    ggplot(aes(sender,ligand_ordered , fill = ligand_expression_scaled_sender)) +
    # geom_point() +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.25, "lines"),
      panel.spacing.y = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Scaled Exprs ligand") + xlab("Ligand Expression") + ylab("Prioritized Ligands")
max_exprs = abs(plot_data$ligand_expression_scaled_sender ) %>% max()
custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
# custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

p_ligands = p1 + custom_scale_fill
# p_ligands: no me interesa este plot porque hace referencia a DE_sender
#
# expresión de los genes diana
# 
# como con los ligandos, pero con targets
targets_oi = prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(target_score >= lfc_cutoff) %>% 
  filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% 
  pull(target) %>% unique()

ordered_targets = prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(target %in% targets_oi) %>% 
  select(niche, receiver, target, target_score) %>% 
  distinct()  %>% arrange(receiver, -target_score)

# si hay duplicados:
ordered_targets = ordered_targets %>% select(-niche) %>% 
  distinct() %>% mutate(niche = receiver)
ordered_targets = ordered_targets %>% 
  mutate(target_ordered = factor(target, ordered = T, 
                                 levels = ordered_targets$target)) %>% 
  distinct(target, target_ordered, niche) %>% 
  rename(niche_prior = niche)

plot_data = output$exprs_tbl_target %>% 
  inner_join(ordered_targets) %>% 
  filter(receiver %in% (prioritization_tables$prioritization_tbl_ligand_target$receiver %>% 
                          unique()))
plot_data = plot_data %>% 
  group_by(target) %>% 
  mutate(target_expression_scaled_myeloid = nichenetr::scaling_zscore(target_expression))

p1 = plot_data %>%
    # ggplot(aes(target_ordered, receiver , color = target_expression_scaled_myeloid, size = target_fraction )) +
    ggplot(aes(target_ordered, receiver , fill = target_expression_scaled_myeloid)) +
    # geom_point() +
    geom_tile(color = "black") +
    # facet_grid(receiver~niche_prior, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x =  element_text(size = 9),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.background.x = element_blank(),
      strip.text.y = element_blank(),
      strip.text.x = element_blank()
    ) + labs(fill = "Scaled Exprs Target") + xlab("Target Expression")
max_exprs = abs(plot_data$target_expression_scaled_myeloid ) %>% max()
custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))
# custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu") %>% rev(),values = c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1),  limits = c(-1*max_exprs, max_exprs))

p_targets = p1 + custom_scale_fill
pdf("FIGURA_10C-plot_targets.pdf")
p_targets
dev.off()
# Ligand-Target heatmap
active_ligand_target_links_df = prioritization_tables$prioritization_tbl_ligand_target %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi) %>% 
  dplyr::select(ligand, target, ligand_target_weight ) %>% 
  dplyr::rename(weight = ligand_target_weight )

active_ligand_target_links_df = active_ligand_target_links_df %>%
  dplyr::filter(!is.na(weight))
if(active_ligand_target_links_df$target %>% 
   unique() %>% length() <= 2){
    cutoff = 0
    } else {
    cutoff = 0.33
}

active_ligand_target_links = 
  nichenetr::prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix, 
    cutoff = cutoff)

order_ligands_ = ordered_ligands$ligand_ordered %>% levels()
order_targets_ = ordered_targets$target_ordered %>% levels()

order_ligands = order_ligands_ %>% make.names()
order_targets = order_targets_ %>% make.names()

rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

if(length(setdiff(order_ligands, colnames(active_ligand_target_links))) > 0){
    removed_ligands = setdiff(order_ligands, colnames(active_ligand_target_links))
    new_lt_tibble = removed_ligands %>% 
      lapply(function(ligand_oi){
        tibble(ligand = ligand_oi, target = order_targets, weight = 0)
        }) %>% 
      bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(ligand, weight)
    active_ligand_target_links = new_lt_tibble %>% select(-target) %>% data.frame() %>% as.matrix(ncol = length(removed_ligands)) %>% cbind(active_ligand_target_links)
}
if(length(setdiff(order_targets, rownames(active_ligand_target_links))) > 0){
    removed_targets = setdiff(order_targets, rownames(active_ligand_target_links))
    new_lt_tibble = removed_targets %>% 
      lapply(function(target_oi){
        tibble(target = target_oi, ligand = order_ligands, weight = 0)
        }) %>% 
      bind_rows
    new_lt_tibble = new_lt_tibble %>% spread(target, weight)
    active_ligand_target_links = new_lt_tibble %>% select(-ligand) %>% data.frame() %>% as.matrix(ncol = length(removed_targets)) %>% t() %>% rbind(active_ligand_target_links)
}
if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = 
      active_ligand_target_links[order_targets,order_ligands] %>%
      matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
    } else {
    vis_ligand_target = 
      active_ligand_target_links[order_targets,order_ligands] %>%
      t()
}

p_ligand_target_network = vis_ligand_target %>% nichenetr::make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory\nPotential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")
pdf("FIGURA_10B-ligand_target_network.pdf")
p_ligand_target_network
dev.off()

if(is.null( prioritization_tables$prioritization_tbl_ligand_target %>% 
            pull(receiver) %>% levels())){
  order_receivers = prioritization_tables$prioritization_tbl_ligand_target %>% pull(receiver) %>% unique()
  } else{
    order_receivers = prioritization_tables$prioritization_tbl_ligand_target %>% pull(receiver) %>% levels()
  }

  # Actividad del ligando (escalada)
  # ligand_pearson_df = prioritization_tables$prioritization_tbl_ligand_target %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi)  %>% dplyr::select(ligand, niche, activity_normalized) %>% dplyr::distinct() %>% tidyr::spread(niche, activity_normalized)
ligand_pearson_df = 
  prioritization_tables$prioritization_tbl_ligand_target %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(ligand %in% best_upstream_ligands)  %>% 
  dplyr::select(ligand, receiver, activity_normalized) %>% 
  dplyr::distinct() %>% 
  tidyr::spread(receiver, activity_normalized)
  # print(ligand_pearson_df)
ligand_pearson_matrix = 
  ligand_pearson_df %>% dplyr::select(-ligand) %>%
  as.matrix() %>% 
  magrittr::set_rownames(ligand_pearson_df$ligand)
rownames(ligand_pearson_matrix) = 
  rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = 
  colnames(ligand_pearson_matrix) %>% make.names()
vis_ligand_pearson = 
  ligand_pearson_matrix[order_ligands %>% 
                          generics::intersect(rownames(ligand_pearson_matrix)), order_receivers  %>% 
                          make.names()] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% 
  nichenetr::make_heatmap_ggplot("Prioritized ligands",
                                 "Scaled Ligand activity", 
                                 color = "purple",
                                 legend_position = "top",
                                 x_axis_position = "top", 
                                 legend_title = "Scaled\nLigand\nActivity") + 
  theme(legend.text = element_text(size = 8))
#limits = c(min(vis_ligand_pearson, na.rm =TRUE), max(vis_ligand_pearson, na.rm =TRUE))
limits = c(-max(abs(vis_ligand_pearson), na.rm = TRUE), 
           max(abs(vis_ligand_pearson), na.rm = TRUE)
           )
# print(limits)
custom_scale_fill = scale_fill_gradientn(colours = c("white",RColorBrewer::brewer.pal(n = 7, name = "PuRd")),
                                         values = c(0, 0.50, 0.55, 0.625, 0.70, 0.80, 0.90, 1),  
                                         limits = limits)
p_ligand_pearson_scaled = p_ligand_pearson + custom_scale_fill
pdf("FIGURA_10A-pearson_ligand_scaled.pdf")
p_ligand_pearson_scaled
dev.off()
# Actividad del ligando
ligand_pearson_df = 
  prioritization_tables$prioritization_tbl_ligand_target %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(ligand %in% best_upstream_ligands)  %>% 
  dplyr::select(ligand, receiver, activity) %>% 
  dplyr::distinct() %>% tidyr::spread(receiver, activity)
ligand_pearson_matrix = 
  ligand_pearson_df %>% 
  dplyr::select(-ligand) %>% as.matrix() %>% 
  magrittr::set_rownames(ligand_pearson_df$ligand)
rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
vis_ligand_pearson = 
  ligand_pearson_matrix[order_ligands %>% 
                          generics::intersect(rownames(ligand_pearson_matrix)), order_receivers %>%
                          make.names()] #%>% as.matrix(ncol = 3) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% 
  nichenetr::make_heatmap_ggplot("Prioritized ligands",
                                 "Ligand activity",
                                 color = "darkorange",
                                 legend_position = "top",
                                 x_axis_position = "top", 
                                 legend_title = "Ligand\nActivity") + 
  theme(legend.text = element_text(size = 8))
custom_scale_fill = 
  scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "Oranges")),
                       values = c(0, 0.30, 0.40, 0.575, 0.70, 0.80, 0.925, 1),  
                       limits = c(min(vis_ligand_pearson, na.rm =TRUE), 
                                  max(vis_ligand_pearson, na.rm =TRUE)))
p_ligand_pearson = p_ligand_pearson + custom_scale_fill
pdf("pearson_ligand.pdf")
p_ligand_pearson
dev.off()
# Combinamos los plots (FIGURA 10 A,B,C)
n_groups = ncol(vis_ligand_pearson)
n_targets = ncol(vis_ligand_target)
n_ligands = nrow(vis_ligand_target)
n_senders = prioritization_tables$prioritization_tbl_ligand_receptor$sender %>% unique() %>% length()

legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 3) %>%
  patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))

heights = c(n_ligands, n_groups + 0.5)
widths = c(n_senders + 0.5, n_groups, n_groups, n_targets)
pdf("FIGURA_10-combined_plot.pdf")
design <- "Aa#B
		###C"
combined_plot = 
  patchwork::wrap_plots(A = p_ligand_pearson_scaled + 
                          theme(legend.position = "none", axis.ticks = element_blank()) + 
                          theme(axis.title.x = element_text()) + ylab(""),
                        a = p_ligand_pearson + 
                          theme(legend.position = "none", axis.ticks = element_blank()) +
                          ylab(""),
                        B = p_ligand_target_network + 
                          theme(legend.position = "none", axis.ticks = element_blank()) + 
                          ylab(""),
                        C = p_targets + 
                          theme(legend.position = "none"),
                        nrow = 2, design = design, 
                        widths = widths, heights = heights)
list(combined_plot = combined_plot, legends = legends)
dev.off()
dev.off()
########### DOTPLOT PARA LA EXPRESIÓN DE RECEPTORES
###########  (FIGURA 11C)

# Seurat dotplot
##celltype_id = "fake_niches"# metadata column name of the cell type of interest
##seurat_proj = SetIdent(assay_seurat, value = assay_seurat[[celltype_id]])

order_receptors_adapted = 
  sort(unique(prioritized_tbl_oi$receptor))
rotated_dotplot = DotPlot(assay_seurat %>% 
                            subset(fake_niches %in% 
                                     c("EMT_Low","EMT_High")), 
                          features = order_receptors_adapted, 
                          cols = "RdYlBu") + 
  coord_flip() + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12)) 
pdf("FIGURA_11B-receptor_expresion_dotplot.pdf")
rotated_dotplot
dev.off()
##CIRCOS PLOT (FIGURA 10 D.)

## Para que se vea mejor, se seleccionan únicamente el TOP-15 ligandos
## para su representación
filtered_ligands = ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  filter(ligand %in% filtered_ligands) %>% 
  top_n(15, prioritization_score) %>% 
  pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
  distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% 
  group_by(ligand) %>% filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% ungroup() 

colors_sender = brewer.pal(n = prioritized_tbl_oi$sender %>% 
                             unique() %>% sort() %>% length(), 
                           name = 'Spectral') %>% 
  magrittr::set_names(prioritized_tbl_oi$sender %>% 
                        unique() %>% sort())
colors_receiver = c("lavender")  %>% 
  magrittr::set_names(prioritized_tbl_oi$receiver %>% 
                        unique() %>% sort())

circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)
pdf("FIGURA_10D-circos_plot.pdf")
circos_output
dev.off()
dev.off()
######################
######################
#### VISUALIZACIÓN DEL SEGUNDO NICHO: Nicho pobre en EMT
#### Se puede obtener redefiniendo la variable receiver_oi
receiver_oi = "EMT_Low"  
# y aplicando el mismo código

filtered_ligands = ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(50, prioritization_score) %>% 
  pull(ligand) %>% unique()

prioritized_tbl_oi = 
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% ungroup() 

# 95 INTERACTIONS LOW
write.table(prioritized_tbl_oi, "prioritized_tbl_oi-EMT_low.tsv", quote=FALSE, row.names = FALSE, sep = "\t")

## desde lfc_plot de antes, hasta el final, 
## podemos copiar y pegar código y generar las mismas figuras
## para mi nicho pobre en EMT
## 
## 
#########################
#########################
####### Para llevar a cabo una corroboración de los resultados obtenidos,
####### se llevo a cabo el mismo análisis con los singletes en lugar de los
####### PICs. Para ello, solo hay que editar la variable fake_niches creada en
####### el step 1, y realizar el mismo análisis:
###
## Definimos los nichos
## setwd("./singlets")
fake_niches <- c()
for (cluster in assay_seurat@meta.data$group){
  if (cluster == "PROLIF_singlet" || cluster == "STRESS_singlet"){
    fake_niches <- c(fake_niches, "EMT_Low")} # nicho "pobre" en EMT
  else{
    if (cluster == "EMT_singlet"){
      fake_niches <- c(fake_niches, "EMT_High")} # nicho "rico" en EMT
    else{
      if(cluster == "AMs"){
        fake_niches <- c(fake_niches, "AMs")} # AMs 
      else{
        fake_niches <- c(fake_niches, "none")} # células que no me interesan
    }
  }
}
assay_seurat@meta.data$fake_niches <- fake_niches
##
## 92 singlet interactions
