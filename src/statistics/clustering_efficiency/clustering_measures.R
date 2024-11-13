library(config)
library(dplyr)
library(reshape2)
library(cluster)

config <- config::get()

# ---- evaluate clustering efficiency ----

source(config$utils$measures)

#pbmc 
dataset = "pbmc"
cell_labels = config$labels$pbmc
cell_label_names = "CellType"
pbmc_clustering_meas <- apply_stat_comp(dataset, cell_labels, cell_label_names)

save(pbmc_clustering_meas, file = config$results$clustering$pbmc)

#breast disease
dataset = "breast"
cell_labels = config$labels$breast$disease
cell_label_names = "Disease"
breast_disease_clustering_meas <- apply_stat_comp(dataset, cell_labels, cell_label_names)

save(breast_disease_clustering_meas, file = config$results$clustering$breast$disease)

#breast cell-types
dataset = "breast"
cell_labels = config$labels$breast$disease
cell_label_names = "group"
breast_groups_clustering_meas <- apply_stat_comp(dataset, cell_labels, cell_label_names)

save(breast_groups_clustering_meas, file = config$results$clustering$breast$group)

#liver
dataset = "liver"
cell_labels = config$labels$liver
cell_label_names = "Cell.type.org"
liver_clustering_meas <- apply_stat_comp(dataset, cell_labels, cell_label_names)

save(liver_clustering_meas, file = config$results$clustering$liver)