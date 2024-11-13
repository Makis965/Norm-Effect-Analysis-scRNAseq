library(cluster)
library(config)

config <- config::get()

source(config$utils$measures)

#pbmc
dataset = "pbmc"
cell_labels = config$labels$pbmc
cell_label_names = "CellType"

pbmc_disp_vals <- compute_dispersions(dataset, cell_labels, cell_label_names)
save(pbmc_disp_vals, file = config$results$reduction$pbmc)


#breast disease
dataset = "breast"
cell_labels = config$labels$breast$disease
cell_label_names = "Disease"

breast_disease_disp_vals <- compute_dispersions(dataset, cell_labels, cell_label_names)
save(breast_disease_disp_vals, file = config$results$reduction$breast$disease)

#breast disease
dataset = "breast"
cell_labels = config$labels$breast$group
cell_label_names = "group"

breast_group_disp_vals <- compute_dispersions(dataset, cell_labels, cell_label_names)
save(breast_group_disp_vals, file = config$results$reduction$breast$group)

#liver
dataset = "liver"
cell_labels = config$labels$liver
cell_label_names = "Cell.type.org"

liver_disp_vals <- compute_dispersions(dataset, cell_labels, cell_label_names)
save(liver_disp_vals, file = config$results$reduction$liver)
