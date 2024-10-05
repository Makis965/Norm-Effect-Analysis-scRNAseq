library(cluster)
library(config)
library(aricode)

config <- config::get()

source(config$utils$statistics)

# ---- pbmc ---- 


# tSNE

load(config$data$reduced$pbmc$tSNE)
load(config$data$raw_data$pbmc)

cell_labels <- meta$CellType
decode_labels <- unlist(config$labels$pbmc)[cell_labels]

sets <- ls()[grepl(".cells", ls())]

pbmc.tSNE.stats <- process_clustered_data(sets, decode_labels, "cell_type", cell_labels)


# UMAP

load(config$data$reduced$pbmc$UMAP)
load(config$data$raw_data$pbmc)

sets <- ls()[grepl(".cells", ls())]

pbmc.UMAP.stats <- process_clustered_data(sets, decode_labels, "cell_type", cell_labels)

rm(list = ls()[grepl("cells", ls())])


# ---- breast cancer ---- 


# tSNE 

load(config$data$reduced$breast$tSNE)
load(config$data$raw_data$breast)

sets <- ls()[grepl(".cells", ls())]

# for disease occurences:

disease_labels <- meta$Disease
decode_labels <- unlist(config$labels$breast$disease)[disease_labels]

breast.tSNE.disease <- process_clustered_data(sets, decode_labels, "disease", disease_labels)

# for disease types:

disease_labels <- meta$group
decode_labels <- unlist(config$labels$breast$group)[disease_labels]

breast.tSNE.group <- process_clustered_data(sets, decode_labels, "group", disease_labels)


# UMAP 

load(config$data$reduced$breast$UMAP)
load(config$data$raw_data$breast)

sets <- ls()[grepl(".cells", ls())]

# for disease occurences:

disease_labels <- meta$Disease
decode_labels <- unlist(config$labels$breast$disease)[disease_labels]

breast.UMAP.stats.disease <- process_clustered_data(sets, decode_labels, "disease", disease_labels)

# for disease types:

disease_labels <- meta$group
decode_labels <- unlist(config$labels$breast$group)[disease_labels]

breast.UMAP.stats.group <- process_clustered_data(sets, decode_labels, "group", disease_labels)

rm(list = ls()[grepl("cells", ls())])


# ---- liver ----- 


# tSNE

load(config$data$reduced$liver$tSNE)
load(config$data$raw_data$liver)
cell_labels <- meta$Cell.type.org

decode_labels <- unlist(config$labels$liver)[cell_labels]

sets <- ls()[grepl(".cells", ls())]

liver.tSNE.stats <- process_clustered_data(sets, decode_labels, "cell_type", cell_labels)


# UMAP

load(config$data$reduced$liver$UMAP)

sets <- ls()[grepl(".cells", ls())]

liver.UMAP.stats <- process_clustered_data(sets, decode_labels, "cell_type", cell_labels)

rm(list = ls()[grepl("cells", ls())])


# ---- save files ----

to_save <- ls()[grepl(".stats", ls())]
dim_reduction_stats <- list()

for(file in to_save){
  dim_reduction_stats[[file]] <- base::get(file)
}


save(dim_reduction_stats, file=config$data$statistics$dim_reduction)
