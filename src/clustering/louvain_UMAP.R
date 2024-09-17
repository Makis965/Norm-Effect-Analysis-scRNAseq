setwd("E:/Normalization_effects")

library(Seurat)
library(SeuratObject)
library(igraph)
library(cluster)

# ---- UMAP ----

# ---- pbmc ----

source(config$utils$clustering)

load(config$data$reduced$pbmc$UMAP)

sets <- ls()[grepl("cells", ls())]

load(config$data$raw_data$pbmc)

var2save <- c()
for(set in sets){
  
  print(paste("set:",set))
  
  data <- base::get(set)
  
  cells.df <- as.data.frame(data)
  cells.df$kmeans <- optimal_k_louvain(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".louvain")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$pbmc$UMAP$louvain)

rm(list = ls())


# ---- breast cancer ----

source(config$utils$clustering)

load(config$data$reduced$breast$UMAP)

sets <- ls()[grepl("cells", ls())]

load(config$data$raw_data$breast)

var2save <- c()
for(set in sets){
  
  print(paste("set:",set))
  
  data <- base::get(set)
  
  cells.df <- as.data.frame(data)
  cells.df$kmeans <- optimal_k_louvain(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".louvain")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$breast$UMAP$louvain)

rm(list = ls())


# ---- liver ----

source(config$utils$clustering)

load(config$data$reduced$liver$UMAP)

sets <- ls()[grepl("cells", ls())]

load(config$data$raw_data$liver)

var2save <- c()
for(set in sets){
  
  print(paste("set:",set))
  
  data <- base::get(set)
  
  cells.df <- as.data.frame(data)
  cells.df$kmeans <- optimal_k_louvain(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".louvain")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$liver$UMAP$louvain)

rm(list = ls())
