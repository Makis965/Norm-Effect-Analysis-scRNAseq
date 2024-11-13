library(Seurat)
library(SeuratObject)
library(igraph)
library(cluster)
library(config)

config <- config::get()

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
  cells.df$louvain <- optimal_k_louvain(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".louvain")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$pbmc$UMAP$louvain)

rm(list = ls())


# ---- breast cancer ----

config <- config::get()

source(config$utils$clustering)

load(config$data$reduced$breast$UMAP)

sets <- ls()[grepl("cells", ls())]

load(config$data$raw_data$breast)

var2save <- c()
for(set in sets){
  
  print(paste("set:",set))
  
  data <- base::get(set)
  
  cells.df <- as.data.frame(data)
  cells.df$louvain <- optimal_k_louvain(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".louvain")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$breast$UMAP$louvain)

rm(list = ls())


# ---- liver ----

config <- config::get()

source(config$utils$clustering)

load(config$data$reduced$liver$UMAP)

sets <- ls()[grepl("cells", ls())]

load(config$data$raw_data$liver)

var2save <- c()
for(set in sets){
  
  print(paste("set:",set))
  
  data <- base::get(set)
  
  cells.df <- as.data.frame(data)
  cells.df$louvain <- optimal_k_louvain(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".louvain")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$liver$UMAP$louvain)

rm(list = ls())
