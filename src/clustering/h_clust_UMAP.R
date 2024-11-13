library(dplyr)
library(cluster)
library(config)

config <- config::get()

# ---- UMAP ----

# ---- pbmc ----

source(config$utils$clustering)

load(config$data$reduced$pbmc$UMAP)

sets <- ls()[grepl("cells.", ls())]

load(config$data$raw_data$pbmc)

var2save <- c()
for(set in sets){
  
  print(paste("set:",set))
  
  data <- base::get(set)
  
  cells.df <- as.data.frame(data)
  cells.df$hclust <- optimal_k_hclust(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".hclust")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$pbmc$UMAP$hclust)
rm(list = ls())


# ---- breast cancer ----

config <- config::get()

source(config$utils$clustering)

load(config$data$reduced$breast$UMAP)

sets <- ls()[grepl("cells.", ls())]

load(config$data$raw_data$breast)

var2save <- c()
for(set in sets){
  
  print(paste("set:",set))
  
  data <- base::get(set)
  
  cells.df <- as.data.frame(data)
  cells.df$hclust <- optimal_k_hclust(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".hclust")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$breast$UMAP$hclust)

rm(list = ls())


# ---- liver ----

config <- config::get()

source(config$utils$clustering)

load(config$data$reduced$liver$UMAP)

sets <- ls()[grepl("cells.", ls())]

load(config$data$raw_data$liver)

var2save <- c()
for(set in sets){
  
  print(paste("set:",set))
  
  data <- base::get(set)
  
  cells.df <- as.data.frame(data)
  cells.df$hclust <- optimal_k_hclust(cells = data)
  cells.df <- bind_cols(cells.df, meta)
  
  new_var <- paste0(set, ".hclust")
  var2save <- c(var2save, new_var)
  
  assign(new_var, cells.df)
  
}

save(list = var2save, file = config$data$clustered$liver$UMAP$hclust)

rm(list = ls())
