library(config)
library(uwot)

# ---- vars ----

n_components <- 2
pca <- 50
metric <- "euclidean"

# ---- pbmc ----

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$pbmc)

sets <- ls()[grepl("cells.", ls())]

var2save <- c()

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  UMAP_results <- uwot::umap(
    X = data, 
    n_components = n_components,
    pca = pca,
    metric = metric
    )
  
  new_var <- paste0(set,".UMAP")
  
  assign(new_var, UMAP_results)
  var2save <- c(var2save, new_var)
}

save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$pbmc$UMAP)
rm(list = ls())

# ---- breast cancer ----

n_components <- 2
pca <- 50
metric <- "euclidean"

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$breast)

sets <- ls()[grepl("cells.", ls())]

var2save <- c()

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  UMAP_results <- uwot::umap(
    X = data, 
    n_components = n_components,
    pca = pca,
    metric = metric
  )
  
  new_var <- paste0(set,".UMAP")
  
  assign(new_var, UMAP_results)
  var2save <- c(var2save, new_var)
}

save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$breast$UMAP)
rm(list = ls())

# ---- liver ----

n_components <- 2
pca <- 50
metric <- "euclidean"

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$liver)

sets <- ls()[grepl("cells.", ls())]

var2save <- c()

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  UMAP_results <- uwot::umap(
    X = data, 
    n_components = n_components,
    pca = pca,
    metric = metric
  )
  
  new_var <- paste0(set,".UMAP")
  
  assign(new_var, UMAP_results)
  var2save <- c(var2save, new_var)
}

save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$liver$UMAP)
rm(list = ls())

