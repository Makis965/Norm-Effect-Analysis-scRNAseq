setwd("E:/Normalization_effects")

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

sets <- setdiff(ls(), omit_vars)

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

save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$pbmc$tSNE)
rm(list = ls())

# ---- breast cancer ----

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$breast)

sets <- setdiff(ls(), omit_vars)

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

save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$breast$tSNE)
rm(list = ls())

# ---- liver ----

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$liver)

sets <- setdiff(ls(), omit_vars)

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

save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$liver$tSNE)
rm(list = ls())
