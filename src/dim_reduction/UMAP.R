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
sets <- sets[-4]
pbmc.times <- c()
var2save <- c()

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  pbmc.times <- c(
    pbmc.times, 
    system.time(
      UMAP_results <- uwot::umap(
        X = data, 
        n_components = n_components,
        pca = pca,
        metric = metric
        )
      )[3]
    )
  
  new_var <- paste0(set,".UMAP")
  
  assign(new_var, UMAP_results)
  var2save <- c(var2save, new_var)
}

pbmc.times <- c(pbmc.times, mean(pbmc.times))
write.csv(as.data.frame(pbmc.times), "UMAP_pbmc_runtime.csv")
# save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$pbmc$UMAP)
rm(list = ls())

# ---- breast cancer ----

n_components <- 2
pca <- 50
metric <- "euclidean"

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$breast)

sets <- ls()[grepl("cells.", ls())]
breast.times <- c()

var2save <- c()

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  breast.times <- c(
    breast.times,
    system.time(
      UMAP_results <- uwot::umap(
      X = data, 
      n_components = n_components,
      pca = pca,
      metric = metric
      )
    )[3]
  )
  
  new_var <- paste0(set,".UMAP")
  
  assign(new_var, UMAP_results)
  var2save <- c(var2save, new_var)
}

breast.times <- c(breast.times, mean(breast.times))
write.csv(as.data.frame(breast.times), "UMAP_breast_runtime.csv")
# save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$breast$UMAP)
rm(list = ls())

# ---- liver ----

n_components <- 2
pca <- 50
metric <- "euclidean"

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$liver)

sets <- ls()[grepl("cells.", ls())]
sets <- sets[-4]
liver.times <- c()
var2save <- c()

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  liver.times <- c(
    liver.times,
    system.time(
      UMAP_results <- uwot::umap(
        X = data, 
        n_components = n_components,
        pca = pca,
        metric = metric
      )
    )[3]
  )
  
  new_var <- paste0(set,".UMAP")
  
  assign(new_var, UMAP_results)
  var2save <- c(var2save, new_var)
}

liver.times <- c(liver.times, mean(liver.times))
write.csv(as.data.frame(liver.times), "UMAP_liver_runtime.csv")
# save(list=ls()[grepl(".UMAP", ls())], file = config$data$reduced$liver$UMAP)

