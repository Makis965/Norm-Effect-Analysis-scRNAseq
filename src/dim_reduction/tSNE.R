library(config)
library(Rtsne)

# ---- pbmc ----

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$pbmc)

sets <- setdiff(ls(), "config")

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  tsne_results <- Rtsne(data)
  assign(paste0(set,'.tsne'), tsne_results)
}

save(list=ls()[grepl(".tsne", ls())], file = config$data$reduced$pbmc$tSNE)
rm(list = ls())

# ---- breast cancer ----

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$breast)

sets <- setdiff(ls(), "config")

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  tsne_results <- Rtsne(data)
  assign(paste0(set,'.tsne'), tsne_results)
}

save(list=ls()[grepl(".tsne", ls())], file = config$data$reduced$breast$tSNE)
rm(list = ls())

# ---- liver ----

config <- config::get()
source(config$utils$dim_reduction)

load("normalization/normalized_liver.RData")

sets <- setdiff(ls(), "config")

for(set in sets){
  
  data <- as.matrix(t(base::get(set)))
  
  tsne_results <- Rtsne(data)
  assign(paste0(set,'.tsne'), tsne_results)
}

save(list=ls()[grepl(".tsne", ls())], file = config$data$reduced$liver$tSNE)
rm(list = ls())
