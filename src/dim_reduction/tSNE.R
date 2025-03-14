library(config)
library(Rtsne)

# ---- pbmc ----

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$pbmc)

sets <- ls()[grepl("cells.", ls())]
pbmc.times <- c()
for(set in sets){
  
  data <- t(as.matrix(base::get(set)))
  
  pbmc.times <- c(pbmc.times, system.time(tsne_results <- Rtsne(data))[3])
  assign(paste0(set,'.tsne'), tsne_results)
}

pbmc.times <- c(pbmc.times, mean(pbmc.times))
write.csv(as.data.frame(pbmc.times), "tSNE_pbmc_runtime.csv")
# save(list=ls()[grepl(".tsne", ls())], file = config$data$reduced$pbmc$tSNE)
rm(list = ls())

# ---- breast cancer ----

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$breast)

sets <- ls()[grepl("cells.", ls())]
breast.times <- c()
for(set in sets){
  
  data <- t(as.matrix(base::get(set)))
  
  breast.times <- c(breast.times, system.time(tsne_results <- Rtsne(data))[3])
  assign(paste0(set,'.tsne'), tsne_results)
}

breast.times <- c(breast.times, mean(breast.times))
write.csv(as.data.frame(breast.times), "tSNE_breast_runtime.csv")
# save(list=ls()[grepl(".tsne", ls())], file = config$data$reduced$breast$tSNE)
rm(list = ls())

# ---- liver ----

config <- config::get()
source(config$utils$dim_reduction)

load(config$data$normalized$liver)

sets <- ls()[grepl("cells.", ls())]
liver.times <- c()
for(set in sets){
  
  data <- t(as.matrix(base::get(set)))
  
  liver.times <- c(liver.times, system.time(tsne_results <- Rtsne(data))[3])
  assign(paste0(set,'.tsne'), tsne_results)
}

liver.times <- c(liver.times, mean(liver.times))
write.csv(as.data.frame(liver.times), "tSNE_liver_runtime.csv")
# save(list=ls()[grepl(".tsne", ls())], file = config$data$reduced$liver$tSNE)
rm(list = ls())

