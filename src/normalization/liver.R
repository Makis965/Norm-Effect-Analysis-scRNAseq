rm(list = ls())

library(config)

config <- config::get()

library(Dino)
library(scran)
library(SCnorm)
library(Seurat)
library(stats)
library(dynamicTreeCut)

source(config$utils$normalization)

# ---- liver ----

load(config$data$raw_data$liver)

# 0. No normalization. 
cells.nonorm <- top_variance(expression_df = cells, gene_names = genes)

# scaling for simple normalizations
cells.scaled<-data.frame(sums=colSums(cells))
med<-median(cells.scaled$sums)
cells.scaled<-t(t(cells*med)/cells.scaled$sums)

# 1. log2 
log2.time <- system.time(cells.log2<-log2(cells.scaled+1))[3]
cells.log2<-top_variance(expression_df = cells.log2, gene_names = genes)

# 2. square root 
ft.time <- system.time(cells.sqrt<-sqrt(cells.scaled)+sqrt(cells.scaled+1))[3]
cells.sqrt<-top_variance(expression_df = cells.sqrt, gene_names = genes)

# ---- scRNA-seq intendeed normalization ----

# 3. dino 
dino.time <- system.time(cells.dino <- Dino(as.matrix(cells), nCores=0))[3]
cells.dino <- top_variance(expression_df = cells.dino, gene_names = genes)

# 4. scran 
scran.time <- system.time(cells.scran <- scran_norm(expression_df = cells, metadata = meta))[3]
cells.scran <- top_variance(expression_df = cells.scran, gene_names = genes)

# 5. SCtransform  
sctransform.time <- system.time(cells.sctransform <- sctransform_norm(expression_df = cells))[3]
gene_names <- rownames(cells.sctransform@assays[["SCT"]]@counts)
cells.sctransform <- top_variance(
  expression_df = cells.sctransform@assays$SCT$counts, 
  gene_names = gene_names,
  genes_num = nrow(cells)
  )

# 6. SCnorm  

# following the documentation, since the Conditions parameter is required, 
# to make this 100% unsupervised all the samples were categorized to single 
# group...

Conditions <- scnorm_conditions(as.matrix(cells))

# cells.scnorm <- SCnorm(cells, Conditions = rep(1, each = ncol(cells)))
scnorm.time <- system.time(cells.scnorm <- SCnorm(cells, Conditions = Conditions, ditherCounts = TRUE))[3]
cells.scnorm <- cells.scnorm@assays@data@listData$normcounts
cells.scnorm <- top_variance(expression_df = cells.scnorm, gene_names = genes)

# save workspace to RData file...

times <- c(
  "log2"=log2.time,
  "ft"=ft.time,
  "dino"=dino.time,
  "scran"=scran.time,
  "sctransform"=sctransform.time,
  "scnorm"=scnorm.time
)

write.csv(as.data.frame(times), "liver_norm_times.csv")

# save(list=ls()[grepl("cells.", ls())], file = config$data$normalized$liver)
