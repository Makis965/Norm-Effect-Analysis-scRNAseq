setwd("E:/Normalization_effects")

library(config)

config <- config::get()

library(Dino)
library(scran)
library(SCnorm)
library(Seurat)
library(stats)

source(config$utils$normalization)

# ---- pbmc ----

load(config$data$raw_data$pbmc)

# 0. No normalization. 

cells.nonorm <- top_variance(expression_df = cells, gene_names = genes)

# scaling for simple normalizations

cells.scaled<-data.frame(sums=colSums(cells))
med<-median(cells.scaled$sums)
cells.scaled<-t(t(cells*med)/cells.scaled$sums)

# 1. log2 
cells.log2<-log2(cells.scaled+1)
cells.log2<-top_variance(expression_df = cells.log2, gene_names = genes)

# 2. square root 
cells.sqrt<-sqrt(cells.scaled)+sqrt(cells.scaled+1)
cells.sqrt<-top_variance(expression_df = cells.sqrt, gene_names = genes)

# ---- scRNA-seq intendeed normalization ----

# 3. dino 
cells.dino <- Dino(as.matrix(cells))
cells.dino <- top_variance(expression_df = cells.dino, gene_names = genes)

# 4. scran 
cells.scran <- scran_norm(expression_df = cells, metadata = meta)
cells.scran <- top_variance(expression_df = cells.scran, gene_names = genes)

# 5. seurat  
cells.seurat <- seurat_norm(expression_df = cells)
cells.seurat <- top_variance(expression_df = cells.seurat, gene_names = genes)

# 6. SCnorm  

# following the documentation, since the Conditions parameter is required, 
# to make this 100% unsupervised all the samples were categorized to single 
# group...

# cells.scnrom <- SCnorm(cells, Conditions = as.vector(meta$CellType))
cells.scnorm <- SCnorm(cells, Conditions = c(1:ncol(cells)))
cells.scnorm <- cells.scnorm@assays@data@listData$normcounts
cells.scnorm <- top_variance(expression_df = cells.scnorm, gene_names = genes)


# save workspace to RData file...

save(list=ls()[grepl("cells.", ls())], file = config$data$normalized$pbmc)
