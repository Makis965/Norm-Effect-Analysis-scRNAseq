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

# 5. SCtransform  
cells.sctransform <- sctransform_norm(expression_df = cells)
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
cells.scnorm <- SCnorm(cells, Conditions = Conditions, ditherCounts = TRUE)
cells.scnorm <- cells.scnorm@assays@data@listData$normcounts
cells.scnorm <- top_variance(expression_df = cells.scnorm, gene_names = genes)

# save workspace to RData file...

save(list=ls()[grepl("cells.", ls())], file = config$data$normalized$liver)
