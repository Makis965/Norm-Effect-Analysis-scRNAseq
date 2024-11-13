library(config)

config <- config::get()

top_var_genes <- 0.2

top_variance <- function(expression_df, gene_names, genes_num = NaN){
  #' function filters genes and keeps those with the highest variance
  #' between samples.
  
  top_var <- data.frame(
    expression_var=rowVars(
      as.matrix(expression_df), 
      useNames=TRUE
    )
  )
  
  top_var$genes <- gene_names
  
  top_var <- expression_df[
    sort.list(
      top_var$expression_var, 
      decreasing=TRUE
    ), 
  ]
  
  if(is.nan(genes_num) == FALSE){
    top_count <- genes_num * top_var_genes
  } else {
    top_count <- nrow(expression_df) * top_var_genes
  }
  
  return(top_var[1:top_count, ])
}


clear_workspace <- function() {
  rm(
    meta, 
    cells, 
    med, 
    data, 
    config, 
    genes, 
    cells.scaled,
    scale_expression,
    scran_norm,
    seurat_norm,
    top_variance,
    clear_workspace,
    envir = .GlobalEnv
    )
}

# ---- data normalization and transformation ---- 

#' All the functions takes the cell expression matrices with
#' samples on columns and gene expression level on rows...

scale_expression <- function(expression_df){
  #' Function proceed median scaling as a pre-processing 
  #' step for simple normalization.
}

scran_norm <- function(expression_df, metadata){
  #' Function proceed scran normalization
  
  sce <- SingleCellExperiment(
    assays=list(
      counts=expression_df,
      logcounts=log2(expression_df+1)
    ),
    metadata=metadata
  )
  
  clusters <- quickCluster(sce)
  
  sce <- computeSumFactors(sce, clusters=clusters)
  sce <- logNormCounts(sce)
  
  return(sce@assays@data@listData[["logcounts"]])
}

sctransform_norm <- function(expression_df){
  options(future.globals.maxSize = 8000 * 1024^2)
  seurat.object <- Seurat::CreateSeuratObject(expression_df)
  sctransform_norm <- Seurat::SCTransform(seurat.object, vst.flavor="v2")
  
  return(sctransform_norm)
}

seurat_norm <- function(expression_df){
  #' Function proceed Seurat normalization
  
  seurat.object <- CreateSeuratObject(expression_df)
  seurat.norm <- NormalizeData(
    seurat.object, 
    normalization.method = "LogNormalize"
  )
  
  return(GetAssayData(seurat.norm))
}

scnorm_conditions <- function(expression_df){
  dist <- as.dist(1 - cor(expression_df, method="spearman"))
  htree <- hclust(dist, method="ward.D")
  
  clusters <- factor(
    unname(
      cutreeDynamic(
        htree,
        minClusterSize=50,
        method="tree",
        respectSmallClusters=FALSE
      )
    )
  )
  
  names(clusters) <- colnames(expression_df)
  
  return(clusters)
}
