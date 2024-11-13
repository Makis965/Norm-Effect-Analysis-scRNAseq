library(mclustcomp)
library(fpc)
library(cluster)
library(config)

# ---- compute dispersion measures ----

compute_dispersions <- function(dataset, cell_labels, cell_label_names){
  meta_data <- load(config$data$raw_data[[dataset]])
  tsne <- load(config$data$reduced[[dataset]]$tSNE)
  umap <- load(config$data$reduced[[dataset]]$UMAP)
  all_sets <- c(tsne, umap)
  
  merged_df <- data.frame()
  for(tsne_set in tsne){
    base::assign(tsne_set, base::get(tsne_set)$Y)
  }
  
  
  for(set in all_sets){
    data <- base::get(set)
    cell_labels <- unlist(cell_labels)[meta[[cell_label_names]]]
    
    dist <- stats::dist(data)
    si <- silhouette(cell_labels, dist)
    splited_text <- strsplit(set, "\\.")[[1]]
    
    norm = splited_text[2]
    reduction = splited_text[3]
    meta_df <- data.frame(
      dataset=dataset, 
      norm=norm, 
      reduction=reduction,
      si=si[, 3]
    )
    
    temp_df <- cbind(data, meta_df, base::get("meta"))
    merged_df <- rbind(merged_df, temp_df)
  }
  
  return(merged_df)
  
}

# ---- compute clustering measures ----

compute_stats <- function(data_list, dataset, cell_labels, cell_label_names){
  output_df = data.frame()
  
  mclust_methods = c("adjrand", "f", "mi")
  data_names <- names(data_list)
  
  for (i in 1:length(data_list)){
    splited_text = strsplit(data_names[i], "\\.")[[1]]
    data <- data_list[[i]]
    
    norm = splited_text[2]
    reduction = splited_text[3]
    clustering = splited_text[4]
    
    cell_labels = unlist(cell_labels)[data[[cell_label_names]]]
    
    dists <- stats::dist(data[, 1:2])
    si <- silhouette(data[[clustering]], dists)
    
    res <- mclustcomp(data[[clustering]], cell_labels, types = mclust_methods)
    res <- rbind(res, data.frame(types="si", scores=mean(si[, 3])))
    temp_df = cbind(data.frame(dataset=dataset, norm=norm, reduction=reduction, clustering=clustering), res)
    output_df = rbind(output_df, temp_df)
  }
  
  rownames(output_df) <- NULL
  return(output_df)
}

apply_stat_comp <- function(dataset, cell_labels, cell_label_names){
  
  source(config$utils$statistics)
  
  merged_df <- data.frame()
  
  clustering_methods <- c("hclust", "kmeans", "louvain")
  
  for(clustering_method in clustering_methods){
    tnse_clust <- load(config$data$clustered[[dataset]]$tSNE[[clustering_method]])
    umap_clust <- load(config$data$clustered[[dataset]]$UMAP[[clustering_method]]) 
    total_clust <- c(tnse_clust, umap_clust)
    
    data_list <- list()
    for(i in 1:length(total_clust)){
      data_list[[total_clust[i]]] <- base::get(total_clust[i])
    }
    
    stats_df <- compute_stats(
      data_list = data_list,
      dataset = dataset,
      cell_labels = cell_labels,
      cell_label_names = cell_label_names
    )
    
    merged_df <- rbind(merged_df, stats_df)
    
  }
  
  return(merged_df)
  
}

