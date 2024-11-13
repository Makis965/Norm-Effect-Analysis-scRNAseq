# ---- hierarchical clustering ----

optimal_k_hclust <- function(cells){
  
  si_avg <- (-1)
  
  distances <- dist(cells)
  clustering <- hclust(distances, method = "complete")
  
  for(k in 2:20){
    clusters <- cutree(clustering, k)
    si_val <- silhouette(clusters, distances)
    
    si_avg_new <- mean(si_val[, 3])
    
    if(si_avg_new > si_avg){
      si_avg <- si_avg_new
      best_clusters <- clusters
    }
  }
  
  return(best_clusters)
}

# ---- k-means ----

optimal_k_kmeans <- function(cells){
  
  si_avg <- (-1)
  
  distances <- dist(cells)
  
  for(k in 2:20){
    clusters <- kmeans(x = cells, centers = k)$cluster
    si_val <- silhouette(clusters, distances)
    
    si_avg_new <- mean(si_val[, 3])
    
    if(si_avg_new > si_avg){
      si_avg <- si_avg_new
      best_clusters <- clusters
    }
  }
  
  return(best_clusters)
  
}

# ---- Louvain community detection ----

custom_seurat_object <- function(cells){
  #' Since dimensionality reduction has been proceed within previous steps, 
  #' the SeuratObject has to be customized...
  
  seurat_obj <- CreateSeuratObject(counts = t(as.matrix(cells)))
  
  seurat_obj@reductions[["reduction"]] <- CreateDimReducObject(
    embeddings = cells, 
    key = "custom_",
    assay = DefaultAssay(seurat_obj)
  )
  
  return(seurat_obj)
}

optimal_k_louvain <- function(cells){
  seurat_object <- custom_seurat_object(cells)
  
  si_avg <- (-1)
  
  distances <- dist(cells)
  
  for(k in seq(5, 100, 5)){
    for(resolution in seq(0.4, 2, 0.1)){
      
      seurat_object <- FindNeighbors(
        seurat_object, 
        reduction = "reduction", 
        k.param = k,
        dims = 1:2)
      
      #TODO: resolution parameters might be verified as well... 
      #Previously I have been optimizing resolution parameter only (!)
      seurat_object <- FindClusters(seurat_object, resolution = resolution)
      
      clusters <- as.numeric(seurat_object@meta.data$seurat_clusters)
      
      si_val <- silhouette(clusters, distances)
      
      si_avg_new <- mean(si_val[, 3])
      if(si_avg_new > si_avg){
        si_avg <- si_avg_new
        best_clusters <- clusters
      }
    }
  }
  return(best_clusters)
}

# ---- autoencoder scCAN ----
#TODO: Works below CUDA 11.8 (current installed version is 12.3)
optimal_k_sccan <- function(cells){
  
}