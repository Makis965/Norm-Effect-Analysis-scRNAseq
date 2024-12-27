#get only nonorm comparisons

get_nonorm_comp <- function(data){
  data <- data[data$Technique1 == "nonorm" | data$Technique2 == "nonorm", ]
  data <- data[c("pval", "p.adj", "conover.d")]
  return(data)
}

get_nonorm_comp(breast.group[breast.group$reduction == "tSNE", ])
stat_results$breast.group$kruskal_tsne

get_nonorm_comp(breast.group[breast.group$reduction == "UMAP", ])
stat_results$breast.group$kruskal_umap

# pbmc 

get_nonorm_comp(pbmc[pbmc$reduction == "tSNE", ])
stat_results$pbmc$kruskal_tsne

get_nonorm_comp(pbmc[pbmc$reduction == "UMAP", ])
stat_results$pbmc$kruskal_umap

# liver 

get_nonorm_comp(liver[liver$reduction == "tSNE", ])
stat_results$liver$kruskal_tsne

get_nonorm_comp(liver[liver$reduction == "UMAP", ])
stat_results$liver$kruskal_umap