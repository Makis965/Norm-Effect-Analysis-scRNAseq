library(config)
library(ggplot2)
library(ggpubr)
library(config)
library(dplyr)

config <- config::get()

source(config$utils$statistics)
source(config$utils$measures)
source(config$utils$data_adjustments)

# clusters dispersion

# ---- load data ----

load(config$results$reduction$pbmc)
load(config$results$reduction$breast$disease)
load(config$results$reduction$breast$group)
load(config$results$reduction$liver)

perform_tests <- function(data) {
  # Run Kruskal-Wallis and Conover test for tSNE
  tsne_data <- data[data$reduction == "tSNE", ]
  kruskal_tsne <- kruskal.test(tsne_data$si, tsne_data$norm)
  tsne_conover <- conover_test(tsne_data$si, tsne_data$norm, kruskal_tsne$statistic)
  
  # Run Kruskal-Wallis and Conover test for UMAP
  umap_data <- data[data$reduction == "UMAP", ]
  kruskal_umap <- kruskal.test(umap_data$si, umap_data$norm)
  umap_conover <- conover_test(umap_data$si, umap_data$norm, kruskal_umap$statistic)
  
  # Return the results as a list
  list(
    tsne_conover = tsne_conover,
    kruskal_tsne = kruskal_tsne,
    umap_conover = umap_conover,
    kruskal_umap = kruskal_umap
  )
}

# pbmc dataset 
pbmc_disp_vals <- adjust_data(pbmc_disp_vals)
pbmc_results <- perform_tests(pbmc_disp_vals)

# breast disease dataset 
breast_disease_disp_vals <- adjust_data(breast_disease_disp_vals)
breast_disease_results <- perform_tests(breast_disease_disp_vals)

# breast disease dataset 
breast_group_disp_vals <- adjust_data(breast_group_disp_vals)
breast_group_results <- perform_tests(breast_group_disp_vals)

# liver dataset 
liver_disp_vals <- adjust_data(liver_disp_vals)
liver_results <- perform_tests(liver_disp_vals)

stat_results <- list(
  pbmc = pbmc_results,
  breast.disease = breast_disease_results,
  breast.group = breast_group_results,
  liver = liver_results
)

save(stat_results, file = config$results$statistics$conover)
