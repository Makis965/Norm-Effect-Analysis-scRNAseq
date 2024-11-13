library(config)
library(ggplot2)
library(ggpubr)
library(config)
library(dplyr)

config <- config::get()

source(config$utils$statistics)
source(config$utils$measures)
source(config$utils$data_adjustments)

# ---- load data ----

load(config$results$reduction$pbmc)
load(config$results$reduction$breast$disease)
load(config$results$reduction$breast$group)
load(config$results$reduction$liver)


# pbmc dataset 
pbmc_disp_vals <- adjust_data(pbmc_disp_vals)

tsne_pbmc <- pbmc_disp_vals[pbmc_disp_vals$reduction == "tSNE", ]
umap_pbmc <- pbmc_disp_vals[pbmc_disp_vals$reduction == "UMAP", ]

kurskal.wallis <- kruskal.test(tsne_pbmc$si, tsne_pbmc$norm)
H <- kurskal.wallis$statistic
pbmc_tsne_conover <- conover_test(tsne_pbmc$si, tsne_pbmc$norm, H)

kurskal.wallis <- kruskal.test(umap_pbmc$si, umap_pbmc$norm)
H <- kurskal.wallis$statistic
pbmc_umap_conover <- conover_test(umap_pbmc$si, umap_pbmc$norm, H)


# breast disease dataset 
breast_disease_disp_vals <- adjust_data(breast_disease_disp_vals)

tsne_breast_disease <- breast_disease_disp_vals[breast_disease_disp_vals$reduction == "tSNE", ]
umap_breast_disease <- breast_disease_disp_vals[breast_disease_disp_vals$reduction == "UMAP", ]

kurskal.wallis <- kruskal.test(tsne_breast_disease$si, tsne_breast_disease$norm)
H <- kurskal.wallis$statistic
breast_disease_tsne_conover <- conover_test(tsne_breast_disease$si, tsne_breast_disease$norm, H)

kurskal.wallis <- kruskal.test(umap_breast_disease$si, umap_breast_disease$norm)
H <- kurskal.wallis$statistic
breast_disease_umap_conover <- conover_test(umap_breast_disease$si, umap_breast_disease$norm, H)


# breast disease dataset 
breast_group_disp_vals <- adjust_data(breast_group_disp_vals)

tsne_breast_group <- breast_group_disp_vals[breast_group_disp_vals$reduction == "tSNE", ]
umap_breast_group <- breast_group_disp_vals[breast_group_disp_vals$reduction == "UMAP", ]

kurskal.wallis <- kruskal.test(tsne_breast_group$si, tsne_breast_group$norm)
H <- kurskal.wallis$statistic
breast_group_tsne_conover <- conover_test(tsne_breast_group$si, tsne_breast_group$norm, H)

kurskal.wallis <- kruskal.test(umap_breast_group$si, umap_breast_group$norm)
H <- kurskal.wallis$statistic
breast_group_umap_conover <- conover_test(umap_breast_group$si, umap_breast_group$norm, H)


# liver dataset 
liver_disp_vals <- adjust_data(liver_disp_vals)

tsne_liver <- liver_disp_vals[liver_disp_vals$reduction == "tSNE", ]
umap_liver <- liver_disp_vals[liver_disp_vals$reduction == "UMAP", ]

kurskal.wallis <- kruskal.test(tsne_liver$si, tsne_liver$norm)
H <- kurskal.wallis$statistic
liver_tsne_conover <- conover_test(tsne_liver$si, tsne_liver$norm, H)

kurskal.wallis <- kruskal.test(umap_liver$si, umap_liver$norm)
H <- kurskal.wallis$statistic
liver_umap_conover <- conover_test(umap_liver$si, umap_liver$norm, H)


matrix(pbmc_umap_conover, nrow=21, ncol=21) -> macierz

heatmap(liver_umap_conover)
