library(config)
library(ggplot2)
library(ggpubr)
library(config)
library(dplyr)

config <- config::get()

source(config$utils$data_adjustments)

# load data
load(config$results$reduction$pbmc)
load(config$results$reduction$breast$disease)
load(config$results$reduction$breast$group)
load(config$results$reduction$liver)

# initial data adjustment
pbmc_disp_vals <- adjust_data(pbmc_disp_vals)
breast_disease_disp_vals <- adjust_data(breast_disease_disp_vals)
breast_group_disp_vals <- adjust_data(breast_group_disp_vals)
liver_disp_vals <- adjust_data(liver_disp_vals)


combinations <- combn(c("nonorm", "log2", "sqrt", "scran", "scnorm", "sctransform", "dino"), 2)
combinations <- as.data.frame(matrix(combinations, ncol=2))
norm_comparisons <- combinations[combinations$V1 == "nonorm", ]

norm_comparisons <- apply(norm_comparisons, 1, function(x) list(as.vector(x)))

  
#calculate the means for each CellType
pbmc_disp_means <- pbmc_disp_vals %>%
  group_by(CellType, norm) %>%
  summarize(SI = mean(si))

pbmc_disp_means <- pbmc_disp_vals %>%
  group_by(CellType, norm) %>%
  summarize(SI = mean(si))

pbmc_disp_means <- pbmc_disp_vals %>%
  group_by(CellType, norm) %>%
  summarize(SI = mean(si))

liver_disp_means <- liver_disp_vals %>%
  group_by(Cell.type.org, norm) %>%
  summarize(SI = mean(si))

ggplot(pbmc_disp_vals)+
  geom_violin(aes(x=norm, y=si, fill=norm), alpha=0.5)+
  geom_point(data=pbmc_disp_means, aes(x=norm, y=SI, color=CellType), alpha=1)+
  facet_wrap(~reduction)+
  ggpubr::stat_compare_means(aes(x=norm, y=si), method = "kruskal.test", vjust=-1, hjust=-1)

  # ggpubr::stat_compare_means(method = "tukey")
  

ggplot(liver_disp_vals)+
  geom_violin(aes(x=norm, y=si, fill=norm), alpha=0.5)+
  geom_point(data=liver_disp_means, aes(x=norm, y=SI, color=Cell.type.org), alpha=1)+
  facet_wrap(~reduction)
  # ggpubr::stat_compare_means(aes(x=norm, y=si), method = "kruskal.test")
