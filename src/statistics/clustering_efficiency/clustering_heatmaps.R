library(config)
library(dplyr)
library(reshape2)
library(cluster)
library(MetBrewer)
library(patchwork)

config <- config::get()

load(config$results$clustering$pbmc)
load(config$results$clustering$breast$disease)
load(config$results$clustering$breast$group)
load(config$results$clustering$liver)

adjust_data <- function(data){
  output_data <- data[!data$norm == "scaled", ]
  levels <- c("nonorm", "log2", "sqrt", "scran", "scnorm", "sctransform", "dino")
  output_data$norm <- factor(output_data$norm, levels=levels)
  output_data[output_data$reduction == "tsne", "reduction"] <- "tSNE"
  return(output_data)
}

# heatmaps

geom_hetmap <- function(data, legend_title){
  
  data[data$reduction == "tsne", "reduction"] = "tSNE"
  
  labs <- labs(x="", y="", fill = legend_title)
  p <- ggplot(data)+
    geom_tile(aes(x=reduction, y=norm, fill=scores))+
    geom_text(aes(x = reduction, y = norm, label = round(scores, 2)), size = 3, color="white", fontface = "bold") +
    facet_wrap(~clustering, ncol = 3) + labs + scale_fill_gradientn(colors = rev(met.brewer("Degas")))
  return(p)
}

theme <- theme(
  plot.background = element_rect(fill = "transparent", color = NA),
  strip.background = element_blank(),
  # strip.text = element_blank(),
  legend.background = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_blank(),            
  legend.key = element_rect(fill = "transparent", color = NA),   
  legend.title = element_text(hjust=0.25),
  # axis.text = element_blank(),                      
  panel.grid = element_blank(),
  panel.background = element_rect(fill = "transparent", color = NA),
  strip.text = element_text(face = "bold", size = 12),
  axis.text.x = element_text(face = "bold", size = 10),
  axis.text.y = element_text(face = "bold", size = 10)
)

# breast_disease_clustering_meas

data <- adjust_data(breast_disease_clustering_meas)

p1 <- geom_hetmap(data[data$types == "adjrand", ], "ARI") + theme
p2 <- geom_hetmap(data[data$types == "sdc", ], "Dice") + theme
p3 <- geom_hetmap(data[data$types == "mi", ], "MI") + theme
p4 <- geom_hetmap(data[data$types == "si", ], "SI") + theme

p_merged_breast_d <- (p1 + p2) / (p3 + p4)
p_merged_breast_d
# ggsave("clustering_hm_breast_d.png", p_merged_breast_d)

# breast_groups_clustering_meas

data <- adjust_data(breast_groups_clustering_meas)

p1 <- geom_hetmap(data[data$types == "adjrand", ], "ARI") + theme
p2 <- geom_hetmap(data[data$types == "sdc", ], "Dice") + theme
p3 <- geom_hetmap(data[data$types == "mi", ], "MI") + theme
p4 <- geom_hetmap(data[data$types == "si", ], "SI") + theme

p_merged_breast_g <- (p1 + p2) / (p3 + p4)
p_merged_breast_g
# ggsave("clustering_hm_breast_g.png", p_merged_breast_g)

# pbmc

data <- adjust_data(pbmc_clustering_meas)

p1 <- geom_hetmap(data[data$types == "adjrand", ], "ARI") + theme
p2 <- geom_hetmap(data[data$types == "sdc", ], "Dice") + theme
p3 <- geom_hetmap(data[data$types == "mi", ], "MI") + theme
p4 <- geom_hetmap(data[data$types == "si", ], "SI") + theme

p_merged_pbmc <- (p1 + p2) / (p3 + p4)
p_merged_pbmc

# ggsave("clustering_hm_pbmc.png", p_merged_pbmc)
# liver

data <- adjust_data(liver_clustering_meas)

p1 <- geom_hetmap(data[data$types == "adjrand", ], "ARI") + theme
p2 <- geom_hetmap(data[data$types == "sdc", ], "Dice") + theme
p3 <- geom_hetmap(data[data$types == "mi", ], "MI") + theme
p4 <- geom_hetmap(data[data$types == "si", ], "SI") + theme

p_merged_liver <- (p1 + p2) / (p3 + p4)
p_merged_liver

# ggsave("clustering_hm_liver.png", p_merged_liver, width = 3.25, height = 4, units = "in", dpi = 300)
