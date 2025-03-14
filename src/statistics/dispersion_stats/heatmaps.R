library(config)
library(ggplot2)
library(ggpubr)
library(config)
library(dplyr)
library(MetBrewer)

config <- config::get()

load(config$results$statistics$conover)

# analyze statistics: 

get_nonorm_comparisons <- function(data){
  data <- as.data.frame(t(data))
  
  method_1 <- c()
  method_2 <- c()
  for(i in rownames(data)){
    methods <- unlist(strsplit(i, ":"))
    method_1 <- c(method_1, methods[1])
    method_2 <- c(method_2, methods[2])
  }
  
  data$Technique1 <- method_1
  data$Technique2 <- method_2
  
  # return(data[data$Technique1 == "nonorm" | data$Technique2 == "nonorm", ])
  return(data)
  
}

#breast disease

breast.dis.tsne <- get_nonorm_comparisons(stat_results$breast.disease$tsne_conover)
breast.dis.tsne$reduction <- "tSNE"
breast.dis.umap <- get_nonorm_comparisons(stat_results$breast.disease$umap_conover)
breast.dis.umap$reduction <- "UMAP"

breast.dis <- bind_rows(breast.dis.tsne, breast.dis.umap)

#breast group

breast.group.tsne <- get_nonorm_comparisons(stat_results$breast.group$tsne_conover)
breast.group.tsne$reduction <- "tSNE"
breast.group.umap <- get_nonorm_comparisons(stat_results$breast.group$umap_conover)
breast.group.umap$reduction <- "UMAP"

breast.group <- bind_rows(breast.group.tsne, breast.group.umap)

# pbmc 

pbmc.tsne <- get_nonorm_comparisons(stat_results$pbmc$tsne_conover)
pbmc.tsne$reduction <- "tSNE"
pbmc.umap <- get_nonorm_comparisons(stat_results$pbmc$umap_conover)
pbmc.umap$reduction <- "UMAP"

pbmc <- bind_rows(pbmc.tsne, pbmc.umap)

# liver

liver.tsne <- get_nonorm_comparisons(stat_results$liver$tsne_conover)
liver.tsne$reduction <- "tSNE"
liver.umap <- get_nonorm_comparisons(stat_results$liver$umap_conover)
liver.umap$reduction <- "UMAP"

liver <- bind_rows(liver.tsne, liver.umap)

# dataset column

breast.group$dataset <- "BC_sub"
breast.dis$dataset <- "BC_dis"
pbmc$dataset <- "PBMC"
liver$dataset <- "Liver"

overall <- bind_rows(breast.group, breast.dis, pbmc, liver)
overall$dataset <- factor(overall$dataset, levels = c("BC_dis", "BC_sub", "PBMC", "Liver"))

#swap elements for visualization

where_nonorm <- which(overall$Technique2 == "nonorm")
overall[where_nonorm, ]$Technique2 <- overall[where_nonorm, ]$Technique1
overall[where_nonorm, ]$Technique1 <- "nonorm"

where_sqrt <- which(overall$Technique2 == "sqrt" & overall$Technique1 == "scran")
overall[where_sqrt, ]$Technique2 <- overall[where_sqrt, ]$Technique1
overall[where_sqrt, ]$Technique1 <- "sqrt"

where_sqrt <- which(overall$Technique2 == "sqrt" & overall$Technique1 == "scnorm")
overall[where_sqrt, ]$Technique2 <- overall[where_sqrt, ]$Technique1
overall[where_sqrt, ]$Technique1 <- "sqrt"

where_sqrt <- which(overall$Technique2 == "sqrt" & overall$Technique1 == "sctransform")
overall[where_sqrt, ]$Technique2 <- overall[where_sqrt, ]$Technique1
overall[where_sqrt, ]$Technique1 <- "sqrt"

where_sqrt <- which(overall$Technique2 == "scran" & overall$Technique1 == "scnorm")
overall[where_sqrt, ]$Technique2 <- overall[where_sqrt, ]$Technique1
overall[where_sqrt, ]$Technique1 <- "scran"

where_nonorm <- which(overall$Technique1 == "dino")
overall[where_nonorm, ]$Technique1 <- overall[where_nonorm, ]$Technique2
overall[where_nonorm, ]$Technique2 <- "dino"

overall$Technique1 <- factor(overall$Technique1, levels = c("nonorm", "log2", "sqrt", "scran", "scnorm", "sctransform"))
overall$Technique2 <- factor(overall$Technique2, levels = c("log2", "sqrt", "scran", "scnorm", "sctransform", "dino"))

# plots

heatmap_conover <- function(data, labs, theme){
  
  color_palette <- rev(MetBrewer::met.brewer("Homer2", type = "continuous", n=6))
  
  p <- ggplot(data, aes(x = Technique1, y = Technique2)) +
    geom_tile(aes(fill = cut(p.adj, 
                             breaks = c(-Inf, 0.00001, 0.0001, 0.001, 0.01, 0.05, 1), 
                             labels = c("<0.00001", "0.00001-0.0001", "0.0001-0.001", "0.001-0.01", "0.01-0.05", ">0.05"))))+
    # geom_point(aes(size = conover.d), 
    #            color = "black", shape = 21, stroke = 1) +
    
    geom_text(aes(label = round(conover.d, 2)), color="white") +
    
    # Define discrete fill scale
    scale_fill_manual(
      name = "Adjusted p-value",
      values = c("<0.00001" = color_palette[1],
                 "0.00001-0.0001" = color_palette[2],
                 "0.0001-0.001" = color_palette[3],
                 "0.001-0.01" = color_palette[4],
                 "0.01-0.05" = color_palette[5],
                 ">0.05" = color_palette[6])
    ) +
    
    guides(
      fill = guide_legend(
        override.aes = list(size = 6)  # Increase the legend point size
      )
    ) +
    
    # Define size scale for Conover's D
    scale_size_continuous(
      range = c(3, 8), 
      name = "Effect Size (Conover's D)",
      breaks = c(0, 0.2, 0.5, 0.8), 
      labels = c("<0.2 - extra small", "0.2-0.5 - small", "0.5-0.8 - medium", ">0.8 - large")
    ) +
    
    theme + labs + facet_wrap(reduction~dataset, nrow = 2)
  
  return(p)
}

labs <- labs(
  title = NULL,
  x = NULL,
  y = NULL
) 

theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
  axis.text.y = element_text(hjust = 1, size = 10, face = "bold"), 
  panel.grid = element_blank(),        
  panel.background = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(face = "bold"),
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  plot.title = element_text(hjust = 0.5, vjust = 0.5),
  plot.background = element_blank(),  
)

output_heatmap <- heatmap_conover(overall, labs, theme)
output_heatmap

width = 32
height = 12

units = "cm"

# ggsave(
#   filename = "article_plot.png",
#   plot = output_heatmap,
#   width = 8,
#   height = 5,
#   units = "in"
# )
