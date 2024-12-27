library(config)
library(ggplot2)
library(ggpubr)
library(config)
library(dplyr)
library(MetBrewer)

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

breast_disease_means <- breast_disease_disp_vals %>%
  group_by(Disease, norm) %>%
  summarize(SI = mean(si))

breast_group_means <- breast_group_disp_vals %>%
  group_by(group, norm) %>%
  summarize(SI = mean(si))

liver_disp_means <- liver_disp_vals %>%
  group_by(Cell.type.org, norm) %>%
  summarize(SI = mean(si))

#Conover outcomes

load(config$results$statistics$conover)

get_pvals <- function(data_tsne, data_umap){
  
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
    rownames(data) <- NULL
    # return(data[data$Technique1 == "nonorm" | data$Technique2 == "nonorm", ])
    return(data)
    
  }
  
  data_tsne <- get_nonorm_comparisons(data_tsne)
  data_tsne$reduction <- "tSNE"
  data_umap <- get_nonorm_comparisons(data_umap)
  data_umap$reduction <- "UMAP"
  merged <- bind_rows(data_tsne, data_umap)
  
  return(merged)
  
}

pairwise_df_pbmc <- get_pvals(stat_results$pbmc$tsne_conover, stat_results$pbmc$umap_conover)
pairwise_df_pbmc <- pairwise_df_pbmc[, c("Technique1", "Technique2", "p.adj")]
colnames(pairwise_df_pbmc) <- c("group1", "group2", "p.value")

# plots

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
  axis.text.x = element_text(face = "bold", size = 9),
  axis.text.y = element_text(face = "bold", size = 9)
)

guides <- guides(fill = "none")

labs <- labs(x="", y="")

plot_violins <- function(data, disp_means, color_col, pairwise){
  
  fill_palette <- rev(MetBrewer::met.brewer("Homer2", type = "continuous", n=7))
  color_palette <- MetBrewer::met.brewer("Redon", type = "continuous", n=13)
    
  ggplot(data)+
    geom_violin(aes(x=norm, y=si, fill=norm), alpha=0.5)+
    geom_point(data=disp_means, aes(x=norm, y=SI, color=!!rlang::sym(color_col)), alpha=1, size=2)+
    facet_wrap(~reduction)+
    ggpubr::stat_compare_means(aes(x=norm, y=si), method = "kruskal.test", vjust=-1, hjust=0)+
    scale_fill_manual(values = fill_palette) +
    scale_color_manual(values = color_palette) + 
    theme + guides + labs 
    
    # stat_pvalue_manual(
    #   pairwise,
    #   label = "p.value",
    #   xmin = "group1",
    #   xmax = "group2",
    #   y.position = seq(max(data$si, na.rm = TRUE) * 1.1, length.out = nrow(pairwise)),
    # )
  
}

width = 32
height = 12

units = "cm"

# pbmc
p1 <- plot_violins(pbmc_disp_vals, pbmc_disp_means, "CellType", pairwise_df_pbmc)
p1

# ggsave(filename = "violin_SI_pbmc.png", p1, width = width, height = height, units = units)

# breast disease
p2 <- plot_violins(breast_disease_disp_vals, breast_disease_means,  "Disease")
p2

# ggsave(filename = "violin_SI_breast_d.png", p2, width = width, height = height, units = units)

# breast groups
p3 <- plot_violins(breast_group_disp_vals, breast_group_means, "group")
p3

# ggsave(filename = "violin_SI_breast_g.png", p3, width = width, height = height, units = units)

# liver
p4 <- plot_violins(liver_disp_vals, liver_disp_means, "Cell.type.org")
p4

# ggsave(filename = "violin_SI_liver.png", p4, width = width, height = height, units = units)
