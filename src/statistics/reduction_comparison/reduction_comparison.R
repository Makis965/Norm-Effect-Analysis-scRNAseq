library(config)
library(ggplot2)
library(ggpubr)
library(config)
library(dplyr)
library(MetBrewer)

config <- config::get()

source(config$utils$statistics)
source(config$utils$measures)
source(config$utils$data_adjustments)

load(config$results$clustering$breast$disease)
load(config$results$clustering$breast$group)
load(config$results$clustering$pbmc)
load(config$results$clustering$liver)

breast_disease_clustering_meas$dataset <- "breast dis"
breast_groups_clustering_meas$dataset <- "breast sub"

overall <- bind_rows(
  breast_disease_clustering_meas,
  breast_groups_clustering_meas,
  pbmc_clustering_meas,
  liver_clustering_meas
)

overall[overall$types == "adjrand", "types"] = "ARI"
overall[overall$types == "mi", "types"] = "MI"
overall[overall$types == "sdc", "types"] = "Dice"
overall[overall$types == "si", "types"] = "SI"

# overall_si <- overall[overall$types == "si", ]

theme <- theme(
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.grid = element_blank(),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(face = "bold")
)

guides <- guides(fill = "none")

labs <- labs(x="", y="")

width = 1300
heigth = 450

#reductions
fill_palette <- rev(MetBrewer::met.brewer("Homer2", type = "continuous", n=2))
color_palette <- rev(MetBrewer::met.brewer("Homer1", type = "continuous", n=4))

ggplot(overall)+
  geom_boxplot(aes(x=reduction, y=scores, fill=reduction))+
  geom_point(aes(x=reduction, y=scores, color=dataset))+
  stat_compare_means(aes(x=reduction, y=scores), method="wilcox.test", vjust=-1, hjust=0, size=3)+
  facet_wrap(~types, scales = "free", ncol = 4)+theme+guides+labs+
  scale_fill_manual(values = fill_palette)+
  scale_color_manual(values = color_palette)

ggsave("reductions_boxplots.png", reductions, height = heigth, width = width, units = "px")

#clusterings
fill_palette <- rev(MetBrewer::met.brewer("Homer2", type = "continuous", n=3))
color_palette <- rev(MetBrewer::met.brewer("Homer1", type = "continuous", n=4))

ggplot(overall)+
  geom_boxplot(aes(x=clustering, y=scores, fill=clustering))+
  geom_point(aes(x=clustering, y=scores, color=dataset))+
  stat_compare_means(aes(x=clustering, y=scores), method="kruskal.test", vjust=-1, hjust=0.2, size=3)+
  facet_wrap(~types, scales = "free", ncol = 4)+theme+guides+labs+
  scale_fill_manual(values = fill_palette)+
  scale_color_manual(values = color_palette)
