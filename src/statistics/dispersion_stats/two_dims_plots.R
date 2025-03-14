library(config)
library(ggplot2)
library(ggpubr)
library(config)
library(dplyr)

config <- config::get()

source(config$utils$data_adjustments)
source(config$utils$visualization)

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


# ---- 2-dimensional, dispersion plots ----
background <- "transparent"


theme <- ggplot2::theme(
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "gray90"),  # Light gray major grid
  panel.grid.minor = element_line(color = "gray90"),  # Lighter gray minor grid
  panel.background = element_rect(fill = background, color = NA),
  plot.background = element_rect(fill = background, color = NA),
  
  strip.background = element_rect(fill = background, color = NA),
  strip.text = element_text(face = "bold", size = 12),
  
  # legend.box = "horizontal",
  # legend.position = c(0.6, 0.15),
  legend.title = element_blank(),
  legend.text = element_text(size = 10, face = "bold"),
  
  legend.background = element_rect(fill = "transparent", color = NA), 
  legend.key = element_rect(fill = "transparent", color = NA)
  
  # legend.text = element_blank()        
)

legend <- guides(color = guide_legend(nrow = 3))

#display plots 

#pbmc 
plot_point(data=pbmc_disp_vals, reduction="tSNE", color="CellType", x_pos=-35, y_pos=40)+
  theme + legend

plot_point(data=pbmc_disp_vals, reduction="UMAP", color="CellType", x_pos=-10, y_pos=10)+
  theme + legend

#breast disease
plot_point(data=breast_disease_disp_vals, reduction="tSNE", color="Disease", x_pos=-10, y_pos=20)+
  theme + legend

plot_point(data=breast_disease_disp_vals, reduction="UMAP", color="Disease", x_pos=-10, y_pos=10)+
  theme + legend

#breast group
plot_point(data=breast_group_disp_vals, reduction="tSNE", color="group", x_pos=-10, y_pos=20)+
  theme + legend

plot_point(data=breast_group_disp_vals, reduction="UMAP", color="group", x_pos=-10, y_pos=10)+
  theme + legend

#liver
legend <- guides(color = guide_legend(nrow = 6))
theme$legend.text <- element_text(size = 8, face = "bold")  
theme$legend.position.inside <- c(0.65, 0.15)

plot_point(data=liver_disp_vals, reduction="tSNE", color="Cell.type.org", x_pos=-35, y_pos=40)+
  theme + legend

plot_point(data=liver_disp_vals, reduction="UMAP", color="Cell.type.org", x_pos=-10, y_pos=10)+
  theme + legend

# Selected by Asia

selected = c("log2", "sqrt", "scran", "sctransform")

# theme$legend.text <- element_text(size = 8, face = "bold")  
# theme$legend.position.inside <- c(0.65, 0.15)
legend <- guides(color = guide_legend(ncol = 1))
plot_point(data=pbmc_disp_vals[pbmc_disp_vals$norm %in% selected, ], reduction="UMAP", color="CellType", x_pos=-10, y_pos=10)+
  theme + legend

# theme$legend.text <- element_text(size = 8, face = "bold")  
# theme$legend.position.inside <- c(0.65, 0.15)
plot_point(data=breast_group_disp_vals[breast_group_disp_vals$norm %in% selected, ], reduction="UMAP", color="group", x_pos=-10, y_pos=10)+
  theme + legend

# theme$legend.text <- element_text(size = 8, face = "bold")  
# theme$legend.position.inside <- c(0.65, 0.15)
plot_point(data=liver_disp_vals[liver_disp_vals$norm %in% selected, ], reduction="UMAP", color="Cell.type.org", x_pos=-10, y_pos=10)+
  theme + legend

# tSNE

theme$legend.text <- element_text(size = 8, face = "bold")
theme$legend.position.inside <- c(0.65, 0.15)
legend <- guides(color = guide_legend(ncol = 1))
plot_point(data=pbmc_disp_vals[pbmc_disp_vals$norm %in% selected, ], reduction="tSNE", color="CellType", x_pos=-30, y_pos=35)+
  theme + legend

# theme$legend.text <- element_text(size = 8, face = "bold")  
# theme$legend.position.inside <- c(0.65, 0.15)
plot_point(data=breast_group_disp_vals[breast_group_disp_vals$norm %in% selected, ], reduction="tSNE", color="group", x_pos=-10, y_pos=20)+
  theme + legend

# theme$legend.text <- element_text(size = 8, face = "bold")  
# theme$legend.position.inside <- c(0.65, 0.15)
plot_point(data=liver_disp_vals[liver_disp_vals$norm %in% selected, ], reduction="tSNE", color="Cell.type.org", x_pos=-40, y_pos=40)+
  theme + legend
