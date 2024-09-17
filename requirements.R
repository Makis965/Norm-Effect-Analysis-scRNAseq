install_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# ---- basic ---- 

packages <- c(
    #base packages
    "config",
    "dplyr",

    #normalization packages
    "Seurat",
    "Dino",
    "scran",
    "SCnorm",

    #dim reduction packages
    "Rtsne",
    "uwot",

    #clustering packages
    "cluster",
    "SeuratObject",
    "igraph",

    #visualization packages
    "ggplot2",
)

sapply(packages, install_missing)
