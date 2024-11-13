install_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

install_missing_biocm <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- basic ---- 

packages <- c(
  #base packages
  "BiocManager",
  "config",
  "dplyr",
  
  #normalization packages
  "Seurat",
  
  #dim reduction packages
  "Rtsne",
  "uwot",
  
  #clustering packages
  "cluster",
  "dynamicTreeCut",
  "SeuratObject",
  "igraph",
  
  #visualization packages
  "ggplot2"
)

packages_bioc <- c(
  #normalization packages
  "Dino",
  "SCnorm",
  "scran"
  )

sapply(packages, install_missing)
sapply(packages_bioc, install_missing_biocm)
