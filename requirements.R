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
  "config",
  "dplyr",
  
  #normalization packages
  "Seurat",
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
  "ggplot2"
)

packages_bioc <- c("Dino")

sapply(packages, install_missing)
sapply(packages_bioc, install_missing_biocm)
