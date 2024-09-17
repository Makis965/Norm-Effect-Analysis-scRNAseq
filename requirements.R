install_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
}

# ---- basic ---- 

packages <- c(
    "config",
    "dplyr",
    "uwot"
)

sapply(packages, install_missing)

# ---- additional ---- 

packages <- c(
    "ggplot2",
)

sapply(packages, install_missing)