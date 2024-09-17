library(config)

clear_workspace <- function() {
  rm(
    cells.nonorm,
    cells.log2,
    cells.sqrt,
    cells.scran,
    cells.scnorm,
    cells.seurat,
    cells.dino,
    config,
    data,
    tsne_results,
    set,
    clear_workspace,
    envir = .GlobalEnv
  )
}
