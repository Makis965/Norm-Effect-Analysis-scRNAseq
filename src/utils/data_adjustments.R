adjust_data <- function(data){
  output_data <- data[!data$norm == "scaled", ]
  
  levels <- c("nonorm", "log2", "sqrt", "scran", "scnorm", "sctransform", "dino")
  
  output_data$norm <- factor(output_data$norm, levels=levels)
  
  new_colnames <- colnames(output_data)
  new_colnames[1:2] <- c("Dim1", "Dim2")
  colnames(output_data) <- new_colnames
  
  output_data[output_data$reduction == "tsne", "reduction"] <- "tSNE"
  
  return(output_data)
}

mean_SI <- function(data, x_pos=-35, y_pos=40){
  output_data <- data %>%
    group_by(norm, reduction) %>%
    summarize(Dim1=x_pos, Dim2=y_pos, SI=mean(si))
  
  return(output_data)
}