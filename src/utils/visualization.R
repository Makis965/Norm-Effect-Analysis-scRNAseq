


# ---- 2dims plots ----

plot_point <- function(data, reduction, color=NULL, fill=NULL, x_pos=-35, y_pos=40){
  
  plot_data <- mean_SI(data, x_pos, y_pos)
  
  plot <- ggplot(data[data$reduction == reduction, ])+
    geom_point(aes(x=Dim1, y=Dim2, color=!!sym(color)), alpha=0.5)+
    geom_text(data = plot_data[plot_data$reduction == reduction, ],
              aes(x = Dim1, y = Dim2, label = paste("SI =",round(SI, 3))),
              color = "black",
              fontface = "italic")+
    facet_wrap(~norm)
  
  return(plot)
}