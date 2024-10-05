process_clustered_data <- function(data_sets, decode_labels, label_type, meta_labels){
  processed_data = list()
  for(set in data_sets){
    
    data <- base::get(set)
    
    SI <- silhouette(decode_labels, stats::dist(data))
    processed_data[[set]] = cbind(data, SI, meta_labels=meta_labels)
  } 
  return(processed_data)
}