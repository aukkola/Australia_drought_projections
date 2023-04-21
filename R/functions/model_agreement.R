model_agreement <- function(x) {
  
  #Calculate fraction of models that indicate a positive change
  #take NA values into account by removing these from the calculation
  if (all(is.na(x))) {
    return(NA)
  } else {
    positive <- length(which(x > 0))
    fraction <- positive / length(which(!is.na(x))) * 100
    return(fraction) 
  }
}
