#Masks pixels where more than 20% of layers have a missing value
mask_data <- function(data) {
  no_na <- length(which(is.na(data)))/length(data)
  if(no_na > 0.2) { return(rep(NA, length(data))) } else { return(data) }
}
