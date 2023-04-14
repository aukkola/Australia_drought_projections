#Calculate percentage of time in drought
freq <- function(data) {
  sum(data) /
    length(data) * 100
}
