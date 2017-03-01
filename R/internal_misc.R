
timeseries_multi <- function(x) {
  data.frame(Year = x@Year, MeanLength = x@MeanLength, ss = x@ss)
}

returnNAoptvalue <- function(x, reps = 1) {
  if(inherits(x, "try-error")) as.numeric(rep(NA, reps)) else getElement(x, "value")
}
