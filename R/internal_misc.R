summary_ML <- function(MLZ_data) {
  data.frame(Year = MLZ_data@Year, MeanLength = MLZ_data@MeanLength, ss = MLZ_data@ss)
}

returnNAobjective <- function(x, reps = 1) {
  if(inherits(x, "try-error")) as.numeric(rep(NA, reps)) else getElement(x, "objective")
}

bin_length.Len_summary <- function(x, y) {
  data.frame(Year = rep(x, length(y$mids)), Length = y$mids, Frequency = y$counts)
}

returnMLCRloglike <- function(x) {
  if(inherits(x, "try-error")) as.numeric(rep(NA, 2)) else x$report()$loglike
}
