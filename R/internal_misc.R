summary_ML <- function(MLZ_data) {
  data.frame(Year = MLZ_data@Year, MeanLength = MLZ_data@MeanLength, ss = MLZ_data@ss)
}

returnNAobjective <- function(x, reps = 1) {
  if(inherits(x, "try-error")) as.numeric(rep(NA, reps)) else getElement(x, "objective")
}

bin_length.Len_summary <- function(x, y) {
  data.frame(Year = rep(x, length(y$mids)), Length = y$mids, Frequency = y$counts)
}

get_M_MSM1S <- function(MLZ_data) {
  if(length(MLZ_data@M) > 0) M <- MLZ_data@M
  else M <- 0.01
  return(M)
}

# For SSM and MSM1, the Z's are estimated.
# For MSM2 and MSM3, the Z's are derived.
produce_MLmulti_warnings <- function(Z, Z.lower, Z.upper) {
  warning.flag <- warning.flag2 <- logical(length = length(Z.lower))
  for(i in 1:length(Z.lower)) warning.flag[i] <- any(Z[, i] <= Z.lower[i])
  if(any(warning.flag)) warning("There are mortality estimates at boundary (Z = 0.01 or M).")
  
  for(i in 1:length(Z.lower)) warning.flag2[i] <- any(Z[, i] >= Z.upper)
  if(any(warning.flag2)) warning(paste0("There are mortality estimates at boundary (Z = ", Z.upper, ")."))
  invisible()
}
