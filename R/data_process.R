#' Calculate modal length from length data.
#'
#' Creates two figures. First, a histogram of length data pooled across all years is produced and the modal
#' length is reported. Second, a plot of the annual modal length is produced and compared to the pooled modal length.
#'
#' @param MLZ_data An object of class \code{MLZ_data}.
#' @examples
#' data(Nephrops)
#' modal_length(Nephrops)
#' @export
modal_length <- function(MLZ_data) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(list = old_par), add = TRUE)

  par(mfrow = c(2,1), mar = c(4, 4, 1, 1))
  Len_matrix <- MLZ_data@Len_matrix
  L <- colSums(Len_matrix, na.rm = TRUE)

  modal <- MLZ_data@Len_bins[which.max(L)]
  nbins <- length(MLZ_data@Len_bins) - 1
  plot(MLZ_data@Len_bins[1:nbins], L, xlab = "Length", ylab = "Frequency", pch = 16, typ = 'o')
  abline(v = modal, col = 'red')

  Len_matrix[is.na(Len_matrix)] <- 0
  annual_modal_index <- apply(Len_matrix, 1, which.max)
  annual_modal <- MLZ_data@Len_bins[annual_modal_index]
  all.na.index <- apply(MLZ_data@Len_matrix, 1, function(x) all(is.na(x)))
  all.zero.index <- rowSums(MLZ_data@Len_matrix) == 0
  annual_modal[all.na.index] <- NA
  plot(MLZ_data@Year, annual_modal, xlab = "Year", ylab = "Modal Length", pch = 16, typ = 'o')
  abline(h = modal, col = 'red')

  return(modal)
}

#' Calculate mean lengths >= Lc
#'
#' Calculates mean lengths from length data and Lc for class \code{\linkS4class{MLZ_data}}.
#'
#' @param MLZ_data An object of class \code{\linkS4class{MLZ_data}}.
#' @param length.slot Name of slot in \code{\linkS4class{MLZ_data}} from which to calculate
#' mean lengths, either: \code{Len_df} or \code{Len_matrix}. Only used if there are data in both slots.
#' @param sample.size If \code{TRUE}, then the annual sample sizes will be calculated by
#' summing the cells in slot \code{Len_matrix}. Otherwise, sample sizes are set to 0 or 1
#' (whether mean lengths are calculated).
#' @return An object of class \code{\linkS4class{MLZ_data}} to fill slots \code{MeanLength}, \code{ss}.
#' @examples
#' data(Nephrops)
#' Nephrops <- calc_ML(Nephrops, sample.size = FALSE)
#' Nephrops@MeanLength
#' plot(Nephrops)
#' @export
calc_ML <- function(MLZ_data, length.slot = c("Len_df", "Len_matrix"), sample.size = TRUE) {
  length.slot <- match.arg(length.slot)
  if(length(MLZ_data@Lc) == 0 || is.null(MLZ_data@Lc)) stop("MLZ_data@Lc is missing.")
  if(nrow(MLZ_data@Len_df) == 0 & ncol(MLZ_data@Len_matrix) == 0 & nrow(MLZ_data@Len_matrix) == 0) stop("No length data available.")
  if(length.slot == "Len_df" && nrow(MLZ_data@Len_df) > 0) {
    message("Data from slot @Len_df used to calculate mean lengths.")

    names(MLZ_data@Len_df) <- c("Year", "Length")
    Len_filter <- filter(MLZ_data@Len_df, Length >= MLZ_data@Lc)
    time.series <- summarise(group_by(Len_filter, Year), MeanLength = mean(Length), ss = n())
  }

  if(length.slot == "Len_matrix" || nrow(MLZ_data@Len_df) == 0) {
    message("Data from slot @Len_matrix used to calculate mean lengths.")
    if(length(MLZ_data@Year) != nrow(MLZ_data@Len_matrix)) stop("Length of vector MLZ_data@Year is not equal to the number of rows in MLZ_data@Len_matrix.")

    truncate.indices <- MLZ_data@Len_bins >= MLZ_data@Lc
    Len_matrix <- MLZ_data@Len_matrix[, truncate.indices]
    num <- sweep(Len_matrix, 2, MLZ_data@Len_bins[truncate.indices], "*")
    denom <- rowSums(Len_matrix)
    if(sample.size) ss = denom
    else ss = ifelse(is.na(denom) | denom == 0, 0, 1)

    time.series <- data.frame(Year = MLZ_data@Year, MeanLength = rowSums(num)/denom, ss = ss)
  }

  MLZ_data@MeanLength <- time.series$MeanLength
  MLZ_data@ss <- time.series$ss

  return(MLZ_data)
}

