#' Modal length from length data
#'
#' Calculates the annual modal length from the length data, which can be used to select Lc. 
#' Note: Modal length can change over time for many reasons, including a change in mortality 
#' (Hordyk et al. 2015), recruitment, or selectivity (Huynh et al. 2017).
#'
#' @param MLZ_data An object of class \code{MLZ_data}.
#' @param length.slot Name of slot in \code{\linkS4class{MLZ_data}} from which to calculate
#' modal lengths, either: \code{Len_df} or \code{Len_matrix}. Only used if there are data in both slots.
#' @param breaks Only used for \code{Len_df}. An optional vector for breaks for \code{\link{bin_length}}.
#' @param figure If TRUE, a plot is also drawn.
#' @return A data frame of plotted values.
#' @details 
#' Length frequency matrix from \code{Len_df} are created by using \code{\link[graphics]{hist}} function.
#' 
#' @examples
#' \dontrun{
#' data(Nephrops)
#' modal_length(Nephrops)
#' 
#' data(SilkSnapper)
#' new.dataset <- new("MLZ_data", Year = 1983:2013, Len_df = SilkSnapper)
#' modal_length(new.dataset)
#' modal_length(new.dataset, breaks = seq(80, 830, 10))
#' }
#' @references 
#' Hordyk, A. Ono, K., Sainsbury, K., Loneragan, N., and Prince, J. 2015. Some explorations of the 
#' life history ratios to describe length composition, spawning-per-recruit, and the 
#' spawning potential ratio. ICES Journal of Marine Science 72:204-216.
#' 
#' Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. 2017. Multispecies Extensions
#' to a Nonequilibrium Length-Based Mortality Estimator. Marine and Coastal Fisheries 9:68-78.
#' 
#' @export
modal_length <- function(MLZ_data, length.slot = c("Len_df", "Len_matrix"), breaks = NULL, figure = TRUE) {
  length.slot <- match.arg(length.slot)
  length.units <- MLZ_data@length.units
  if(length(length.units) > 0) length.units <- paste0("(", length.units, ")") else length.units <- NULL
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(list = old_par), add = TRUE)
  
  if(nrow(MLZ_data@Len_df) == 0 & ncol(MLZ_data@Len_matrix) == 0 & nrow(MLZ_data@Len_matrix) == 0) stop("No length data available.")
  if(length.slot == "Len_df" && nrow(MLZ_data@Len_df) > 0) {
    message("Data from slot @Len_df used to plot modal lengths.")
    
    Len_df <- MLZ_data@Len_df
    names(Len_df) <- c("Year", "Length")
    if(is.null(breaks)) breaks <- hist(Len_df$Length, plot = FALSE, right = FALSE)$breaks
    
    Len.output <- bin_length(Len_df, breaks = breaks)
    Len_matrix <- Len.output$Len_matrix
    Len_bins <- Len.output$Len_bins
    Year <- as.numeric(Len.output$Year)
  }
  
  if(length.slot == "Len_matrix" || nrow(MLZ_data@Len_df) == 0) {
    message("Data from slot @Len_matrix used to calculate mean lengths.")
    Len_matrix <- MLZ_data@Len_matrix
    Len_bins <- MLZ_data@Len_bins
    Year <- MLZ_data@Year
  }
  
  all.na.index <- apply(Len_matrix, 1, function(x) all(is.na(x)))
  all.zero.index <- rowSums(Len_matrix) == 0
  
  Len_matrix[is.na(Len_matrix)] <- 0
  annual_modal_index <- apply(Len_matrix, 1, which.max)
  modal <- Len_bins[annual_modal_index]
  modal[all.na.index] <- NA
  modal[all.zero.index] <- NA
  modal.df <- data.frame(Year = Year, ModalLength = modal)
  
  if(figure) {
    plot(Year, modal, ylab = paste("Modal Length", length.units), pch = 16, typ = "o")
    if(length(MLZ_data@Lc) == 0 || is.null(MLZ_data@Lc)) Lc <- NULL else Lc <- MLZ_data@Lc
    if(!is.null(Lc)) {
      abline(h = Lc, col = "red")
      legend("topleft", expression(L[c]), col = "red", lty = 1)
    }
  }
  
  return(modal.df)
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
#' \dontrun{
#' data(Nephrops)
#' Nephrops <- calc_ML(Nephrops, sample.size = FALSE)
#' Nephrops@MeanLength
#' plot(Nephrops)
#' }
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
    if(length(MLZ_data@Year) == 0) MLZ_data@Year <- time.series$Year
  }

  if(length.slot == "Len_matrix" || nrow(MLZ_data@Len_df) == 0) {
    message("Data from slot @Len_matrix used to calculate mean lengths.")
    if(length(MLZ_data@Year) != nrow(MLZ_data@Len_matrix)) stop("Length of vector MLZ_data@Year is not equal to the number of rows in MLZ_data@Len_matrix.")

    truncate.indices <- MLZ_data@Len_bins >= MLZ_data@Lc
    Len_matrix <- MLZ_data@Len_matrix[, truncate.indices]
    num <- sweep(Len_matrix, 2, MLZ_data@Len_bins[truncate.indices], "*")
    denom <- rowSums(Len_matrix)
    if(sample.size) ss <- denom
    else ss <- ifelse(is.na(denom) | denom == 0, 0, 1)

    time.series <- data.frame(Year = MLZ_data@Year, MeanLength = rowSums(num)/denom, ss = ss)
  }

  MLZ_data@MeanLength <- ifelse(is.nan(time.series$MeanLength), NA, time.series$MeanLength)
  MLZ_data@ss <- time.series$ss

  return(MLZ_data)
}

#' Bin length data
#'
#' A tool to bin raw length observations into a length frequency matrix.
#'
#' @param df A data frame or matrix of length observations. The first column should be named 
#' 'Year' and the second column should be named 'Length'.
#' @param breaks An optional vector for breaks for \code{\link[graphics]{hist}}.
#' 
#' @details 
#' Length frequencies from \code{Len_df} are created by using \code{\link[graphics]{hist}} function.
#' 
#' @return A list with length bins, years, and frequency matrix.
#' @examples
#' \dontrun{
#' data(SilkSnapper)
#' Silk.matrix <- bin_length(SilkSnapper, breaks = seq(80, 830, 10))
#' Silk.matrix <- bin_length(SilkSnapper)
#' new.dataset <- new("MLZ_data", Year = Silk.matrix$Year, Len_bins = Silk.matrix$Len_bins,
#' Len_matrix = Silk.matrix$Len_matrix)
#' } 
#' @export
bin_length <- function(df, breaks = NULL) {
  df <- as.data.frame(df)
  names(df) <- c("Year", "Length")
  if(is.null(breaks)) breaks <- hist(df$Length, plot = FALSE)$breaks
  Len_filter <- split(df$Length, df$Year)
  Len_hist <- lapply(Len_filter, hist, plot = FALSE, breaks = breaks, right = FALSE)
  
  Year <- new("list", as.numeric(names(Len_filter)))
  Len_summary <- Map(bin_length.Len_summary, x = Year, y = Len_hist)
  Len_summary <- do.call(rbind, Len_summary)
  Len_matrix <- acast(Len_summary, Year ~ Length, value.var = "Frequency")
  output <- list(Len_bins = unique(Len_summary$Length), 
                 Year = unique(Len_summary$Year),
                 Len_matrix = Len_matrix)
  return(output)
}
