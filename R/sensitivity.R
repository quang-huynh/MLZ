#' Sensitivity to Lc
#' 
#' The function re-calculates mean lengths for each alternative value of Lc and 
#' re-estimates mortality. Currently supports only the ML estimator.
#' 
#' @param MLZ_data An object of class \code{\linkS4class{MLZ_data}} containing mean lengths and
#' life history data of stock. Must contain length composition data.
#' @param MLZ_model An object of class \code{\linkS4class{MLZ_model}} with base value of Lc.
#' @param Lc.vec A vector of alternative Lc values.
#' @param grid.search Whether a grid search is performed or not. By default, the starting values 
#' in the sensitivity analysis are the estimates from object \code{MLZ_model}.
#' @param figure Whether a figure will be produced, similar to Figure 6 of Huynh et al. (2017).
#' 
#' @return A matrix of mortality and change point estimates with each value Lc.
#' @references 
#' Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. 2017. Multispecies Extensions
#' to a Nonequilibrium Length-Based Mortality Estimator. Marine and Coastal Fisheries 9:68-78.
#'
#' @examples
#' data(SilkSnapper)
#' new.dataset <- new("MLZ_data", Year = 1983:2013, Len_df = SilkSnapper, length.units = "mm",
#' vbLinf = 794, vbK = 0.1)
#' 
#' new.dataset@Lc <- 310
#' new.dataset <- calc_ML(new.dataset)
#' 
#' first.MLZmodel <- ML(new.dataset, 1)
#' Lc.vec <- seq(240, 340, 5)
#' 
#' sensitivity_Lc(new.dataset, first.MLZmodel, Lc.vec)
#' 
#' @seealso \code{\link{ML}}
#' @export
sensitivity_Lc <- function(MLZ_data, MLZ_model, Lc.vec, grid.search = FALSE, figure = TRUE) {
  if(length(MLZ_data@Len_df) == 0 & (nrow(MLZ_data@Len_matrix) == 0 & ncol(MLZ_data@Len_matrix) == 0)) 
    stop("No length composition data available to re-calculate mean lengths with Lc.")
  #if(MLZ_data@Stock != MLZ_model@Stock) warning("Stock names do not match between data and model.")
  Lc.base <- MLZ_data@Lc
  #if(Lc.base %in% Lc.vec) Lc.vec <- Lc.vec[-which(Lc.vec == Lc.base)]
  if(MLZ_model@Model != "ML") stop("Functions other than ML currently not supported yet.")
  
  ncp <- MLZ_model@n.changepoint
  Z.ind <- startsWith(rownames(MLZ_model@estimates), "Z")
  if(ncp > 0) yearZ.ind <- startsWith(rownames(MLZ_model@estimates), "y")
  if(ncp == 0 || grid.search) start <- NULL
  if(ncp > 0 & !grid.search) 
    start <- list(Z = MLZ_model@estimates[Z.ind, 1], 
                  yearZ = MLZ_model@estimates[yearZ.ind, 1])
  
  res <- list()
  for(i in 1:length(Lc.vec)) {
    MLZ.new <- MLZ_data
    MLZ.new@Lc <- Lc.vec[i]
    MLZ.new <- calc_ML(MLZ.new)
    res[[i]] <- ML(MLZ.new, ncp = ncp, start = start, spawn = attr(MLZ_model, "spawn"),
              grid.search = grid.search, figure = FALSE)
  }
  
  output <- lapply(res, function(x) x@estimates[, 1])
  output <- do.call(rbind, output)
  
  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(list = old_par), add = TRUE)
    par(las = 1)
    if(ncp > 0) par(mfrow = c(1,2))
    
    Z <- matrix(output[, Z.ind], nrow = length(Lc.vec))
    Z.max <- range(Z, na.rm = TRUE, finite = TRUE)[2]
    color.vec <- rich.colors(ncp+1)
    
    if(length(MLZ_data@length.units) == 0) Lc.label <- expression(L[c])
    else Lc.label <- parse(text = paste0("L[c]~(", MLZ_data@length.units, ")"))
    
    plot(Lc.vec, Z[, 1], typ = 'o', pch = 16, ylim = c(0, 1.1 * Z.max), lwd = 2, las = 1, 
         xlab = Lc.label, ylab = "Total Mortality Z", col = color.vec[1])
    abline(v = Lc.base, lty = 2)
    if(ncp > 0) {
      for(i in 2:(ncp+1)) lines(Lc.vec, Z[, i], typ = 'o', pch = 16, col = color.vec[i], lwd = 2)
      legend("topright", paste0("Z[", 1:(ncp+1), "]"), lwd = 2, lty = 1, pch = 16, col = color.vec, bty = "n")
      
      yearZ <- matrix(output[, yearZ.ind], nrow = length(Lc.vec))
      yearZ.range <- range(yearZ, na.rm = TRUE, finite = TRUE)
      plot(Lc.vec, yearZ[, 1], typ = 'o', pch = 16, ylim = yearZ.range, lwd = 2,
           las = 1, xlab = Lc.label, ylab = "Change point", col = color.vec[1])
      abline(v = Lc.base, lty = 2)
      if(ncp > 1) {
        for(i in 2:ncp) lines(Lc.vec, yearZ[, i], typ = 'o', pch = 16, col = color.vec[i], lwd = 2)
      }
      legend("topright", paste0("yearZ[", 1:(ncp), "]"), lwd = 2, lty = 1, pch = 16, col = color.vec, bty = "n")
    }
  }
  output <- cbind(Lc = Lc.vec, output)
  return(output)
}
