#' MLZ_data
#'
#' An S4 class for storing data and life history parameters for a single stock.
#' Method functions \code{summary} and \code{plot} are available for this class (see examples).
#'
#' @slot Stock Name of stock.
#' @slot Year A vector of years to be considered in the model. Need to be increments of a single year.
#' @slot Len_df A data frame containing individual length observations. The first column should be the Year and the second column should be the length.
#' @slot Len_bins A vector of midpoints of length bins for \code{Len_matrix}.
#' @slot Len_matrix A matrix of size data. Rows index \code{Year} and columns index \code{Len_bins}.
#' @slot vbLinf L-infinity from the von Bertalanffy growth function.
#' @slot vbK Parameter K from the von Bertalanffy growth function.
#' @slot vbt0 Parameter t0 from the von Bertalanffy growth function.
#' @slot Lc Length of full selectivity.
#' @slot M Natural mortality rate.
#' @slot lwb Exponent \code{b} from the allometric length-weight function \eqn{W = aL^b}.
#' @slot MeanLength Vector of mean lengths of animals larger than Lc.
#' @slot ss Annual sample sizes for MeanLength.
#' @slot CPUE Vector of catch-per-unit-effort data.
#' @slot Effort Vector of effort data.
#' @slot length.units Unit of measurement for lengths, i.e. "cm" or "mm".
#'
#' @examples
#' data(Goosefish); Goosefish
#' summary(Goosefish)
#' plot(Goosefish)
#'
#' new("MLZ_data")
#' @export
#' @import ggplot2
#' @import dplyr
#' @import parallel
#' @import graphics
#' @importFrom reshape2 acast melt
#' @importFrom methods new
#' @importFrom gplots rich.colors
#' @importFrom Rcpp sourceCpp
#' @importFrom stats optim cov2cor
#' @importFrom grDevices cm.colors
#' @importFrom numDeriv hessian grad
#' @importFrom msm deltamethod
#' @useDynLib MLZ
setClass("MLZ_data", slots = c(Stock = "character", Year = "vector", Len_df = "data.frame", Len_bins = "vector",
                               Len_matrix = "matrix", vbLinf = "numeric", vbK = "numeric", vbt0 = "numeric",
                               M = "numeric", Lc = "numeric", lwb = "numeric", MeanLength = "numeric", ss = "numeric",
                               CPUE = "numeric", Effort = "numeric", length.units = "character"))

#' MLZ_model
#'
#' An S4 class for storing model results.
#' Method functions \code{summary} and \code{plot} are available for this class (see examples).
#'
#' @slot Stock Name of stock (obtained from an object of class \code{MLZ_data}).
#' @slot Model Name of model used for mortality estimation.
#' @slot time.series A data frame summarizing observed time series data and predicted values from model.
#' @slot estimates A matrix of parameter estimates and derived values and their standard errors. Standard errors
#' of derived quantities are obtained from call to \code{\link[msm]{deltamethod}}.
#' @slot negLL The negative log-likelihood from the model.
#' @slot n.changepoint The number of change points in the model.
#' @slot n.species The number of species/stocks in the model.
#' @slot grid.search A data frame reporting the log-likelihood values from a grid search over change points.
#' See \code{\link{profile_ML}}, \code{\link{profile_MLCR}}, and \code{\link{profile_MLmulti}}.
#' @slot opt A list with components from \code{\link[stats]{optim}}, including
#' correlation matrix and gradient vector at estimated minimum from call to \code{\link[numDeriv]{hessian}} and
#' \code{\link[numDeriv]{grad}}, respectively.
#' @slot length.units Unit of measurement for lengths, i.e. "cm" or "mm".
#'
#' @examples
#' data(Goosefish)
#' goose.model <- ML(Goosefish, ncp = 2, grid.search = FALSE, figure = FALSE)
#' class(goose.model)
#'
#' summary(goose.model)
#' plot(goose.model, residuals = FALSE)
#'
#' @export
setClass("MLZ_model", slots = c(Stock = "character", Model = "character", time.series = "data.frame",
                                estimates = "matrix", negLL = "numeric",
                                n.changepoint = "integer", n.species = "integer", grid.search = "data.frame",
                                opt = "list", length.units = "character"))

#' \code{summary} method for S4 class \code{MLZ_data}
#'
#' @param object An object of class \code{MLZ_data}.
#' @examples
#' data(MuttonSnapper)
#' summary(MuttonSnapper)
#'
#' @export
setMethod("summary", signature(object = "MLZ_data"), function(object) {
  MLZ_data <- object
  time.series.df <- try(data.frame(Year = MLZ_data@Year, MeanLength = MLZ_data@MeanLength,
                          ss = MLZ_data@ss))
  if(inherits(time.series.df, "try-error")) {
    plot(MLZ_data)
    stop("Could not print a data.frame summary. Check length of vectors in slots @Year,
           @MeanLength, and @ss.")
  } else {
    if(length(MLZ_data@CPUE) > 0) time.series.df$CPUE <- MLZ_data@CPUE
    if(length(MLZ_data@Effort) > 0) time.series.df$Effort <- MLZ_data@Effort
    return(time.series.df)
  }

})

#' \code{summary} method for S4 class \code{MLZ_model}
#'
#' @param object An object of class \code{MLZ_model}.
#' @examples
#' /dontrun{
#' data(Goosefish)
#' goose.model <- ML(Goosefish, ncp = 2, grid.search = FALSE)
#' summary(goose.model)
#' }
#'
#' @export
setMethod("summary", signature(object = "MLZ_model"), function(object)
  list(Stock = object@Stock, Model = object@Model, Estimates = object@estimates))







#' \code{plot} method for S4 class \code{MLZ_data}
#'
#' Plots annual length frequencies from slot \code{Len_matrix}.
#'
#' @param x An object of class \code{MLZ_data}.
#' @param ggplot_layer Layers to add to ggplot2 plot.
#' @examples
#' data(Nephrops)
#' plot(Nephrops)
#'
#' @export
setMethod("plot", signature(x = "MLZ_data"), function(x, ggplot_layer = NULL) {
  MLZ_data <- x

  old_par <- par(no.readonly = TRUE)
  on.exit(par(list = old_par), add = TRUE)

  length.units <- MLZ_data@length.units
  length.units <- paste0("(", length.units, ")")

  no.Len_matrix <- nrow(MLZ_data@Len_matrix) == 0 & ncol(MLZ_data@Len_matrix) == 0
  no.Len_df <- nrow(MLZ_data@Len_df) == 0
  no.ML <- length(MLZ_data@MeanLength) == 0

  if(no.Len_matrix & no.Len_df & no.ML) stop("No length data available.")

  if(!no.Len_matrix) {
    Len_matrix <- MLZ_data@Len_matrix
    rownames(Len_matrix) <- MLZ_data@Year
    colnames(Len_matrix) <- MLZ_data@Len_bins
    Len_matrix <- melt(Len_matrix)
    names(Len_matrix) <- c("Year", "Length", "Frequency")

    if(length(MLZ_data@Lc) != 0 || !is.null(MLZ_data@Lc)) z <- geom_vline(xintercept = MLZ_data@Lc, colour = "red")
    else z <- NULL
    ggplot_layer <- sapply(ggplot_layer, eval)

    theme_clean <- theme_bw() + theme(panel.spacing = unit(0, "inches"),
                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank())

    zz <- ggplot(Len_matrix, aes(x = Length, y = Frequency)) + geom_line() + facet_wrap(~ Year) + theme_clean + z +
      ggplot_layer

  }
  if(!no.Len_df) {
    if(!no.Len_matrix) message("Data from slot @Len_matrix was plotted.")
    if(no.Len_matrix) {
      Len_df <- MLZ_data@Len_df
      names(Len_df) <- c("Year", "Length")

      if(length(MLZ_data@Lc) != 0 || !is.null(MLZ_data@Lc)) z <- geom_vline(xintercept = MLZ_data@Lc, colour = "red")
      else z <- NULL
      ggplot_layer <- sapply(ggplot_layer, eval)

      theme_clean <- theme_bw() + theme(panel.spacing = unit(0, "inches"),
                                        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

      zz <- ggplot(Len_df, aes(x = Length)) + geom_histogram() + facet_wrap(~ Year) + theme_clean + z +
        ggplot_layer
    }
  }
  if(!no.ML) {
    nplots <- 1
    summary.MLZ <- summary(MLZ_data)
    if("Effort" %in% names(summary.MLZ)) nplots <- nplots + 1
    if("CPUE" %in% names(summary.MLZ)) nplots <- nplots + 1
    if(nplots == 3) layout(matrix(c(1,1,1,1,2,3), nrow = 2))
    if(nplots < 3) par(mfrow = c(nplots, 1))
    plot(MeanLength ~ Year, data = summary.MLZ, pch = 16, typ = "o",
         xlab = "Year", ylab = paste("Mean Length (> Lc)", length.units))
    if("CPUE" %in% names(summary.MLZ)) plot(CPUE ~ Year, data = summary.MLZ, pch = 16, typ = "o",
                                            xlab = "Year", ylab = "CPUE")
    if("Effort" %in% names(summary.MLZ)) plot(Effort ~ Year, data = summary.MLZ, pch = 16, typ = "o",
                                              xlab = "Year", ylab = "Effort")
  }
  if(exists("zz")) return(zz) else invisible()
}
)

#' \code{plot} method for S4 class \code{MLZ_model}
#'
#' Plots time series of observed and predicted data from an object of class \code{MLZ_model}.
#'
#' @param x An object of class \code{MLZ_model}.
#' @param residuals logical; whether a plot of residuals will also be produced.
#' @examples
#' data(Goosefish)
#' goose.model <- ML(Goosefish, ncp = 2, grid.search = FALSE, figure = FALSE)
#' plot(goose.model)
#' @export
setMethod("plot", signature(x = "MLZ_model"), function(x, residuals = TRUE) {
  MLZ_model <- x
  length.units <- MLZ_model@length.units

  old_par <- par(no.readonly = TRUE)
  on.exit(par(list = old_par), add = TRUE)

  if(MLZ_model@Model == "ML") {
    if(residuals) par(mfrow = c(1, 2))
    ts <- MLZ_model@time.series
    plot(MeanLength ~ Year, ts, ylab = paste("Mean Length (> Lc)", length.units), typ = "o", pch = 16)
    lines(Predicted ~ Year, ts, col = "red", lwd = 2)
    if(residuals) {
      plot(Residual ~ Year, ts, typ = "o", pch = 16)
      abline(h = 0, lwd = 2)
    }
  }

  if(MLZ_model@Model == "MLCR") {
    if(residuals) par(mfrow = c(2, 2), mar = c(5,4,1,1))
    else par(mfrow = c(2,1))
    ts <- MLZ_model@time.series
    plot(MeanLength ~ Year, ts, ylab = paste("Mean Length (> Lc)", length.units), typ = "o", pch = 16)
    lines(Predicted.ML ~ Year, ts, col = "red", lwd = 2)
    if(residuals) {
      plot(Residual.ML ~ Year, ts, typ = "o", pch = 16, ylab = "Mean Length Residual")
      abline(h = 0, lwd = 2)
    }
    plot(CPUE ~ Year, ts, typ = "o", pch = 16)
    lines(Predicted.CPUE ~ Year, ts, col = "red", lwd = 2)
    if(residuals) {
      plot(Residual.CPUE ~ Year, ts, typ = "o", pch = 16, ylab = "CPUE Residual")
      abline(h = 0, lwd = 2)
    }
  }
  if(MLZ_model@Model == "MLeffort") {
    if(residuals) layout(matrix(c(1,1,1,1,2,3), nrow = 2))
    else par(mfrow = c(1,2))
    par(mar = c(4, 4, 1, 1))
    ts <- MLZ_model@time.series

    plot(MeanLength ~ Year, ts, ylab = paste("Mean Length (> Lc)", length.units), typ = "o", pch = 16)
    lines(Predicted ~ Year, ts, col = "red", lwd = 2)
    if(residuals) {
      plot(Residual ~ Year, ts, typ = "o", pch = 16)
      abline(h = 0, lwd = 2)
    }
    max.mort <- 1.1 * max(c(ts$F, ts$Z))
    plot(F ~ Year, ts, ylab = "Mortality", lty = 2, typ = "o", pch = 16, ylim = c(0, max.mort), col = "red")
    lines(Z ~ Year, ts, typ = "o", pch = 16, col = "red")
    legend("topright", c("Z", "F"), lty = c(1:2), pch = 16, col = "red")

  }

  if(MLZ_model@Model == "MLmulti") {
    par(mfrow = c(2,2))
    ts <- MLZ_model@time.series
    nspec <- MLZ_model@n.species
    nyrs <- nrow(ts)/nspec
    if(length(MLZ_model@Stock) == nspec) spec.name <- MLZ_model@Stock
    else {
      warning("Number of stock names in MLZ_model@Stock not equal to number of species.")
      spec.name <- paste0("Species_", 1:nspec)
    }
    for(i in 1:nspec) {
      plot(MeanLength ~ Year, ts[((i-1)*nyrs+1):(i*nyrs), ],
           ylab = paste("Mean Length (> Lc)", length.units), typ = "o", pch = 16)
      lines(Predicted ~ Year, ts[((i-1)*nyrs+1):(i*nyrs), ], lwd = 2)
      title(spec.name[i])
      if(residuals) {
        plot(Residual ~ Year, ts[((i-1)*nyrs+1):(i*nyrs), ],
             typ = "o", pch = 16)
        abline(h = 0, lwd = 2)
      }
    }

  }

  invisible()

}
)
