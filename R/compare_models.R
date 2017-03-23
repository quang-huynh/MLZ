#' Model selection
#'
#' Produces a matrix of AIC for model selection.
#'
#' @param MLZ_model.list A list containing objects of class \code{MLZ_model}, all from the same mortality
#' estimator and same data set.
#' @param figure If \code{TRUE}, produces a figure of model fits to the observed data.
#' @param color Optional vector of colors for the figure each representing a separate model
#' in \code{MLZ_model.list}. If \code{NULL}, colors from \code{\link[gplots]{rich.colors}}
#' will be used.
#'
#' @examples
#' data(Goosefish)
#' goose <- ML(Goosefish, ncp = 0)
#' goose1 <- ML(Goosefish, ncp = 1)
#' goose2 <- ML(Goosefish, ncp = 2, grid.search = TRUE, figure = FALSE)
#'
#' compare_models(list(goose, goose1, goose2))
#'
#' data(PRSnapper)
#' ssm <- MLmulti(PRSnapper, ncp = 1, model = "SSM")
#' msm1 <- MLmulti(PRSnapper, ncp = 1, model = "MSM1")
#' msm2 <- MLmulti(PRSnapper, ncp = 1, model = "MSM2")
#' msm3 <- MLmulti(PRSnapper, ncp = 1, model = "MSM3")
#'
#' compare_models(list(ssm, msm1, msm2, msm3))
#'
#' @export
compare_models <- function(MLZ_model.list, figure = TRUE, color = NULL) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(list = old_par), add = TRUE)

  model.names <- vapply(MLZ_model.list, getElement, c("x"), "Model")
  model <- unique(model.names)
  if(length(model) > 1) stop("More than one model identified in MLZ_model.list")
  length.units <- vapply(MLZ_model.list, getElement, c("x"), "length.units")
  length.units <- unique(length.units)
  if(length(length.units) != 0) length.units <- paste0("(", length.units, ")") else length.units <- NULL
  
  nspec <- vapply(MLZ_model.list, getElement, integer(1), "n.species")
  nspec <- unique(nspec)
  if(length(nspec) > 1) stop("Different number of species among models.")
  ncp <- vapply(MLZ_model.list, getElement, integer(1), "n.changepoint")
  ncp.text <- paste0(ncp, "-change point")
  negLL <- vapply(MLZ_model.list, getElement, numeric(1), "negLL")

  if(model == "ML") {
    npar = 2 * (ncp + 1)
    AIC = 2 * (negLL + npar)
    delta.AIC = AIC - min(AIC)
    output <- matrix(c(negLL, npar, AIC, delta.AIC), nrow = length(MLZ_model.list))
    dimnames(output) <- list(ncp.text, c("negLL", "npar", "AIC", "delta.AIC"))
    if(figure) {
      par(las = 1)
      if(is.null(color)) color <- rich.colors(length(MLZ_model.list))
      plot(MeanLength ~ Year, MLZ_model.list[[1]]@time.series, ylab = paste("Mean Length (> Lc)", length.units), 
           typ = "o", pch = 16)
      for(i in 1:length(MLZ_model.list)) {
        lines(Predicted ~ Year, MLZ_model.list[[i]]@time.series, col = color[i], lwd = 2)
      }
      legend("topright", legend = ncp.text, lwd = 2, col = color, bty = "n")
    }
  }
  if(model == "MLCR") {
    npar = 2 * (ncp + 1) + 2
    AIC = 2 * (negLL + npar)
    delta.AIC = AIC - min(AIC)
    output <- matrix(c(negLL, npar, AIC, delta.AIC), nrow = length(MLZ_model.list))
    dimnames(output) <- list(ncp.text, c("negLL", "npar", "AIC", "delta.AIC"))
    if(figure) {
      if(is.null(color)) color <- rich.colors(length(MLZ_model.list))
      par(mfrow = c(1,2), mar = c(5,4,1,1), las = 1)
      plot(MeanLength ~ Year, MLZ_model.list[[1]]@time.series, ylab = paste("Mean Length (> Lc)", length.units),
           typ = "o", pch = 16)
      for(i in 1:length(MLZ_model.list)) {
        lines(Predicted.ML ~ Year, MLZ_model.list[[i]]@time.series, col = color[i], lwd = 2)
      }
      legend("topright", legend = ncp.text, lwd = 2, col = color, bty = "n")

      plot(CPUE ~ Year, MLZ_model.list[[1]]@time.series, typ = "o", pch = 16)
      for(i in 1:length(MLZ_model.list)) {
        lines(Predicted.CPUE ~ Year, MLZ_model.list[[i]]@time.series, col = color[i], lwd = 2)
      }
      legend("topright", legend = ncp.text, lwd = 2, col = color, bty = "n")
    }
  }
  if(model == "MLmulti") {
    msm.model <- vapply(MLZ_model.list, getElement, c("x"), "multimodel")

    ssm.ind <- grep("SSM", msm.model)
    msm1.ind <- grep("MSM1", msm.model)
    msm2.ind <- grep("MSM2", msm.model)
    msm3.ind <- grep("MSM3", msm.model)

    npar <- numeric(length(MLZ_model.list))
    npar[ssm.ind] <- 2 * nspec * (ncp[ssm.ind] + 1)
    npar[msm1.ind] <- nspec + ncp[msm1.ind] + nspec * (ncp[msm1.ind] + 1)
    npar[msm2.ind] <- 3 * nspec + 2 * ncp[msm2.ind] - 1
    npar[msm3.ind] <- 2 * (nspec + ncp[msm3.ind])

    msm.text <- vector(length = length(MLZ_model.list))
    msm.text[ssm.ind] <- "Single Species Model"
    msm.text[msm1.ind] <- "Multispecies Model 1"
    msm.text[msm2.ind] <- "Multispecies Model 2"
    msm.text[msm3.ind] <- "Multispecies Model 3"

    AIC = 2 * (negLL + npar)
    delta.AIC = round(AIC - min(AIC), 2)

    output <- matrix(c(ncp, negLL, npar, AIC, delta.AIC), nrow = length(MLZ_model.list))
    rownames(output) <- msm.text
    colnames(output) <- c("n.changepoint", "negLL", "npar", "AIC", "delta.AIC")

    if(figure) {
      par(mfrow = c(2,2), las = 1)
      if(is.null(color)) color <- rich.colors(length(MLZ_model.list))
      nyrs <- nrow(MLZ_model.list[[1]]@time.series)/nspec

      for(i in 1:nspec) {
        plot(MeanLength ~ Year, MLZ_model.list[[1]]@time.series[((i-1)*nyrs+1):(i*nyrs), ],
             ylab = paste("Mean Length (> Lc)", length.units), typ = "o", pch = 16)
        for(j in 1:length(MLZ_model.list)) {
          lines(Predicted ~ Year, MLZ_model.list[[j]]@time.series[((i-1)*nyrs+1):(i*nyrs), ],
                col = color[j], lwd = 2)
        }
        title(MLZ_model.list[[1]]@Stock[i])
      }
      plot.new()
      legend("topright", legend = msm.text, lwd = 2, lty = 2, col = color, bty = "n")
    }
  }
  return(output[order(delta.AIC), ])
}
