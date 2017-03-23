#' Mean length-based mortality estimator
#'
#' Estimator of instantaneous total mortality (Z) from a time series of mean length data.
#'
#' @param MLZ_data An object of class \code{\linkS4class{MLZ_data}} containing mean lengths and
#' life history data of stock.
#' @param ncp The number of change points in total mortality in the time series. \code{ncp + 1} total
#' mortality rates will be estimated.
#' @param start An optional list of starting values. See details.
#' @param spawn Whether 'continuous' or 'annual' (pulse) spawning is modeled.
#' @param grid.search If \code{TRUE}, a grid search will be performed using the \code{\link{profile_ML}}
#' function to find the best starting values for the change points (the years when mortality changes).
#' Ignored if \code{ncp = 0} or if \code{start} is provided. 
#' @param parallel Whether grid search is performed with parallel processing. Ignored if \code{grid.search = FALSE}.
#' @param min.time The minimum number of years between each change point for the grid search, passed
#' to \code{\link{profile_ML}}. Not used if \code{grid.search = FALSE}.
#' @param figure If \code{TRUE}, a call to \code{plot} of observed and predicted mean lengths will be produced.
#' @details For a model with \code{I} change points, the starting values in
#' \code{start} is a list with the following entries:
#' Z a vector of \code{length = I+1}.
#' yearZ a vector of \code{length = I}.
#'
#' \code{start} can be \code{NULL}, in which case, the supplied starting values depend on
#' the value of \code{grid.search}. If \code{grid.search = TRUE}, starting values will use the
#' values for \code{yearZ} which minimize the negative log-likelihood from the grid search.
#' Otherwise, the starting values for \code{yearZ} evenly divide the time series.
#'
#' @references Gedamke, T. and Hoenig, J.M. 2006. Estimating mortality from mean length data in
#' nonequilibrium situations, with application to the assessment of goosefish.
#' Transactions of the American Fisheries Society 135:476-487.
#' @return An object of class \code{\linkS4class{MLZ_model}}.
#' @examples
#' data(Goosefish)
#' ML(Goosefish, ncp = 2)
#' ML(Goosefish, ncp = 2, start = list(Z = c(0.1, 0.3, 0.5), yearZ = c(1978, 1988)))
#' ML(Goosefish, ncp = 2, grid.search = TRUE)
#'
#' @seealso \code{\link{profile_ML}}
#' @export
ML <- function(MLZ_data, ncp, start = NULL, spawn = c("continuous", "annual"),
               grid.search = TRUE, parallel = ifelse(ncp > 2, TRUE, FALSE), min.time = 3, figure = TRUE) {
  if(!inherits(MLZ_data, "MLZ_data")) stop("No object of class 'MLZ_data' found.")
  ncp <- as.integer(ncp)
  parallel <- as.logical(parallel)
  spawn <- match.arg(spawn)
  if(spawn == "continuous") spawn.cont <- 1L else spawn.cont <- 0L
  
  years <- MLZ_data@Year
  max.years <- data.frame(Year = min(years):max(years))
  data.df <- data.frame(Year = MLZ_data@Year, MeanLength = MLZ_data@MeanLength, ss = MLZ_data@ss)
  full <- left_join(max.years, data.df, by = "Year")

  tmb.dat <- list(LH = c(MLZ_data@vbLinf, MLZ_data@vbK), Lc = MLZ_data@Lc, nbreaks = ncp,
                  Lbar = full$MeanLength, ss = full$ss)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss <= 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0 | tmb.dat$ss <= 0] <- 0

  if(ncp == 0) {
    if(!is.null(start)) {
      if(!"Z" %in% names(start)) stop("Error in start list. Entry of name 'Z' not found.")
      if(length(start$Z) != 1)
        stop("Entry in name 'Z' of start list is not a numeric of length 1.")
    }
    else {
      start <- list(Z = MLZ_data@vbK * (MLZ_data@vbLinf - MLZ_data@MeanLength[1]) /
                      (MLZ_data@MeanLength[1] - MLZ_data@Lc))
      start$Z[is.na(start$Z) | start$Z <= 0 | is.infinite(start$Z)] <- 0.5
      start$Z[start$Z > 1] <- 1
    }
    opt <- optim(start$Z, MLeqnegLL, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                 LH = tmb.dat$LH, Lc = tmb.dat$Lc, spCont = spawn.cont, 
                 method = "BFGS", control = list(maxit = 1e7))
    data.pred <- MLeqpred(opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, LH = tmb.dat$LH, 
                          Lc = tmb.dat$Lc, spCont = spawn.cont)
    opt$Lpred <- rep(data.pred[[1]], length(tmb.dat$Lbar))
    sigmaL <- data.pred[[2]]
    opt$par <- c(opt$par, sigmaL)
    MLhessian <- hessian(MLeqfullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, spCont = spawn.cont)
    covariance <- solve(MLhessian)
    opt$corr <- cov2cor(covariance)
    opt$gradient <- grad(MLeqfullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, spCont = spawn.cont)

    results.matrix <- matrix(c(opt$par, sqrt(diag(covariance))), ncol = 2)
    rownames(results.matrix) <- rownames(opt$corr) <- colnames(opt$corr) <- names(opt$par) <- names(opt$gradient) <- c("Z", "sigma")
  }
  if(ncp > 0) {
    if(!is.null(start)) {
      if(!"Z" %in% names(start) || !"yearZ" %in% names(start)) stop("Error in setup of start list. See help file.")
      if(length(start$Z) != ncp + 1) stop("Length of starting Z vector not equal to ncp + 1.")
      if(length(start$yearZ) != ncp) stop("Length of starting yearZ vector not equal to ncp.")
      start$yearZ <- start$yearZ - MLZ_data@Year[1] + 1
    }
    else {
      if(grid.search) {
        grid.output <- profile_ML(MLZ_data, ncp, min.time = min.time, parallel = parallel, spawn = spawn, figure = FALSE)
        index.min <- which.min(grid.output$negLL)
        styearZ <- as.numeric(grid.output[index.min, 1:ncp])
        styearZ <- styearZ - MLZ_data@Year[1] + 1
      }
      else {
        styearZ <- length(tmb.dat$Lbar) * (1:ncp) / (ncp+1)
      }
      stZ <- MLZ_data@vbK * (MLZ_data@vbLinf - MLZ_data@MeanLength[c(1,styearZ)]) /
        (MLZ_data@MeanLength[c(1,styearZ)] - MLZ_data@Lc)
      stZ[is.na(stZ) | stZ <= 0 | is.infinite(stZ)] <- 0.5
      stZ[stZ > 1] <- 1
      start <- list(Z = stZ, yearZ = styearZ)
    }
    opt <- optim(c(start$Z, start$yearZ), MLnegLL, Lbar = tmb.dat$Lbar,
                 ss = tmb.dat$ss, LH = tmb.dat$LH, Lc = tmb.dat$Lc,
                 nbreaks = tmb.dat$nbreaks, spCont = spawn.cont, 
                 method = "BFGS", control = list(maxit = 1e7))
    data.pred <- MLpred(opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, LH = tmb.dat$LH, 
                        Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks, spCont = spawn.cont)
    opt$Lpred <- data.pred[[1]]
    sigmaL <- data.pred[[2]]
    opt$par <- c(opt$par, sigmaL)
    MLhessian <- hessian(MLfullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks, 
                         spCont = spawn.cont)
    covariance <- solve(MLhessian)
    opt$corr <- cov2cor(covariance)
    opt$gradient <- grad(MLfullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks, 
                         spCont = spawn.cont)

    results.matrix <- matrix(c(opt$par, sqrt(diag(covariance))), ncol = 2)
    Z.name <- paste0("Z[", 1:(ncp+1), "]")
    yearZ.name <- paste0("yearZ[", 1:ncp, "]")
    rownames(results.matrix) <- rownames(opt$corr) <- colnames(opt$corr) <- names(opt$par) <- names(opt$gradient) <- c(Z.name, yearZ.name, "sigma")
    year.ind <- grep("yearZ", rownames(results.matrix))
    results.matrix[year.ind, 1] <- results.matrix[year.ind, 1] + MLZ_data@Year[1] - 1
  }
  colnames(results.matrix) <- c("Estimate", "Std. Error")
  time.series <- full
  time.series$Predicted <- opt$Lpred
  time.series$Residual <- time.series$MeanLength - time.series$Predicted
  opt$message <- c(opt$message, paste0("Spawning is assumed to be ", spawn, " in model."))
  MLZ_model <- new("MLZ_model", Stock = MLZ_data@Stock, Model = "ML", time.series = time.series,
                   estimates = results.matrix, negLL = opt$value, n.changepoint = ncp, n.species = 1L,
                   opt = opt, length.units = MLZ_data@length.units)
  attr(MLZ_model, "spawn") <- spawn
  if(exists("grid.output")) MLZ_model@grid.search <- grid.output
  if(figure) plot(MLZ_model)
  return(MLZ_model)
}

#' Mean length with catch rate mortality estimator
#'
#' Estimator of instantaneous total mortality (Z) from a time series of mean length data.
#'
#' @param MLZ_data An object of class \code{\linkS4class{MLZ_data}} containing mean lengths and
#' life history data of stock.
#' @param ncp The number of change points in total mortality in the time series. \code{ncp + 1} total
#' mortality rates will be estimated.
#' @param CPUE.type Indicates whether CPUE time series is abundance or biomass based.
#' @param loglikeCPUE Indicates whether the log-likelihood for the CPUE will be lognormally or
#' normally distributed.
#' @param start An optional list of starting values. See details.
#' @param spawn Whether 'continuous' or 'annual' (pulse) spawning is modeled.
#' @param grid.search If \code{TRUE}, a grid search will be performed using the \code{\link{profile_MLCR}}
#' function to find the best starting values for the change points (the years when mortality changes).
#' Ignored if \code{ncp = 0} or if \code{start} is provided.
#' @param parallel Whether grid search is performed with parallel processing. Ignored if \code{grid.search = FALSE}.
#' @param min.time The minimum number of years between each change point for the grid search, passed
#' to \code{\link{profile_MLCR}}. Not used if \code{grid.search = FALSE}.
#' @param figure If \code{TRUE}, a call to \code{plot} of observed and predicted mean lengths will be produced.
#'
#' @references Huynh, Q.C., Gedamke, T., Porch, C.E., Hoenig, J.M., Walter, J.F, Bryan, M., and
#' Brodziak, J. In revision. Estimating Total Mortality Rates from Mean Lengths and
#' Catch Rates in Non-equilibrium Situations. Transactions of the American Fisheries Society.
#'
#' @return An object of class \code{\linkS4class{MLZ_model}}.
#'
#' @details For a model with \code{I} change points, the starting values in
#' \code{start} is a list with the following entries:
#' Z a vector of \code{length = I+1}.
#' yearZ a vector of \code{length = I}.
#'
#' \code{start} can be \code{NULL}, in which case, the supplied starting values depend on
#' the value of \code{grid.search}. If \code{grid.search = TRUE}, starting values will use the
#' values for \code{yearZ} which minimize the negative log-likelihood from the grid search.
#' Otherwise, the starting values for \code{yearZ} evenly divide the time series.
#'
#' @examples
#' data(MuttonSnapper)
#' MLCR(MuttonSnapper, ncp = 2, CPUE.type = "WPUE", grid.search = TRUE)
#'
#' @seealso \code{\link{profile_MLCR}}
#'
#' @export
MLCR <- function(MLZ_data, ncp, CPUE.type = c(NA, "WPUE", "NPUE"), loglikeCPUE = c("lognormal", "normal"),
                 spawn = c("continuous", "annual"), start = NULL, grid.search = TRUE, 
                 parallel = ifelse(ncp > 2, TRUE, FALSE), min.time = 3, figure = TRUE) {

  if(!inherits(MLZ_data, "MLZ_data")) stop("No object of class 'MLZ_data' found.")
  if(ncp == 0L) stop("Zero changepoints is not currently supported. Use ML().")
  parallel <- as.logical(parallel)
  
  CPUE.type <- match.arg(CPUE.type)
  if(is.na(CPUE.type)) stop("Argument CPUE.type must be identified in either: weight 'WPUE' or abundance 'NPUE'.")
  loglikeCPUE <- match.arg(loglikeCPUE)
  
  spawn <- match.arg(spawn)
  if(CPUE.type == "WPUE" & spawn == "annual") {
    message("Annual spawning currently is not supported with WPUE. Continuous spawning assumed.")
    spawn <- "continuous"
  }
  if(spawn == "continuous") spawn.cont <- 1L else spawn.cont <- 0L
  
  years <- MLZ_data@Year
  max.years <- data.frame(Year = min(years):max(years))
  data.df <- data.frame(Year = MLZ_data@Year, MeanLength = MLZ_data@MeanLength, ss = MLZ_data@ss,
                        CPUE = MLZ_data@CPUE)
  full <- left_join(max.years, data.df, by = "Year")

  ncp <- as.integer(ncp)
  ll.int <- ifelse(loglikeCPUE == "lognormal", 0L, 1L)
  tmb.dat <- list(LH = c(MLZ_data@vbLinf, MLZ_data@vbK), Lc = MLZ_data@Lc, nbreaks = ncp,
                  Lbar = full$MeanLength, ss = full$ss, CPUE = full$CPUE, loglikeCPUE = ll.int)
  if(CPUE.type == "WPUE") tmb.dat$LH <- c(tmb.dat$LH, MLZ_data@lwb)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss <= 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0 | tmb.dat$ss <= 0] <- 0
  tmb.dat$CPUE[is.na(tmb.dat$CPUE) | tmb.dat$CPUE < 0] <- -99

  if(!is.null(start)) {
    if(!"Z" %in% names(start) || !"yearZ" %in% names(start)) stop("Error in setup of start list. See help file.")
    if(length(start$Z) != ncp + 1) stop("Length of starting Z vector not equal to ncp + 1.")
    if(length(start$yearZ) != ncp) stop("Length of starting yearZ vector not equal to ncp.")
    start$yearZ <- start$yearZ - MLZ_data@Year[1] + 1
  }
  else {
    if(grid.search) {
      grid.output <- profile_MLCR(MLZ_data, ncp, CPUE.type, loglikeCPUE, spawn = spawn, 
                                  min.time = min.time, parallel = parallel, figure = FALSE)
      index.min <- which.min(grid.output$negLL)
      styearZ <- as.numeric(grid.output[index.min, 1:ncp])
      styearZ <- styearZ - MLZ_data@Year[1] + 1
    }
    else {
      styearZ <- length(tmb.dat$Lbar) * (1:ncp) / (ncp+1)
    }
    stZ <- MLZ_data@vbK * (MLZ_data@vbLinf - MLZ_data@MeanLength[c(1,styearZ)]) /
      (MLZ_data@MeanLength[c(1,styearZ)] - MLZ_data@Lc)
    stZ[is.na(stZ) | stZ <= 0] <- 0.5
    stZ[stZ > 1] <- 1
    start <- list(Z = stZ, yearZ = styearZ)
  }

  fx <- get(paste0("ML", CPUE.type, "negLL"))
  gx <- get(paste0("ML", CPUE.type, "pred"))
  hx <- get(paste0("ML", CPUE.type, "fullnegLL"))

  opt <- optim(c(start$Z, start$yearZ), fx, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
               CPUE = tmb.dat$CPUE, LH = tmb.dat$LH, Lc = tmb.dat$Lc,
               nbreaks = tmb.dat$nbreaks, loglikeCPUE = tmb.dat$loglikeCPUE, spCont = spawn.cont,
               method = "BFGS", control = list(maxit = 1e7))
  data.pred <- gx(opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, CPUE = tmb.dat$CPUE,
                  LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks, 
                  loglikeCPUE = tmb.dat$loglikeCPUE, spCont = spawn.cont)
  opt$Lpred <- data.pred[[1]]
  opt$Ipred <- data.pred[[2]]
  q <- data.pred[[3]]
  sigmav <- data.pred[[4]]
  opt$par <- c(opt$par, q, sigmav)
  MLCRhessian <- hessian(hx, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, CPUE = tmb.dat$CPUE, 
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks, 
                         loglikeCPUE = tmb.dat$loglikeCPUE, spCont = spawn.cont)
  covariance <- solve(MLCRhessian)
  opt$corr <- cov2cor(covariance)
  opt$gradient <- grad(hx, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, CPUE = tmb.dat$CPUE, 
                       LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks, 
                       loglikeCPUE = tmb.dat$loglikeCPUE, spCont = spawn.cont)

  results.matrix <- matrix(c(opt$par, sqrt(diag(covariance))), ncol = 2)
  Z.name <- paste0("Z[", 1:(ncp+1), "]")
  yearZ.name <- paste0("yearZ[", 1:ncp, "]")
  rownames(results.matrix) <- rownames(opt$corr) <- colnames(opt$corr) <- names(opt$par) <- names(opt$gradient) <- c(Z.name, yearZ.name, "q", "sigmaL", "sigmaI")
  colnames(results.matrix) <- c("Estimate", "Std. Error")
  year.ind <- grep("yearZ", rownames(results.matrix))
  results.matrix[year.ind, 1] <- results.matrix[year.ind, 1] + MLZ_data@Year[1] - 1

  time.series <- data.frame(Predicted.ML = opt$Lpred, Predicted.CPUE = opt$Ipred)
  time.series <- cbind(full, time.series)
  time.series$Residual.ML <- time.series$MeanLength - time.series$Predicted.ML
  time.series$Residual.CPUE <- time.series$CPUE - time.series$Predicted.CPUE

  opt$message <- c(opt$message, paste0("Catch rate assumed to be ", CPUE.type), 
                   paste0("Catch rate likelihood used ", loglikeCPUE, " distribution."),
                   paste0("Spawning is assumed to be ", spawn, " in model."))
  
  MLZ_model <- new("MLZ_model", Stock = MLZ_data@Stock, Model = "MLCR", time.series = time.series,
                   estimates = results.matrix, negLL = opt$value, n.changepoint = ncp, n.species = 1L,
                   opt = opt, length.units = MLZ_data@length.units)
  attr(MLZ_model, "spawn") <- spawn
  if(exists("grid.output")) MLZ_model@grid.search <- grid.output
  if(figure) plot(MLZ_model)
  return(MLZ_model)
}

#' Multispecies mean length mortality estimator
#'
#' Estimator of instantaneous total mortality (Z) from a time series of mean length data for a suite of
#' stocks that are fished together.
#'
#' @param MLZ.list A list containing objects of class \code{\linkS4class{MLZ_data}}.
#' @param ncp The number of change points in total mortality in the time series. \code{ncp + 1} total
#' mortality rates will be estimated.
#' @param model The multispecies model to be used.
#' @param start An optional list of starting values. See details.
#' @param spawn Whether 'continuous' or 'annual' (pulse) spawning is modeled.
#' @param grid.search If \code{TRUE}, a grid search will be performed using the \code{\link{profile_MLmulti}}
#' function to find the best starting values for the change points (the years when mortality changes).
#' Ignored if \code{start} is provided.
#' @param parallel Whether grid search is performed in parallel. Ignored if \code{grid.search = FALSE}.
#' @param min.time The minimum number of years between each change point for the grid search, passed
#' to \code{\link{profile_MLmulti}}. Not used if \code{grid.search = FALSE}.
#' @param figure If \code{TRUE}, a call to \code{plot} of observed and predicted mean lengths will be produced.
#'
#' @return An object of class \code{\linkS4class{MLZ_model}}.
#'
#' @details For a model with \code{I} change points and \code{N} species, the starting values in
#' \code{start} is a list with the following entries:
#'
#' Single Species Model (SSM):
#' \tabular{ll}{
#' \code{Z} \tab a matrix with \code{nrow = I+1} and \code{ncol = N}.\cr
#' \code{yearZ} \tab a matrix with \code{nrow = I} and \code{ncol = N}.\cr
#' }
#'
#' Multispecies Model 1 (MSM1):
#' \tabular{ll}{
#' \code{Z} \tab a matrix with \code{nrow = I+1} and \code{ncol = N}.\cr
#' \code{yearZ} \tab a vector with \code{length = I}.\cr
#' }
#'
#' Multispecies Model 2 (MSM2):
#' \tabular{ll}{
#' \code{Z1} \tab a vector with \code{length = N}.\cr
#' \code{yearZ} \tab a vector with \code{length = I}.\cr
#' \code{delta} \tab a vector with \code{length = I}.\cr
#' \code{epsilon} \tab a vector with \code{length = N-1}.\cr
#' }
#'
#' Multispecies Model 3 (MSM3):
#' \tabular{ll}{
#' \code{Z1} \tab a vector with \code{length = N}.\cr
#' \code{yearZ} \tab a vector with \code{length = I}.\cr
#' \code{delta} \tab a vector with \code{length = I}.\cr
#' }
#'
#' \code{start} can be \code{NULL}, in which case, the supplied starting values depend on
#' the value of \code{grid.search}. If \code{grid.search = TRUE}, starting values will use the
#' values for \code{yearZ} which minimize the negative log-likelihood from the grid search.
#' Otherwise, the starting values for \code{yearZ} evenly divide the time series.
#'
#' @references Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. 2017. Multispecies Extensions
#' to a Nonequilibrium Length-Based Mortality Estimator. Marine and Coastal Fisheries 9:68-78.
#'
#' @seealso \code{\link{profile_MLmulti}}
#'
#' @examples
#' data(PRSnapper)
#' MLmulti(PRSnapper, ncp = 1, model = "SSM")
#'
#' MSM1.start.Z <- matrix(0.5, nrow = 3, ncol = 2)
#' MSM1.start.yearZ <- 1990
#' MLmulti(PRSnapper, ncp = 1, model = "MSM1", start = list(Z = MSM1.start.Z, yearZ = MSM1.start.yearZ),
#' grid.search = FALSE)
#'
#' MLmulti(PRSnapper, ncp = 1, model = "MSM2")
#'
#' st.Z1 <- rep(0.5, 3)
#' st.yearZ <- 1990
#' st.delta <- 1
#' MLmulti(PRSnapper, ncp = 1, model = "MSM3", start = list(Z1 = st.Z1, yearZ = st.yearZ, delta = st.delta))
#' @export
MLmulti <- function(MLZ.list, ncp, model = c("SSM", "MSM1", "MSM2", "MSM3"), start = NULL,
                    spawn = c("continuous", "annual"), grid.search = TRUE, parallel = ifelse(ncp > 2, TRUE, FALSE), 
                    min.time = 3, figure = TRUE) {

  model <- match.arg(model)
  parallel <- as.logical(parallel)
  if(!is.list(MLZ.list)) stop("No list found.")
  class.test <- vapply(MLZ.list, inherits, TRUE, "MLZ_data")
  if(!all(class.test)) stop("Not all entries in list are of class 'MLZ_data'")
  length.units <- vapply(MLZ.list, getElement, c("x"), "length.units")
  length.units <- unique(length.units)
  if(length(length.units) > 1) stop("Different length units among species. All lengths should be the same units, e.g. cm.")
  
  ncp <- as.integer(ncp)
  if(ncp == 0) stop("Zero change points not supported. Use ML().")
  spawn <- match.arg(spawn)
  if(spawn == "continuous") spawn.cont <- 1L else spawn.cont <- 0L
  
  nspec <- length(MLZ.list)
  years <- lapply(MLZ.list, getElement, "Year")
  years <- unique(do.call(c, years))

  max.years <- data.frame(Year = min(years):max(years))
  Lbar.df <- lapply(MLZ.list, timeseries_multi)
  Lbar.df <- do.call(rbind, Lbar.df)
  full <- left_join(max.years, Lbar.df, by = "Year")
  full$Sp <- rep(1:nspec, nrow(max.years))
  Lbar <- acast(full, Sp ~ Year, value.var = "MeanLength")
  ss <- acast(full, Sp ~ Year, value.var = "ss")

  Linf <- vapply(MLZ.list, getElement, numeric(1), "vbLinf")
  K <- vapply(MLZ.list, getElement, numeric(1), "vbK")
  if("MSM2" %in% model || "MSM3" %in% model) {
    M <- vapply(MLZ.list, getElement, numeric(1), "M")
    LH <- matrix(c(Linf, K, M), ncol = 3)
  } else LH <- matrix(c(Linf, K), ncol = 2)
  Lc <- vapply(MLZ.list, getElement, numeric(1), "Lc")

  tmb.dat <- list(LH = LH, Lc = Lc, nbreaks = ncp, Lbar = Lbar, ss = ss, nspec = nspec)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss == 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0] <- 0

  if(!is.null(start)) {
    if("SSM" %in% model) {
      if(!"Z" %in% names(start) || !"yearZ" %in% names(start)) stop("Error in names of start list.")
      if(nrow(start$Z) != ncp+1 || ncol(start$Z) != nspec || nrow(start$yearZ) != ncp ||
         ncol(start$yearZ) != nspec) stop("Error in dimension of entries of start list.")
    }
    if("MSM1" %in% model) {
      if(!"Z" %in% names(start) || !"yearZ" %in% names(start)) stop("Error in names of start list.")
      if(nrow(start$Z) != ncp+1 || ncol(start$Z) != nspec ||
         length(start$yearZ) != ncp) stop("Error in dimension of entries of start list.")
    }
    if("MSM2" %in% model) {
      if(!"Z1" %in% names(start) || !"yearZ" %in% names(start) || !"delta" %in% names(start) ||
         !"epsilon" %in% names(start)) stop("Error in names of start list.")
      if(length(start$Z1) != nspec || length(start$yearZ) != ncp || length(start$delta) != ncp ||
         length(start$epsilon) != nspec-1) stop("Error in dimension of entries of start list.")
    }
    if("MSM3" %in% model) {
      if(!"Z1" %in% names(start) || !"yearZ" %in% names(start) ||
         !"delta" %in% names(start)) stop("Error in names of start list.")
      if(length(start$Z1) != nspec || length(start$yearZ) != ncp ||
         length(start$delta) != ncp) stop("Error in dimension of entries of start list.")
    }
    start$yearZ <- start$yearZ - max.years[1, 1] + 1
  } else {
    if(grid.search) {
      year.range <- min(years):max(years) - min(years) + 1
      year.range <- year.range[2:(length(year.range) - 1)]
      ygrid <- list()
      for(i in 1:ncp) ygrid[[i]] <- year.range
      ygrid <- expand.grid(ygrid)
      if(ncp > 1) {
        for(i in 1:(ncp-1)) {
          truncate.index <- (ygrid[, i] < ygrid[, i+1]) & (ygrid[, i+1] - ygrid[, i] >= min.time)
          ygrid <- ygrid[truncate.index, ]
        }
      }
      if("SSM" %in% model || "MSM1" %in% model) pr <- .SSM_MSM1_profile(tmb.dat, 0.5, parallel, ygrid, spawn.cont)
      if("MSM2" %in% model || "MSM3" %in% model) pr <- .MSM23_profile(tmb.dat, 0.5, parallel, ygrid, model, spawn.cont)

      if(model == "SSM") {
        index.min <- apply(pr$loglikesp, 2, which.min)
        styearZ <- t(as.matrix(ygrid[index.min, ]))
      } else {
        index.min <- which.min(pr$negLL)
        styearZ <- as.numeric(ygrid[index.min, ])
      }
    }
    else {
      styearZ <- nrow(max.years) * (1:ncp) / (ncp+1)
    }
    if("SSM" %in% model || "MSM1" %in% model) {
      stZ <- K * (Linf - tmb.dat$Lbar[, styearZ]) / (tmb.dat$Lbar[, styearZ] - Lc)
      stZ[is.na(stZ) | stZ <= 0] <- 0.5
      stZ[stZ > 1] <- 1
    }
    if("MSM2" %in% model || "MSM3" %in% model) {
      stZ1 <- K * (Linf - tmb.dat$Lbar[, 1]) / (tmb.dat$Lbar[, 1] - Lc)
      stZ1[is.na(stZ1) | stZ1 <= 0] <- 0.5
      stZ[stZ > 1] <- 1
    }
    if(model == "SSM") {
      stZ.new <- matrix(NA, nrow = nspec, ncol = ncp)
      for(i in 1:ncp) stZ.new[, i] <- diag(stZ[, ((i-1)*nspec+1):(i*nspec)])
      stZ <- cbind(stZ1, stZ.new)
      start <- list(Z = t(stZ), yearZ = styearZ)
    }
    if(model == "MSM1") {
      stZ <- cbind(stZ1, stZ)
      start <- list(Z = t(stZ), yearZ = styearZ)
    }
    if(model == "MSM2") start <- list(Z1 = stZ1, delta = rep(1, ncp), epsilon = rep(1, nspec - 1), yearZ = styearZ)
    if(model == "MSM3") start <- list(Z1 = stZ1, delta = rep(1, ncp), yearZ = styearZ)
  }

  if("SSM" %in% model || "MSM1" %in% model) est <- .SSM_MSM1_est(tmb.dat, start, model, spawn.cont)
  if("MSM2" %in% model || "MSM3" %in% model) est <- .MSM23_est(tmb.dat, start, model, spawn.cont)

  results.matrix <- est$results.matrix
  year.ind <- grep("yearZ", rownames(results.matrix))
  results.matrix[year.ind, 1] <- results.matrix[year.ind, 1] + min(years) - 1

  opt <- est$opt
  Lbar.df$Predicted <- as.numeric(t(opt$Lpred))
  Lbar.df$Residual <- Lbar.df$MeanLength - Lbar.df$Predicted

  sp.name <- lapply(MLZ.list, getElement, "Stock")
  sp.ind <- vapply(sp.name, length, numeric(1))
  num <- which(sp.ind == 0L)
  sp.name[num] <- paste0("Species_", num)
  Lbar.df$Stock <- rep(do.call(c, sp.name), each = length(years))
  
  opt$message <- c(opt$message, paste0("Spawning is assumed to be ", spawn, " in model."))

  MLZ_model <- new("MLZ_model", Stock = do.call(c, sp.name), Model = "MLmulti",
                   time.series = Lbar.df, estimates = results.matrix,
                   negLL = opt$value, n.changepoint = ncp, n.species = nspec,
                   opt = opt, length.units = length.units)
  attr(MLZ_model, "multimodel") <- model
  attr(MLZ_model, "spawn") <- spawn
  if(grid.search) {
    grid.output <- as.data.frame(ygrid + min(years) - 1)
    grid.output <- cbind(grid.output, pr$negLL, pr$loglikesp)
    names(grid.output) <- c(paste0("Year", 1:ncp), "negLL", paste0("Species_", 1:nspec))
    rownames(grid.output) <- 1:nrow(grid.output)
    MLZ_model@grid.search <- grid.output
  }

  if(figure) plot(MLZ_model)
  return(MLZ_model)
}

#' Mean length with effort mortality estimator
#'
#' Estimator of fishing and natural mortality from a time series of mean length and effort data.
#'
#' @param MLZ_data An object of class \code{\linkS4class{MLZ_data}} containing mean lengths and
#' life history data of stock.
#' @param start A list of starting values. Names of start list must contain \code{q} and \code{M}.
#' @param n_age The number of ages above age tc in the model.
#' @param estimate.M If \code{TRUE}, natural mortality (M) will be estimated. Otherwise, the value of M
#' is obtained from slot \code{MLZ_data@M}.
#' @param log.par Whether parameters are estimated in logspace (\code{TRUE}) or untransformed space (\code{FALSE}).
#' @param eff_init The assumed equilibrium effort prior to the first year of the model (0 = virgin conditions).
#' @param n_season The number of seasons modeled in a year.
#' @param obs_season The season corresponding to the observed mean lengths.
#' @param timing The fraction of time (i.e., between 0 - 1) within \code{obs_season} that mean lengths are observed.
#' @param figure If \code{TRUE}, a call to \code{plot} of observed and predicted mean lengths will be produced.
#'
#' @references Then, A.Y, Hoenig, J.M, and Huynh, Q.C. In revision. Estimating fishing and natural
#' mortality rates, and catchability coefficient, from a series of observations on mean length and
#' fishing effort. ICES Journal of Marine Science.
#'
#' @return An object of class \code{\linkS4class{MLZ_model}}.
#'
#' @examples
#' data(Nephrops)
#' Nephrops <- calc_ML(Nephrops, sample.size = FALSE)
#' MLeffort(Nephrops, start = list(q = 0.1, M = 0.2),
#'          n_age = 24, eff_init = Nephrops@Effort[1])
#' @export
MLeffort <- function(MLZ_data, start, n_age, estimate.M = TRUE, log.par = FALSE,
                     eff_init = 0, n_season = 1L, obs_season = 1L, timing = 0, figure = TRUE) {

  if(!inherits(MLZ_data, "MLZ_data")) stop("No object of class 'MLZ_data' found.")
  if(!"q" %in% names(start)) stop("Starting value for q not found.")
  if(!"M" %in% names(start) & estimate.M) stop("Starting value for M not found.")
  if(length(MLZ_data@vbt0) == 0) stop("vbt0 is missing in MLZ_data.")
  if(!estimate.M & length(MLZ_data@M) == 0) stop("M is missing in MLZ_data.")

  n_age <- as.integer(n_age)
  n_season <- as.integer(n_season)
  obs_season <- as.integer(obs_season)
  if(obs_season > n_season) stop("obs_season is greater than n_season.")
  if(timing < 0 || timing > 1) stop("timing must be between 0 or 1.")

  tmb.dat <- list(LH = c(MLZ_data@vbLinf, MLZ_data@vbK, MLZ_data@vbt0, MLZ_data@M),
                  Lc = MLZ_data@Lc, Lbar = MLZ_data@MeanLength, ss = MLZ_data@ss,
                  Effort = MLZ_data@Effort, n_age = n_age, eff_init = eff_init)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss == 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0] <- 0
  if(log.par) {
    tmb.dat$logpar <- 1L
    start <- lapply(start, log)
  } else tmb.dat$logpar <- 0L

  if(estimate.M) {
    opt <- optim(c(start$q, start$M), MLe, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, eff = tmb.dat$Effort,
                 LH = tmb.dat$LH, Lc = tmb.dat$Lc, eff_init = tmb.dat$eff_init,
                 n_age = tmb.dat$n_age, n_season = n_season, obs_season = obs_season,
                 timing = timing, logpar = tmb.dat$logpar,
                 method = "BFGS", control = list(maxit = 1e7))
    data.pred <- MLepred(opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, eff = tmb.dat$Effort,
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, eff_init = tmb.dat$eff_init,
                         n_age = tmb.dat$n_age, n_season = n_season, obs_season = obs_season,
                         timing = timing, logpar = tmb.dat$logpar)
    opt$Lpred <- data.pred[[1]]
    opt$par <- c(opt$par, data.pred[[2]])
    MLehess <- hessian(MLefullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, eff = tmb.dat$Effort,
                       LH = tmb.dat$LH, Lc = tmb.dat$Lc, eff_init = tmb.dat$eff_init,
                       n_age = tmb.dat$n_age, n_season = n_season, obs_season = obs_season,
                       timing = timing, logpar = tmb.dat$logpar)
    covariance <- solve(MLehess)
    opt$corr <- cov2cor(covariance)
    opt$gradient <- grad(MLefullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, eff = tmb.dat$Effort,
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, eff_init = tmb.dat$eff_init,
                         n_age = tmb.dat$n_age, n_season = n_season, obs_season = obs_season,
                         timing = timing, logpar = tmb.dat$logpar)

    if(!log.par) std.err <- sqrt(diag(covariance))
    if(log.par) {
      variance.normal <- deltamethod(list(~exp(x1), ~exp(x2)), opt$par[1:2], covariance[1:2, 1:2],
                                     ses = FALSE)
      std.err <- sqrt(c(diag(variance.normal)[1:2], diag(covariance)[3]))
    }
    results.matrix <- matrix(c(opt$par, std.err), ncol = 2)
    rownames(results.matrix) <- rownames(opt$corr) <- colnames(opt$corr) <- names(opt$par) <- names(opt$gradient) <- c("q", "M", "sigma")
  } else {
    opt <- optim(start$q, MLefixM, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, eff = tmb.dat$Effort,
                 LH = tmb.dat$LH, Lc = tmb.dat$Lc, eff_init = tmb.dat$eff_init,
                 n_age = tmb.dat$n_age, n_season = n_season, obs_season = obs_season,
                 timing = timing, logpar = tmb.dat$logpar,
                 method = "BFGS", control = list(maxit = 1e7))
    data.pred <- MLepred(c(opt$par, tmb.dat$LH[4]), Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, eff = tmb.dat$Effort,
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, eff_init = tmb.dat$eff_init,
                         n_age = tmb.dat$n_age, n_season = n_season, obs_season = obs_season,
                         timing = timing, logpar = tmb.dat$logpar)
    opt$Lpred <- data.pred[[1]]
    opt$par <- c(opt$par, data.pred[[2]])
    MLehess <- hessian(MLefixMfullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, eff = tmb.dat$Effort,
                       LH = tmb.dat$LH, Lc = tmb.dat$Lc, eff_init = tmb.dat$eff_init,
                       n_age = tmb.dat$n_age, n_season = n_season, obs_season = obs_season,
                       timing = timing, logpar = tmb.dat$logpar)
    covariance <- solve(MLehess)
    opt$corr <- cov2cor(covariance)
    opt$gradient <- grad(MLefixMfullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, eff = tmb.dat$Effort,
                         LH = tmb.dat$LH, Lc = tmb.dat$Lc, eff_init = tmb.dat$eff_init,
                         n_age = tmb.dat$n_age, n_season = n_season, obs_season = obs_season,
                         timing = timing, logpar = tmb.dat$logpar)

    if(!log.par) std.err <- sqrt(diag(covariance))
    if(log.par) {
      variance.normal <- deltamethod(~exp(x1), opt$par[1], covariance[1,1], ses = FALSE)
      std.err <- sqrt(c(variance.normal, diag(covariance)[2]))
    }
    results.matrix <- matrix(c(opt$par, std.err), ncol = 2)
    rownames(results.matrix) <- rownames(opt$corr) <- colnames(opt$corr) <- names(opt$par) <- names(opt$gradient) <- c("q", "sigma")
  }
  colnames(results.matrix) <- c("Estimate", "Std. Err.")
  time.series <- data.frame(Year = MLZ_data@Year, MeanLength = MLZ_data@MeanLength,
                            ss = MLZ_data@ss, Predicted = opt$Lpred)
  time.series$Residual <- time.series$MeanLength - time.series$Predicted
  time.series$F <- results.matrix[1,1] * tmb.dat$Effort
  if(estimate.M) time.series$M <- rep(results.matrix[2,1], nrow(time.series))
  if(!estimate.M) time.series$M <- rep(tmb.dat$LH[4], nrow(time.series))
  time.series$Z <- time.series$F + time.series$M

  MLZ_model <- new("MLZ_model", Stock = MLZ_data@Stock, Model = "MLeffort", time.series = time.series,
                   estimates = results.matrix, negLL = opt$value, n.species = 1L, opt = opt,
                   length.units = MLZ_data@length.units)
  if(figure) plot(MLZ_model)
  return(MLZ_model)
}
