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
  spawn <- match.arg(spawn)
  if(spawn == "continuous") spawn.cont <- 1L else spawn.cont <- 0L
  
  years <- MLZ_data@Year
  max.years <- data.frame(Year = min(years):max(years))
  data.df <- summary_ML(MLZ_data)
  full <- left_join(max.years, data.df, by = "Year")

  tmb.dat <- list(Linf = MLZ_data@vbLinf, K = MLZ_data@vbK, Lc = MLZ_data@Lc, nbreaks = ncp,
                  Lbar = full$MeanLength, ss = full$ss, spCont = spawn.cont)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss <= 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0 | tmb.dat$ss <= 0] <- 0
  
  if(length(MLZ_data@M) > 0) Z.limit <- MLZ_data@M else Z.limit <- 0.01

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
    obj <- MakeADFun(data = tmb.dat, parameters = start, hessian = TRUE, DLL = "MLeq", silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, lower = Z.limit)
    sdrep <- sdreport(obj)
    
    results.matrix <- summary(sdrep)
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
        grid.output <- profile_ML(MLZ_data, ncp, min.time = min.time, parallel = as.logical(parallel), 
                                  spawn = spawn, figure = FALSE)
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
    obj <- MakeADFun(data = tmb.dat, parameters = start, hessian = TRUE, DLL = "ML", silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, lower = c(rep(Z.limit, ncp + 1), rep(-Inf, ncp)))
    sdrep <- sdreport(obj)

    results.matrix <- summary(sdrep)
    Z.name <- paste0("Z[", 1:(ncp+1), "]")
    yearZ.name <- paste0("yearZ[", 1:ncp, "]")
    rownames(results.matrix) <- c(Z.name, yearZ.name, "sigma")
    year.ind <- grep("yearZ", rownames(results.matrix))
    results.matrix[year.ind, 1] <- results.matrix[year.ind, 1] + MLZ_data@Year[1] - 1
  }
  if(any(results.matrix[1:(ncp+1), 1] == Z.limit)) {
    warning("There are mortality estimates at boundary (Z = 0.01 or M).")
  }
  time.series <- full
  time.series$Predicted <- obj$report()$Lpred
  time.series$Residual <- time.series$MeanLength - time.series$Predicted
  opt$message <- c(opt$message, paste0("Spawning is assumed to be ", spawn, " in model."))
  class(sdrep) <- "list"
  MLZ_model <- new("MLZ_model", Stock = MLZ_data@Stock, Model = "ML", time.series = time.series,
                   estimates = results.matrix, negLL = opt$objective, n.changepoint = ncp, n.species = 1L,
                   obj = obj, opt = opt, sdrep = sdrep, length.units = MLZ_data@length.units)
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
  
  CPUE.type <- match.arg(CPUE.type)
  if(is.na(CPUE.type)) stop("Argument CPUE.type must be either 'WPUE' (weight-based) or 'NPUE' (abundance/numbers-based)")
  loglikeCPUE <- match.arg(loglikeCPUE)
  ll.int <- ifelse(loglikeCPUE == "lognormal", 0L, 1L)
  
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
  tmb.dat <- list(Linf = MLZ_data@vbLinf, K = MLZ_data@vbK, Lc = MLZ_data@Lc, nbreaks = ncp,
                  Lbar = full$MeanLength, ss = full$ss, CPUE = full$CPUE, 
                  loglikeCPUE = ll.int, spCont = spawn.cont)
  if(CPUE.type == "WPUE") tmb.dat$b <- MLZ_data@lwb
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss <= 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0 | tmb.dat$ss <= 0] <- 0
  tmb.dat$CPUE[is.na(tmb.dat$CPUE) | tmb.dat$CPUE < 0] <- -99
  
  if(length(MLZ_data@M) > 0) Z.limit <- MLZ_data@M else Z.limit <- 0.01
  
  if(ncp == 0) {
    if(CPUE.type == "NPUE") {
      tmb.dat$b <- 3
      tmb.dat$isNPUE <- 1L
    } else tmb.dat$isNPUE <- 0L
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
    obj <- MakeADFun(data = tmb.dat, parameters = start, hessian = TRUE, DLL = "MLCReq", silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, lower = Z.limit)
    sdrep <- sdreport(obj)
    
    results.matrix <- summary(sdrep)
    rownames(results.matrix) <- c("Z", "q", "sigmaL", "sigmaI")
    
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
        grid.output <- profile_MLCR(MLZ_data, ncp, CPUE.type, loglikeCPUE, spawn = spawn, 
                                    min.time = min.time, parallel = as.logical(parallel), figure = FALSE)
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
    
    fx <- paste0("ML", CPUE.type)
    obj <- MakeADFun(data = tmb.dat, parameters = start, hessian = TRUE, DLL = fx, silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, lower = c(rep(Z.limit, ncp + 1), rep(-Inf, ncp)))
    sdrep <- sdreport(obj)
    
    results.matrix <- summary(sdrep)
    Z.name <- paste0("Z[", 1:(ncp+1), "]")
    yearZ.name <- paste0("yearZ[", 1:ncp, "]")
    rownames(results.matrix) <- c(Z.name, yearZ.name, "q", "sigmaL", "sigmaI")
    year.ind <- grep("yearZ", rownames(results.matrix))
    results.matrix[year.ind, 1] <- results.matrix[year.ind, 1] + MLZ_data@Year[1] - 1
    
  }
  
  if(any(results.matrix[1:(ncp+1), 1] == Z.limit)) {
    warning("There are mortality estimates at boundary (Z = 0.01 or M).")
  }
  time.series <- data.frame(Predicted.ML = obj$report()$Lpred, Predicted.CPUE = obj$report()$Ipred)
  time.series <- cbind(full, time.series)
  time.series$Residual.ML <- time.series$MeanLength - time.series$Predicted.ML
  time.series$Residual.CPUE <- time.series$CPUE - time.series$Predicted.CPUE

  opt$message <- c(opt$message, paste0("Catch rate assumed to be ", CPUE.type, "."), 
                   paste0("Catch rate likelihood used ", loglikeCPUE, " distribution."),
                   paste0("Spawning is assumed to be ", spawn, " in model."))
  class(sdrep) <- "list"
  MLZ_model <- new("MLZ_model", Stock = MLZ_data@Stock, Model = "MLCR", time.series = time.series,
                   estimates = results.matrix, negLL = opt$objective, n.changepoint = ncp, n.species = 1L,
                   obj = obj, opt = opt, sdrep = sdrep, length.units = MLZ_data@length.units)
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
#' If \code{ncp = 0} change points is specified, then the method simplifies to the Single Species Model.
#' The \code{start} list should contain a single entry:
#' \tabular{ll}{
#' \code{Z} \tab a vector with \code{length = N}.\cr
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
#' MLmulti(PRSnapper, ncp = 0, start = list(Z = rep(0.5, 3)))
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
  if(!is.list(MLZ.list)) stop("No list found.")
  class.test <- vapply(MLZ.list, inherits, TRUE, "MLZ_data")
  if(!all(class.test)) stop("Not all entries in list are of class 'MLZ_data'")
  length.units <- vapply(MLZ.list, getElement, c("x"), "length.units")
  length.units <- unique(length.units)
  if(length(length.units) > 1) stop("Different length units among species. All lengths should be the same units, e.g. cm.")
  
  ncp <- as.integer(ncp)
  spawn <- match.arg(spawn)
  if(spawn == "continuous") spawn.cont <- 1L else spawn.cont <- 0L
  
  nspec <- length(MLZ.list)
  years <- lapply(MLZ.list, getElement, "Year")
  years <- unique(do.call(c, years))

  max.years <- data.frame(Year = min(years):max(years))
  Lbar.df <- lapply(MLZ.list, summary_ML)
  Lbar.df <- do.call(rbind, Lbar.df)
  full <- left_join(max.years, Lbar.df, by = "Year")
  full$Sp <- rep(1:nspec, nrow(max.years))
  Lbar <- acast(full, Year ~ Sp, value.var = "MeanLength")
  ss <- acast(full, Year ~ Sp, value.var = "ss")
  Lbar[is.na(ss) | ss == 0] <- -99
  ss[is.na(ss) | Lbar < 0] <- 0

  Linf <- vapply(MLZ.list, getElement, numeric(1), "vbLinf")
  K <- vapply(MLZ.list, getElement, numeric(1), "vbK")
  Lc <- vapply(MLZ.list, getElement, numeric(1), "Lc")
  
  if(ncp == 0) {
    grid.search <- FALSE
    model <- "SSM"
    if(!is.null(start)) {
      if(!"Z" %in% names(start)) stop("Error in start list. Entry of name 'Z' not found.")
      if(length(start$Z) != nspec)
        stop("Length of 'Z' in start list is not equal to the number of species.")
    }
    else {
      Z <- K * (Linf - Lbar[1, ]) / (Lbar[1, ] - Lc)
      Z[is.na(Z) | Z <= 0 | is.infinite(Z)] <- 0.5
      Z[Z > 1] <- 1
      start <- list(Z = Z)
    }
    
    start <- Map(list, Z = start$Z)
    tmb.dat <- Map(list, Linf = Linf, K = K, Lc = Lc, Lbar = split(Lbar, rep(1:nspec, each = nrow(Lbar))),
                   ss = split(ss, rep(1:nspec, each = nrow(ss))), spCont = spawn.cont)
    obj <- Map(MakeADFun, data = tmb.dat, parameters = start, 
               MoreArgs = list(hessian = TRUE, DLL = "MLeq", silent = TRUE))
    opt <- lapply(obj, function(x) nlminb(x$par, x$fn, x$gr, x$he))
    sdrep <- lapply(obj, sdreport)
    
    res <- vapply(sdrep, function(x) summary(x)[1, ], numeric(2))
    res2 <- vapply(sdrep, function(x) summary(x)[2, ], numeric(2))
    results.matrix <- rbind(t(res), t(res2))
    rownames(results.matrix) <- c(paste0("Z[", 1:nspec, "]"), paste0("sigma[", 1:nspec, "]"))
    
    Lbar.df$Predicted <- as.numeric(vapply(obj, function(x) x$report()$Lpred, numeric(nrow(Lbar))))
  }
  if(ncp > 0) {
    tmb.dat <- list(Linf = Linf, K = K, Lc = Lc, nbreaks = ncp, nspec = nspec, Lbar = Lbar, ss = ss, 
                    spCont = spawn.cont)
    if("MSM2" %in% model || "MSM3" %in% model) tmb.dat$M <- vapply(MLZ.list, getElement, numeric(1), "M")

    
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
        new.yearZ <- matrix(start$yearZ, nrow = ncp, ncol = nspec)
        start$yearZ <- new.yearZ
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
        start$epsilon <- rep(1, nspec - 1)
      }
      start$yearZ <- start$yearZ - max.years[1, 1] + 1
    } else {
      if(grid.search) {
        pr <- profile_MLmulti(MLZ.list, ncp, model, spawn, parallel = as.logical(parallel), 
                              min.time = min.time, figure = FALSE)
        if(model == "SSM") {
          index.min <- apply(pr[, -c(1:(ncp+1))], 2, which.min)
          styearZ <- t(as.matrix(pr[index.min, 1:ncp]))
        } else {
          index.min <- which.min(pr$negLL)
          styearZ <- as.numeric(pr[index.min, 1:ncp])
        }
        styearZ <- styearZ - max.years[1, 1] + 1
      }
      else {
        styearZ <- nrow(max.years) * (1:ncp) / (ncp+1)
      }
      stZ1 <- K * (Linf - tmb.dat$Lbar[1, ]) / (tmb.dat$Lbar[1, ] - Lc)
      stZ1[is.na(stZ1) | stZ1 <= 0] <- 0.5
      stZ1[stZ1 > 1] <- 1
      if("SSM" %in% model || "MSM1" %in% model) {
        styearZ <- matrix(styearZ, nrow = ncp, ncol = nspec)
        stZ <- matrix(NA, nrow = ncp, ncol = nspec)
        for(i in 1:ncp) {
          tmp <- K * (Linf - tmb.dat$Lbar[styearZ[i, ], ]) / (tmb.dat$Lbar[styearZ[i, ], ] - Lc)
          stZ[i, ] <- diag(tmp)
        }
        stZ[is.na(stZ) | stZ <= 0] <- 0.5
        stZ[stZ > 1] <- 1
      }
      if("SSM" %in% model || "MSM1" %in% model) {
        stZ <- rbind(stZ1, stZ)
        start <- list(Z = stZ, yearZ = matrix(styearZ, nrow = ncp, ncol = nspec))
      }
      if("MSM2" %in% model || "MSM3" %in% model) {
        start <- list(Z1 = stZ1, delta = rep(1, ncp), epsilon = rep(1, nspec - 1), yearZ = styearZ)
      }
    }
    
    if("SSM" %in% model) {
      fx <- "MSM1S"    
      map <- list()
    }
    if("MSM1" %in% model) {
      fx <- "MSM1S"
      map <- list(yearZ = factor(rep(1:ncp, nspec)))
    }
    if("MSM2" %in% model) {
      fx <- "MSM23"
      map <- list()
    }
    if("MSM3" %in% model) {
      fx <- "MSM23"
      map <- list(epsilon = factor(rep(NA, nspec - 1)))
    }
    obj <- MakeADFun(data = tmb.dat, parameters = start, hessian = TRUE, 
                     map = map, DLL = fx, silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
    sdrep <- sdreport(obj)
    
    results.matrix <- summary(sdrep)
    year.ind <- grep("yearZ", rownames(results.matrix))
    results.matrix[year.ind, 1] <- results.matrix[year.ind, 1] + min(years) - 1
    
    if("SSM" %in% model || "MSM1" %in% model) {
      Z.name <- paste0(rep(paste0("Z[",1:(ncp+1)), nspec), rep(paste0(",", 1:nspec, "]"), each = ncp+1))
      if(model == "SSM") yearZ.name <- paste0(rep(paste0("yearZ[", 1:ncp), nspec), rep(paste0(",", 1:nspec, "]"), each = ncp))
      if(model == "MSM1") yearZ.name <- paste0("yearZ[", 1:ncp, "]")
      sigma.name <- paste0("sigma[", 1:nspec, "]")
      rownames(results.matrix) <- c(Z.name, yearZ.name, sigma.name)
    }
    
    if("MSM2" %in% model || "MSM3" %in% model) {
      results.matrix <- results.matrix[-c(1:nspec), ]
      delta.name <- paste0("delta[", 1:ncp, ",1]")
      yearZ.name <- paste0("yearZ[", 1:ncp, "]")
      Z.name <- paste0(rep(paste0("Z[",1:(ncp+1)), nspec), rep(paste0(",", 1:nspec, "]"), each = ncp+1))
      sigma.name <- paste0("sigma[", 1:nspec, "]")
      if(model == "MSM2") {
        epsilon.name <- paste0("epsilon[", 2:nspec, "]")
        rownames(results.matrix) <- c(delta.name, epsilon.name, yearZ.name, Z.name, sigma.name)
      }
      if(model == "MSM3") rownames(results.matrix) <- c(delta.name, yearZ.name, Z.name, sigma.name)
    }
    
    Lbar.df$Predicted <- as.numeric(obj$report()$Lpred)
    class(sdrep) <- "list"
  }
  
  if(any(results.matrix[grep("Z", rownames(results.matrix)), 1] < 0)) warning("There are negative mortality estimates from model.")
  
  Lbar.df$Residual <- Lbar.df$MeanLength - Lbar.df$Predicted
  sp.name <- lapply(MLZ.list, getElement, "Stock")
  sp.ind <- vapply(sp.name, length, numeric(1))
  num <- which(sp.ind == 0L)
  sp.name[num] <- paste0("Species_", num)
  Lbar.df$Stock <- rep(do.call(c, sp.name), each = length(years))
  
  opt$message <- c(opt$message, paste0("Spawning is assumed to be ", spawn, " in model."))
  
  if(ncp == 0) negLL <- sum(vapply(opt[1:nspec], getElement, numeric(1), "objective"))
  if(ncp > 0) negLL <- opt$objective
  MLZ_model <- new("MLZ_model", Stock = do.call(c, sp.name), Model = "MLmulti",
                   time.series = Lbar.df, estimates = results.matrix,
                   negLL = negLL, n.changepoint = ncp, n.species = nspec,
                   obj = obj, opt = opt, sdrep = sdrep, length.units = length.units)
  attr(MLZ_model, "multimodel") <- model
  attr(MLZ_model, "spawn") <- spawn
  if(grid.search) MLZ_model@grid.search <- pr
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
  
  years <- MLZ_data@Year
  max.years <- data.frame(Year = min(years):max(years))
  data.df <- data.frame(Year = MLZ_data@Year, MeanLength = MLZ_data@MeanLength, ss = MLZ_data@ss,
                        Effort = MLZ_data@Effort)
  full <- left_join(max.years, data.df, by = "Year")
  if(any(is.na(full$Effort))) {
    yr.missing <- full$Year[is.na(full$Effort)]
    stop(paste("Missing effort data in Year(s): ", yr.missing))
  }

  tmb.dat <- list(Linf = MLZ_data@vbLinf, K = MLZ_data@vbK, a0 = MLZ_data@vbt0, Lc = MLZ_data@Lc, 
                  Lbar = full$MeanLength, ss = full$ss, eff = full$Effort, 
                  eff_init = eff_init, n_age = n_age, n_season = n_season, obs_season = obs_season,
                  timing = timing)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss == 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0] <- 0

  if(estimate.M) {
    map <- list()
    start <- list(logq = start$q, logM = start$M)
  } else {    
    map <- list(logM = factor(NA))
    start <- list(logq = start$q, logM = MLZ_data@M)
  }
  
  if(log.par) {
    tmb.dat$logpar <- 1L
    start <- lapply(start, log)
  } else tmb.dat$logpar <- 0L
  
  obj <- MakeADFun(data = tmb.dat, parameters = start, hessian = TRUE, map = map, 
                   DLL = "MLe", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he)
  sdrep <- sdreport(obj)
  
  results.matrix <- summary(sdrep)
  if(!log.par) {
    if(estimate.M) ind.remove <- 1:2
    if(!estimate.M) ind.remove <- 1
    results.matrix <- results.matrix[-ind.remove, ]
    if(any(results.matrix[, 1] < 0)) warning("There are negative estimates from model.")
  }
  
  full$Predicted <- obj$report()$Lpred
  full$Residual <- full$MeanLength - full$Predicted
  full$Z <- obj$report()$Z
  if(estimate.M) {
    ind.M <- which(rownames(results.matrix) == "M")
    full$M <- rep(results.matrix[ind.M, 1], nrow(full))
  }
  if(!estimate.M) full$M <- rep(MLZ_data@M, nrow(full))
  full$F <- full$Z - full$M

  class(sdrep) <- "list"
  MLZ_model <- new("MLZ_model", Stock = MLZ_data@Stock, Model = "MLeffort", time.series = full,
                   estimates = results.matrix, negLL = opt$objective, n.species = 1L, 
                   obj = obj, opt = opt, sdrep = sdrep, length.units = MLZ_data@length.units)
  if(figure) plot(MLZ_model)
  return(MLZ_model)
}
