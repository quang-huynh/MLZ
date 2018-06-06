#' Grid search for the mean length estimator
#'
#' A grid search is performed over the time series, which can be used to identify local and global minima. A
#' plot of the likelihood surface is also created similar to Figure 6 of Gedamke and Hoenig (2006) or
#' Figure 3 of Huynh et al. (2017).
#'
#' @param MLZ_data An object of class \code{MLZ_data}.
#' @param ncp The number of change points.
#' @param startZ A vector of length \code{ncp+1} as the starting value of total mortality rate used in the grid search.
#' @param min.time The minimum number of years between change points. Only used if \code{ncp > 1}.
#' @param parallel Whether grid search is performed using parallel processing.
#' @param figure If \code{TRUE}, creates a plot of the likelihood over the grid search. Only used
#' if \code{ncp = 1} or \code{2}.
#' @param color If \code{TRUE}, creates a color plot for the likelihood surface. Only used if
#' \code{ncp = 2}.
#' @return A matrix of change points with the negative log-likelihood values.
#' @examples
#' \dontrun{
#' data(Goosefish)
#' profile_ML(Goosefish, ncp = 1)
#' profile_ML(Goosefish, ncp = 2)
#' }
#' @references 
#' Gedamke, T. and Hoenig, J.M. 2006. Estimating mortality from mean length data in
#' nonequilibrium situations, with application to the assessment of goosefish.
#' Transactions of the American Fisheries Society 135:476-487.
#' 
#' Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. 2017. Multispecies Extensions
#' to a Nonequilibrium Length-Based Mortality Estimator. Marine and Coastal Fisheries 9:68-78.
#'
#' @export
profile_ML <- function(MLZ_data, ncp, startZ = rep(0.5, ncp+1), min.time = 3, 
                       parallel = ifelse(ncp > 2, TRUE, FALSE), figure = TRUE, color = TRUE) {
  
  parallel <- as.logical(parallel)
  ncp <- as.integer(ncp)
  if(ncp == 0) stop("profile_ML is not needed with zero change points.")
  
  year.range <- MLZ_data@Year - MLZ_data@Year[1] + 1
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
  
  years <- MLZ_data@Year
  max.years <- data.frame(Year = min(years):max(years))
  data.df <- summary_ML(MLZ_data)
  full <- left_join(max.years, data.df, by = "Year")
  
  tmb.dat <- list(model = "ML", Linf = MLZ_data@vbLinf, K = MLZ_data@vbK, Lc = MLZ_data@Lc, nbreaks = ncp,
                  Lbar = full$MeanLength, ss = full$ss)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss <= 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0 | tmb.dat$ss <= 0] <- 0
  
  if(length(MLZ_data@M) > 0) Z.limit <- MLZ_data@M else Z.limit <- 0.01
  
  start.yearZ <- split(ygrid, 1:nrow(ygrid))
  tmb.start <- Map(function(x, y) list(Z = x, yearZ = as.numeric(y)), y = start.yearZ, 
                   MoreArgs = list(x = startZ))
  obj <- Map(MakeADFun, parameters = tmb.start, 
             MoreArgs = list(data = tmb.dat, map = list(yearZ = factor(rep(NA, ncp))), 
                             DLL = "MLZ", silent = TRUE))
  if(!parallel) opt <- lapply(obj, function(x) try(nlminb(x$par, x$fn, x$gr, lower = rep(Z.limit, ncp + 1))))
  if(parallel) {
    cl <- makeCluster(detectCores())
    clusterExport(cl, c('Z.limit', 'ncp'), envir = environment())
    opt <- parLapply(cl, obj, function(x) try(nlminb(x$par, x$fn, x$gr, lower = rep(Z.limit, ncp + 1))))
    stopCluster(cl)
  }
  nll.vec <- vapply(opt, returnNAobjective, numeric(1))
  
  output <- as.data.frame(ygrid + MLZ_data@Year[1] - 1)
  output$negLL <- nll.vec
  names(output) <- c(paste0("Year", 1:ncp), "negLL")
  rownames(output) <- 1:nrow(output)
  
  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(list = old_par), add = TRUE)
    par(las = 1)

    if(ncp == 1) plot(negLL ~ Year1, output, xlab = "Change Point",
                      ylab = "Negative log likelihood", typ = "o", pch = 16, las = 1)
    if(ncp == 2) {
      z.matrix <- acast(output, Year1 ~ Year2, value.var = "negLL")
      min.index <- which.min(nll.vec)
      x.vec <- as.numeric(rownames(z.matrix))
      y.vec <- as.numeric(colnames(z.matrix))
      if(color) filled.contour(x = x.vec, y = y.vec, z = z.matrix,
                               xlab = "First Change Point", ylab = "Second Change Point",
                               plot.axes = {
                                 axis(1);
                                 axis(2);
                                 lines(x.vec, x.vec + min.time);
                                 points(output[min.index, 1], output[min.index, 2], pch = 16, cex = 1.5);
                                 segments(x0 = x.vec[1], x1 = output[min.index, 1], y0 = output[min.index, 2], lty = 2, lwd = 2);
                                 segments(x0 = output[min.index, 1], y0 = y.vec[1], y1 = output[min.index, 2], lty = 2, lwd = 2)
                               }
                               )
      if(!color) {
        contour(x = x.vec, y = y.vec, z = z.matrix,
                xlab = "First Change Point", ylab = "Second Change Point", labcex = 1, las = 1)
        lines(x.vec, x.vec + min.time)
        points(output[min.index, 1], output[min.index, 2], pch = 16, cex = 1.5)
        segments(x0 = x.vec[1], x1 = output[min.index, 1], y0 = output[min.index, 2], lty = 2, lwd = 2)
        segments(x0 = output[min.index, 1], y0 = y.vec[1], y1 = output[min.index, 2], lty = 2, lwd = 2)
      }
      title("Negative log-likelihood surface")
    }
  }
  return(output)

}

#' Grid search for the mean length with catch rate estimator
#'
#' A grid search is performed over the time series, which can be used to identify local and global minima. A
#' plot of the likelihood surface is also created similar to Figure 6 of Gedamke and Hoenig (2006) or
#' Figure 3 of Huynh et al. (2017).
#' 
#' @param MLZ_data An object of class \code{MLZ_data}.
#' @param ncp The number of change points.
#' @param CPUE.type Indicates whether CPUE time series is abundance or biomass based.
#' @param loglikeCPUE Indicates whether the log-likelihood for the CPUE will be lognormally or
#' normally distributed.
#' @param startZ A vector of length \code{ncp+1} as the starting value of total mortality rate used in the grid search.
#' @param parallel Whether the grid search is performed with parallel processing.
#' @param min.time The minimum number of years between change points. Only used if \code{ncp > 1}.
#' @param figure If \code{TRUE}, creates a plot of the likelihood over the grid search. Only used
#' if \code{ncp = 1} or \code{2}.
#' @param color If \code{TRUE}, creates a color plot for the likelihood surface. Only used if
#' \code{ncp = 2}.
#' @return A matrix of change points with the total negative log-likelihood values and values 
#' from the mean lengths and catch rates.
#' 
#' @references 
#' Gedamke, T. and Hoenig, J.M. 2006. Estimating mortality from mean length data in
#' nonequilibrium situations, with application to the assessment of goosefish.
#' Transactions of the American Fisheries Society 135:476-487.
#' 
#' Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. 2017. Multispecies Extensions
#' to a Nonequilibrium Length-Based Mortality Estimator. Marine and Coastal Fisheries 9:68-78.
#'
#' @examples
#' \dontrun{
#' data(MuttonSnapper)
#' profile_MLCR(MuttonSnapper, ncp = 1, CPUE.type = 'WPUE')
#' }
#' @export
profile_MLCR <- function(MLZ_data, ncp, CPUE.type = c(NA, "NPUE", "WPUE"), 
                         loglikeCPUE = c("normal", "lognormal"), startZ = rep(0.5, ncp+1), min.time = 3, 
                         parallel = ifelse(ncp > 2, TRUE, FALSE), figure = TRUE, color = TRUE) {
  
  CPUE.type <- match.arg(CPUE.type)
  if(is.na(CPUE.type)) stop("Argument 'CPUE.type' must be identified in either: weight 'WPUE' or abundance 'NPUE'.")
  loglikeCPUE <- match.arg(loglikeCPUE)
  parallel <- as.logical(parallel)
  ncp <- as.integer(ncp)
  if(ncp == 0) stop("profile_MLCR is not needed with zero change points.")

  year.range <- MLZ_data@Year - MLZ_data@Year[1] + 1
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

  
  years <- MLZ_data@Year
  max.years <- data.frame(Year = min(years):max(years))
  data.df <- data.frame(Year = MLZ_data@Year, MeanLength = MLZ_data@MeanLength, ss = MLZ_data@ss,
                        CPUE = MLZ_data@CPUE)
  full <- left_join(max.years, data.df, by = "Year")
  
  tmb.dat <- list(model = "MLCR", Linf = MLZ_data@vbLinf, K = MLZ_data@vbK, Lc = MLZ_data@Lc, nbreaks = ncp,
                  Lbar = full$MeanLength, ss = full$ss, CPUE = full$CPUE, 
                  CPUEisnormal = ifelse(loglikeCPUE == "normal", 1L, 0L))
  if(CPUE.type == "WPUE") {
    tmb.dat$b <- MLZ_data@lwb
    tmb.dat$isWPUE <- 1L
  } else {
    tmb.dat$b <- 1e-4
    tmb.dat$isWPUE <- 0L
  }
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss <= 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0 | tmb.dat$ss <= 0] <- 0
  tmb.dat$CPUE[is.na(tmb.dat$CPUE) | tmb.dat$CPUE < 0] <- -99
  
  if(length(MLZ_data@M) > 0) Z.limit <- MLZ_data@M else Z.limit <- 0.01
  
  start.yearZ <- split(ygrid, 1:nrow(ygrid))
  tmb.start <- Map(function(x, y) list(Z = x, yearZ = as.numeric(y)), y = start.yearZ, 
                   MoreArgs = list(x = startZ))
  obj <- Map(MakeADFun, parameters = tmb.start, 
             MoreArgs = list(data = tmb.dat, map = list(yearZ = factor(rep(NA, ncp))), 
                             DLL = "MLZ", silent = TRUE))
  if(!parallel) opt <- lapply(obj, function(x) try(nlminb(x$par, x$fn, x$gr, lower = rep(Z.limit, ncp + 1))))
  if(parallel) {
    cl <- makeCluster(detectCores())
    opt <- parLapply(cl, obj, function(x) try(nlminb(x$par, x$fn, x$gr, lower = rep(Z.limit, ncp + 1))))
    stopCluster(cl)
  }
  nll.vec <- vapply(opt, returnNAobjective, numeric(1))
  nll.data <- t(vapply(obj, function(x) x$report(x$env$last.par.best)$nllc, numeric(2)))
  
  output <- as.data.frame(ygrid + MLZ_data@Year[1] - 1)
  output <- cbind(output, nll.vec, nll.data)
  names(output) <- c(paste0("Year", 1:ncp), "negLL", "ML", "CR")
  rownames(output) <- 1:nrow(output)
  
  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(list = old_par), add = TRUE)
    par(las = 1)

    if(ncp == 1) {
      nll.vec2 <- nll.vec - min(nll.vec, na.rm = TRUE)
      nll.data.min <- apply(nll.data, 2, min, na.rm = TRUE)
      nll.data2 <- t(t(nll.data) - nll.data.min)
      range.y <- range(c(nll.vec2, nll.data2), finite = TRUE) * c(0.9, 1.1)

      plot(output$Year1, nll.vec2, xlab = "Change Point", ylab = "Change in negative log likelihood",
           typ = "l", pch = 16, las = 1, lwd = 3, ylim = range.y)
      lines(output$Year1, nll.data2[, 1], lty = 2, lwd = 2)
      lines(output$Year1, nll.data2[, 2], lty = 3, lwd = 2)
      legend("topright", c("Mean Length", "CPUE", "Total Likelihood"), lwd = c(2,2,3), lty = c(2,3,1))
    }

    if(ncp == 2) {
      if(color) plot_mlcr_2ncp_filled.contour(output, min.time)
      if(!color) {
        layout(matrix(c(1,2,3,3,3,3), nrow = 2))
        plot_mlcr_2ncp_contour(output, min.time, value.var = "ML")
        plot_mlcr_2ncp_contour(output, min.time, value.var = "CR")
        plot_mlcr_2ncp_contour(output, min.time, value.var = "negLL")
      }
    }
  }
  return(output)
}

#' Grid search for the multispecies mean length estimator
#'
#' A grid search is performed over the time series, which can be used to identify local and global minima. A
#' plot of the likelihood surface is also created similar to Figure 6 of Gedamke and Hoenig (2006) or
#' Figure 3 of Huynh et al. (2017).
#' 
#' @param MLZ.list A list containing an object of class \code{MLZ_data} for each species or stock.
#' @param ncp The number of change points.
#' @param model The name of the multispecies model for the grid search.
#' @param startZ1 A vector of length \code{ncp+1} as the starting value of total mortality rate used in the grid search.
#' @param parallel Whether the grid search is performed with parallel processing.
#' @param min.time The minimum number of years between change points. Only used if \code{ncp > 1}.
#' @param figure If \code{TRUE}, creates a plot of the likelihood over the grid search. Only used
#' if \code{ncp = 1} or \code{2}.
#' @param color If \code{TRUE}, creates a color plot for the likelihood surface. Only used if
#' \code{ncp = 2}.#' 
#' @return A matrix of change points with the total negative log-likelihood values and values 
#' from the each species.
#' 
#' @references 
#' Gedamke, T. and Hoenig, J.M. 2006. Estimating mortality from mean length data in
#' nonequilibrium situations, with application to the assessment of goosefish.
#' Transactions of the American Fisheries Society 135:476-487.
#' 
#' Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. 2017. Multispecies Extensions
#' to a Nonequilibrium Length-Based Mortality Estimator. Marine and Coastal Fisheries 9:68-78.
#'
#' @examples
#' \dontrun{
#' data(PRSnapper)
#' profile_MLmulti(PRSnapper, ncp = 1, model = "MSM1")
#' }
#' @export
profile_MLmulti <- function(MLZ.list, ncp, model = c("SSM", "MSM1", "MSM2", "MSM3"), 
                            startZ1 = rep(0.5, length(MLZ.list)), 
                            parallel = ifelse(ncp > 2, TRUE, FALSE), min.time = 3, figure = TRUE, color = TRUE) {

  model <- match.arg(model)
  parallel <- as.logical(parallel)
  ncp <- as.integer(ncp)
  if(ncp == 0) stop("profile_MLmulti is not needed with zero changepoints.")
  
  nspec <- length(MLZ.list)
  years <- lapply(MLZ.list, getElement, "Year")
  years <- unique(do.call(c, years))

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

  max.years <- data.frame(Year = min(years):max(years))
  Lbar.df <- lapply(MLZ.list, summary)
  Lbar.df <- do.call(rbind, Lbar.df)
  full <- left_join(max.years, Lbar.df, by = "Year")
  full$Sp <- rep(1:nspec, nrow(max.years))
  Lbar <- acast(full, Year ~ Sp, value.var = "MeanLength")
  ss <- acast(full, Year ~ Sp, value.var = "ss")

  Linf <- vapply(MLZ.list, getElement, numeric(1), "vbLinf")
  K <- vapply(MLZ.list, getElement, numeric(1), "vbK")
  Lc <- vapply(MLZ.list, getElement, numeric(1), "Lc")

  tmb.dat <- list(Linf = Linf, K = K, Lc = Lc, nbreaks = ncp, nspec = nspec, Lbar = Lbar, ss = ss)
  if("MSM2" %in% model || "MSM3" %in% model) tmb.dat$M <- vapply(MLZ.list, getElement, numeric(1), "M")
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss == 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0] <- 0
  
  start.yearZ <- split(ygrid, 1:nrow(ygrid))
  if("SSM" %in% model || "MSM1" %in% model) {
    tmb.dat$model <- "MSM1S" 
    lower.limit <- rep(vapply(MLZ.list, get_M_MSM1S, numeric(1)), each = ncp + 1)
    
    startZ <- matrix(startZ1, ncol = nspec, nrow = ncp+1, byrow = TRUE)
    start.yearZ <- lapply(start.yearZ, function(x) matrix(as.numeric(x), ncol = nspec, nrow = ncp))
    tmb.start <- Map(function(x, y) list(Z = x, yearZ = y), y = start.yearZ, MoreArgs = list(x = startZ))
    obj <- Map(MakeADFun, parameters = tmb.start, DLL = "MLZ", silent = TRUE,
               MoreArgs = list(data = tmb.dat, map = list(yearZ = factor(rep(NA, nspec * ncp)))))
    lower.limit <- rep(vapply(MLZ.list, get_M_MSM1S, numeric(1)), each = ncp)
  }
  if("MSM2" %in% model) {
    tmb.dat$model <- "MSM23"
    Z1.limit <- tmb.dat$M
    delta.limit <- rep(0, ncp)
    eps.limit <- rep(0, nspec - 1)
    lower.limit <- c(Z1.limit, delta.limit, eps.limit)
    
    tmb.start <- Map(function(x, y, z, zz) list(Z1 = x, yearZ = as.numeric(y), delta = z, epsilon = zz), 
                     y = start.yearZ, MoreArgs = list(x = startZ1, z = rep(1, ncp), zz = rep(1, nspec-1)))
    obj <- Map(MakeADFun, parameters = tmb.start, DLL = "MLZ", silent = TRUE,
               MoreArgs = list(data = tmb.dat, map = list(yearZ = factor(rep(NA, ncp)))))
  }
  if("MSM3" %in% model) {
    tmb.dat$model <- "MSM23"
    Z1.limit <- tmb.dat$M
    delta.limit <- rep(0, ncp)
    lower.limit <- c(Z1.limit, delta.limit)
    
    tmb.start <- Map(function(x, y, z, zz) list(Z1 = x, yearZ = as.numeric(y), delta = z, epsilon = zz), 
                     y = start.yearZ, MoreArgs = list(x = startZ1, z = rep(1, ncp), zz = rep(1, nspec-1)))
    obj <- Map(MakeADFun, parameters = tmb.start, 
               MoreArgs = list(data = tmb.dat, DLL = "MLZ", silent = TRUE,
                               map = list(yearZ = factor(rep(NA, ncp)), epsilon = factor(rep(NA, nspec-1)))))
  }
  if(!parallel) opt <- lapply(obj, function(x) try(nlminb(x$par, x$fn, x$gr, lower = lower.limit)))
  if(parallel) {
    cl <- makeCluster(detectCores())
    clusterExport(cl, list("lower.limit"))
    opt <- parLapply(cl, obj, function(x) try(nlminb(x$par, x$fn, x$gr, lower = lower.limit)))
    stopCluster(cl)
  }
  nll.vec <- vapply(opt, returnNAobjective, numeric(1))
  nll.data <- t(vapply(obj, function(x) x$report(x$env$last.par.best)$nllc, numeric(nspec)))
  
  output <- as.data.frame(ygrid + min(years) - 1)
  output <- cbind(output, nll.vec, nll.data)
  names(output) <- c(paste0("Year", 1:ncp), "negLL", paste0("Species_", 1:nspec))
  rownames(output) <- 1:nrow(output)
  
  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(list = old_par), add = TRUE)
    par(las = 1)

    sp.name <- lapply(MLZ.list, getElement, "Stock")
    sp.ind <- vapply(sp.name, length, numeric(1))
    num <- which(sp.ind == 0L)
    sp.name[num] <- paste0("Species_", num)

    if(ncp == 1) {
      color.vec <- rich.colors(nspec)
      loglikesp <- nll.data
      nll.vec2 <- output$negLL - min(output$negLL)
      min.loglikesp <- apply(loglikesp, 2, min)
      loglikesp2 <- t(t(loglikesp) - min.loglikesp)

      plot(output$Year1, nll.vec2, xlab = "Change Point", ylim = c(0, ceiling(max(cbind(nll.vec2, loglikesp2), na.rm = TRUE))),
           ylab = "Change in negative log likelihood", typ = "o", pch = 16, las = 1, lwd = 2)
      for(i in 1:nspec) lines(output$Year1, loglikesp2[, i], typ = "o", pch = 16, col = color.vec[i])
      legend("topright", c("Total Likelihood", do.call(c, sp.name)), col = c("black", color.vec),
             pch = 16, lty = 1)
    }
    if(ncp == 2) {
      for(i in 1:(nspec+1)) {
        if(i < nspec+1) {
          z.matrix <- acast(output[, c(1,2,3+i)], Year1 ~ Year2, value.var = paste0("Species_", i))
          min.index <- which.min(output[, 3+i])
        }
        if(i == nspec+1) {
          z.matrix <- acast(output[, 1:3], Year1 ~ Year2, value.var = "negLL")
          min.index <- which.min(output[, 3])
        }
        x.vec <- as.numeric(rownames(z.matrix))
        y.vec <- as.numeric(colnames(z.matrix))
        xyline <- union(x.vec, y.vec)
        if(color) filled.contour(x = x.vec, y = y.vec, z = z.matrix,
                                 xlab = "First Change Point", ylab = "Second Change Point",
                                 plot.axes = {
                                   axis(1);
                                   axis(2);
                                   lines(x.vec, x.vec + min.time);
                                   points(output[min.index, 1], output[min.index, 2], pch = 16, cex = 1.5);
                                   segments(x0 = x.vec[1], x1 = output[min.index, 1], y0 = output[min.index, 2], lty = 2, lwd = 2);
                                   segments(x0 = output[min.index, 1], y0 = y.vec[1], y1 = output[min.index, 2], lty = 2, lwd = 2);
                                 }
        )
        if(!color) {
          contour(x = x.vec, y = y.vec, z = z.matrix,
                  xlab = "First Change Point", ylab = "Second Change Point", labcex = 1, las = 1)
          lines(x.vec, x.vec + min.time)
          points(output[min.index, 1], output[min.index, 2], pch = 16, cex = 1.5)
          segments(x0 = x.vec[1], x1 = output[min.index, 1], y0 = output[min.index, 2], lty = 2, lwd = 2)
          segments(x0 = output[min.index, 1], y0 = y.vec[1], y1 = output[min.index, 2], lty = 2, lwd = 2)
        }
        if(i < nspec+1) title(paste0("Negative log-likelihood surface for ", sp.name[i]))
        if(i == nspec+1) title("Total negative log-likelihood surface")
      }
    }
  }

  return(output)
}


