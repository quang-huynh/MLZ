#' Grid search for the mean length estimator
#'
#' A grid search is performed over the time series, which can be used to identify local and global minima. A
#' plot of the likelihood surface is also created similar to Figure 6 of Gedamke and Hoenig (2006).
#'
#' @param MLZ_data An object of class \code{MLZ_data}.
#' @param ncp The number of change points.
#' @param stZ The starting value of total mortality rate used in the grid search.
#' @param min.time The minimum number of years between change points. Only used if \code{ncp > 1}.
#' @param parallel Whether grid search is performed using parallel processing.
#' @param figure If \code{TRUE}, creates a plot of the likelihood over the grid search. Only used
#' if \code{ncp = 1} or \code{2}.
#' @param color If \code{TRUE}, creates a color plot for the likelihood surface. Only used if
#' \code{ncp = 2}.
#' @examples
#' data(Goosefish)
#' profile_ML(Goosefish, ncp = 1)
#' profile_ML(Goosefish, ncp = 2)
#'
#' @export
profile_ML <- function(MLZ_data, ncp, stZ = 0.5, min.time = 3, parallel = FALSE, figure = TRUE, color = TRUE) {
  ncp <- as.integer(ncp)
  startZ <- rep(stZ, ncp + 1)
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

  tmb.dat <- list(LH = c(MLZ_data@vbLinf, MLZ_data@vbK), Lc = MLZ_data@Lc, nbreaks = ncp,
                  Lbar = MLZ_data@MeanLength, ss = MLZ_data@time.series$ss)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss <= 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0 | tmb.dat$ss <= 0] <- 0

  if(!parallel) {
    nll.vec <- numeric(nrow(ygrid))
    for(i in 1:nrow(ygrid)) {
      opt <- try(optim(startZ, MLprofile, yearZ = as.numeric(ygrid[i, ]), Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                       LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                       method = "BFGS", control = list(maxit = 1e7)))
      if(inherits(opt, "try-error")) nll.vec[i] <- NA else nll.vec[i] <- opt$value
    }
  }
  if(parallel) {
    start.yearZ <- split(ygrid, 1:nrow(ygrid))
    cl <- makeCluster(detectCores())
    clusterExport(cl, ls(), envir = environment())
    start.yearZ <- parLapply(cl, start.yearZ, as.numeric)
    opt <- parLapply(cl, start.yearZ, function(x) try(optim(startZ, MLprofile, yearZ = x,
                                                        Lbar = tmb.dat$Lbar, ss = tmb.dat$ss, LH = tmb.dat$LH,
                                                        Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                                                        method = "BFGS", control = list(maxit = 1e7))))
    stopCluster(cl)
    nll.vec <- vapply(opt, returnNAoptvalue, numeric(1))
  }

  output <- as.data.frame(ygrid + MLZ_data@Year[1] - 1)
  output$negLL <- nll.vec
  names(output) <- c(paste0("Year", 1:ncp), "negLL")
  rownames(output) <- 1:nrow(output)

  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(list = old_par), add = TRUE)

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
                                 lines(x.vec, y.vec);
                                 points(output[min.index, 1], output[min.index, 2], pch = 16, cex = 1.5);
                                 segments(x0 = x.vec[1], x1 = output[min.index, 1], y0 = output[min.index, 2], lty = 2, lwd = 2);
                                 segments(x0 = output[min.index, 1], y0 = y.vec[1], y1 = output[min.index, 2], lty = 2, lwd = 2)
                               }
                               )
      if(!color) {
        contour(x = x.vec, y = y.vec, z = z.matrix,
                xlab = "First Change Point", ylab = "Second Change Point", labcex = 1, las = 1)
        lines(x.vec, y.vec)
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
#' A grid search is performed over the time series. Can be used to identify local and global minima.
#'
#' @param MLZ_data An object of class \code{MLZ_data}.
#' @param ncp The number of change points.
#' @param CPUE.type Indicates whether CPUE time series is abundance or biomass based.
#' @param loglikeCPUE Indicates whether the log-likelihood for the CPUE will be lognormally or
#' normally distributed.
#' @param stZ The starting value of total mortality rate used in the grid search.
#' @param parallel Whether the grid search is performed with parallel processing.
#' @param min.time The minimum number of years between change points. Only used if \code{ncp > 1}.
#' @param figure If \code{TRUE}, creates a plot of the likelihood over the grid search. Only used
#' if \code{ncp = 1} or \code{2}.
#' @param color If \code{TRUE}, creates a color plot for the likelihood surface. Only used if
#' \code{ncp = 2}.
#' @examples
#' data(MuttonSnapper)
#' profile_MLCR(MuttonSnapper, ncp = 1)
#' @export
profile_MLCR <- function(MLZ_data, ncp, CPUE.type = c("NPUE", "WPUE"), loglikeCPUE = c("lognormal", "normal"),
                         stZ = 0.5, min.time = 3, parallel = FALSE, figure = TRUE, color = TRUE) {

  CPUE.type <- match.arg(CPUE.type)
  loglikeCPUE <- match.arg(loglikeCPUE)
  ncp <- as.integer(ncp)
  startZ <- rep(stZ, ncp + 1)
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

  ll.int <- ifelse(loglikeCPUE == "lognormal", 0L, 1L)
  tmb.dat <- list(LH = c(MLZ_data@vbLinf, MLZ_data@vbK), Lc = MLZ_data@Lc, nbreaks = ncp,
                  Lbar = MLZ_data@MeanLength, ss = MLZ_data@ss,
                  CPUE = MLZ_data@CPUE, loglikeCPUE = ll.int)
  if(CPUE.type == "WPUE") tmb.dat$LH <- c(tmb.dat$LH, MLZ_data@lwb)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss <= 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0 | tmb.dat$ss <= 0] <- 0
  tmb.dat$CPUE[is.na(tmb.dat$CPUE) | tmb.dat$CPUE < 0] <- -99

  fx <- get(paste0("ML", CPUE.type, "profile"))
  gx <- get(paste0("ML", CPUE.type, "pred"))
  if(!parallel) {
    nll.vec <- numeric(nrow(ygrid))
    nll.data <- matrix(0, nrow = nrow(ygrid), ncol = 2)
    for(i in 1:nrow(ygrid)) {
      opt <- try(optim(startZ, fx, yearZ = as.numeric(ygrid[i, ]), Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                   CPUE = tmb.dat$CPUE, LH = tmb.dat$LH, Lc = tmb.dat$Lc,
                   nbreaks = tmb.dat$nbreaks, loglikeCPUE = tmb.dat$loglikeCPUE,
                   method = "BFGS", control = list(maxit = 1e7)))
      if(inherits(opt, "try-error")) {
        nll.vec[i] <- NA
        nll.data[i, ] <- rep(NA, 2)
      } else {
        data.pred <- gx(c(opt$par, as.numeric(ygrid[i, ])), Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                        CPUE = tmb.dat$CPUE, LH = tmb.dat$LH, Lc = tmb.dat$Lc,
                        nbreaks = tmb.dat$nbreaks, loglikeCPUE = tmb.dat$loglikeCPUE)
        nll.vec[i] <- opt$value
        nll.data[i, ] <- -1 * data.pred[[5]]
      }
    }
  }
  if(parallel) {
    start.yearZ <- split(ygrid, 1:nrow(ygrid))
    cl <- makeCluster(detectCores())
    clusterExport(cl, ls(), envir = environment())
    start.yearZ <- parLapply(cl, start.yearZ, as.numeric)
    opt <- parLapply(cl, start.yearZ, function(x) try(optim(startZ, fx, yearZ = x, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                                                            CPUE = tmb.dat$CPUE, LH = tmb.dat$LH, Lc = tmb.dat$Lc,
                                                            nbreaks = tmb.dat$nbreaks, loglikeCPUE = tmb.dat$loglikeCPUE,
                                                            method = "BFGS", control = list(maxit = 1e7))))
    data.pred <- clusterMap(cl, function(x, y) if(!inherits(x, "try-error")) gx(c(getElement(x, "par"), y), Lbar = tmb.dat$Lbar,
                                                      ss = tmb.dat$ss, CPUE = tmb.dat$CPUE,
                                                      LH = tmb.dat$LH, Lc = tmb.dat$Lc,
                                                      nbreaks = tmb.dat$nbreaks, loglikeCPUE = tmb.dat$loglikeCPUE),
                            x = opt, y = start.yearZ)
    stopCluster(cl)
    nll.vec <- vapply(opt, returnNAoptvalue, numeric(1))
    nll.data <- t(vapply(data.pred, function(x) if(!is.null(x)) -x[[5]] else as.numeric(rep(NA, 2)), numeric(2)))
  }

  output <- as.data.frame(ygrid + MLZ_data@Year[1] - 1)
  output <- cbind(output, nll.vec, nll.data)
  names(output) <- c(paste0("Year", 1:ncp), "negLL", "ML", "CR")
  rownames(output) <- 1:nrow(output)

  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(list = old_par), add = TRUE)

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
      if(color) plot_mlcr_2ncp_filled.contour(output)
      if(!color) {
        layout(matrix(c(1,2,3,3,3,3), nrow = 2))
        plot_mlcr_2ncp_contour(output, value.var = "ML")
        plot_mlcr_2ncp_contour(output, value.var = "CR")
        plot_mlcr_2ncp_contour(output, value.var = "negLL")
      }
    }
  }
  return(output)
}

#' Grid search for the multispecies mean length estimator
#'
#' A grid search is performed over the time series. Can be used to identify local and global minima.
#'
#' @param MLZ.list A list containing an object of class \code{MLZ_data} for each species or stock.
#' @param ncp The number of change points.
#' @param model The name of the multispecies model for the grid search.
#' @param stZ The starting value of total mortality rate used in the grid search.
#' @param parallel Whether the grid search is performed with parallel processing.
#' @param min.time The minimum number of years between change points. Only used if \code{ncp > 1}.
#' @param figure If \code{TRUE}, creates a plot of the likelihood over the grid search. Only used
#' if \code{ncp = 1} or \code{2}.
#' @param color If \code{TRUE}, creates a color plot for the likelihood surface. Only used if
#' \code{ncp = 2}.
#' @examples
#' data(PRSnapper)
#' profile_MLmulti(PRSnapper, ncp = 1, model = "MSM1")
#'
#' @export
profile_MLmulti <- function(MLZ.list, ncp, model = c("SSM", "MSM1", "MSM2", "MSM3"),
                          stZ = 0.5, parallel = FALSE, min.time = 3, figure = TRUE, color = TRUE) {

  model <- match.arg(model)

  ncp <- as.integer(ncp)
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

  tmb.dat <- list(LH = LH, Lc = Lc, nbreaks = ncp,
                  Lbar = Lbar, ss = ss, nspec = nspec)
  tmb.dat$Lbar[is.na(tmb.dat$ss) | tmb.dat$ss == 0] <- -99
  tmb.dat$ss[is.na(tmb.dat$ss) | tmb.dat$Lbar < 0] <- 0

  if("SSM" %in% model || "MSM1" %in% model) pr <- .SSM_MSM1_profile(tmb.dat, stZ, parallel, ygrid)
  if("MSM2" %in% model || "MSM3" %in% model) pr <- .MSM23_profile(tmb.dat, stZ, parallel, ygrid, model)

  output <- as.data.frame(ygrid + min(years) - 1)
  output <- cbind(output, pr$negLL, pr$loglikesp)
  names(output) <- c(paste0("Year", 1:ncp), "negLL", paste0("Species_", 1:nspec))
  rownames(output) <- 1:nrow(output)

  if(figure) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(list = old_par), add = TRUE)

    sp.name <- lapply(MLZ.list, getElement, "Stock")
    sp.ind <- vapply(sp.name, length, numeric(1))
    num <- which(sp.ind == 0L)
    sp.name[num] <- paste0("Species_", num)

    if(ncp == 1) {
      color.vec <- rich.colors(nspec)
      loglikesp <- pr$loglikesp
      nll.vec2 <- output$negLL - min(output$negLL)
      min.loglikesp <- apply(loglikesp, 2, min)
      loglikesp2 <- t(t(loglikesp) - min.loglikesp)

      plot(output$Year1, nll.vec2, xlab = "Change Point", ylim = c(0, ceiling(max(cbind(nll.vec2, loglikesp2)))),
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
        if(color) filled.contour(x = x.vec, y = y.vec, z = z.matrix,
                                 xlab = "First Change Point", ylab = "Second Change Point",
                                 plot.axes = {
                                   axis(1);
                                   axis(2);
                                   lines(x.vec, y.vec);
                                   points(output[min.index, 1], output[min.index, 2], pch = 16, cex = 1.5);
                                   segments(x0 = x.vec[1], x1 = output[min.index, 1], y0 = output[min.index, 2], lty = 2, lwd = 2);
                                   segments(x0 = output[min.index, 1], y0 = y.vec[1], y1 = output[min.index, 2], lty = 2, lwd = 2);
                                 }
        )
        if(!color) {
          contour(x = x.vec, y = y.vec, z = z.matrix,
                  xlab = "First Change Point", ylab = "Second Change Point", labcex = 1, las = 1)
          lines(x.vec, y.vec)
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


