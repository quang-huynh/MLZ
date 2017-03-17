.SSM_MSM1_est <- function(tmb.dat, start, model, spCont) {
  start <- lapply(start, as.numeric)
  stpar <- as.numeric(do.call(c, start))
  ncp <- tmb.dat$nbreaks
  nspec <- tmb.dat$nspec
  if(model == "MSM1") MSM1.switch <- 1L else MSM1.switch <- 0L

  opt <- optim(stpar, SSM_MSM1_negLL, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
               LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
               isMSM1 = MSM1.switch, spCont = spCont, 
               method = "BFGS", control = list(maxit = 1e7), hessian = T)
  data.pred <- SSM_MSM1_pred(opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                             LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                             isMSM1 = MSM1.switch, spCont = spCont)
  opt$Lpred <- data.pred[[1]]
  sigmaL <- data.pred[[2]]
  opt$par <- c(opt$par, sigmaL)
  MLmhessian <- hessian(SSM_MSM1_fullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                        LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                        isMSM1 = MSM1.switch, spCont = spCont)
  covariance <- solve(MLmhessian)
  opt$corr <- cov2cor(covariance)
  opt$gradient <- grad(SSM_MSM1_fullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                       LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                       isMSM1 = MSM1.switch, spCont = spCont)

  results.matrix <- matrix(c(opt$par, sqrt(diag(covariance))), ncol = 2)
  Z.name <- paste0(rep(paste0("Z[",1:(ncp+1)), nspec), rep(paste0(",", 1:nspec, "]"), each = ncp+1))
  if(model == "SSM") yearZ.name <- paste0(rep(paste0("yearZ[", 1:ncp), nspec), rep(paste0(",", 1:nspec, "]"), each = ncp))
  if(model == "MSM1") yearZ.name <- paste0("yearZ[", 1:ncp, "]")
  sigma.name <- paste0("sigma[", 1:nspec, "]")

  rownames(results.matrix) <- rownames(opt$corr) <- colnames(opt$corr) <- c(Z.name, yearZ.name, sigma.name)
  colnames(results.matrix) <- c("Estimate", "Std. Error")
  return(list(opt = opt, results.matrix = results.matrix))
}

.MSM23_est <- function(tmb.dat, start, model, spCont) {
  start <- lapply(start, as.numeric)
  stpar <- as.numeric(do.call(c, start))
  ncp <- tmb.dat$nbreaks
  nspec <- tmb.dat$nspec
  if(model == "MSM3") MSM3.switch <- 1L else MSM3.switch <- 0L

  opt <- optim(stpar, MSM23_negLL, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
               LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
               isMSM3 = MSM3.switch, spCont = spCont, method = "BFGS", control = list(maxit = 1e7))
  data.pred <- MSM23_pred(opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                             LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                             isMSM3 = MSM3.switch, spCont = spCont)
  opt$Lpred <- data.pred[[1]]
  sigmaL <- data.pred[[2]]
  opt$par <- c(opt$par, sigmaL)
  MLmhessian <- hessian(MSM23_fullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                        LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                        isMSM3 = MSM3.switch, spCont = spCont)
  covariance <- solve(MLmhessian)
  opt$corr <- cov2cor(covariance)
  opt$gradient <- grad(MSM23_fullnegLL, opt$par, Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                       LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                       isMSM3 = MSM3.switch, spCont = spCont)

  Z.name <- paste0(rep(paste0("Z[",1:(ncp+1)), nspec), rep(paste0(",", 1:nspec, "]"), each = ncp+1))
  Z1.name <- paste0("Z1[", 1:nspec, "]")
  delta.name <- paste0("delta[", 1:ncp, ",1]")
  yearZ.name <- paste0("yearZ[", 1:ncp, "]")
  sigma.name <- paste0("sigma[", 1:nspec, "]")
  if(model == "MSM2") {
    eps.name <- paste0("epsilon[", 2:nspec, "]")
    rownames(opt$corr) <- colnames(opt$corr) <- c(Z1.name, delta.name, eps.name, yearZ.name, sigma.name)
    Z.var <- deltamethod_MSM2(Z = data.pred[[4]], delta = opt$par[(nspec+1):(nspec+ncp)],
                              epsilon = opt$par[(nspec+ncp+1):(nspec+ncp+nspec-1)],
                              M = tmb.dat$LH[, 3], cov = covariance[1:(nspec+ncp+nspec-1), 1:(nspec+ncp+nspec-1)])
    results.rownames <- c(Z.name, delta.name, eps.name, yearZ.name, sigma.name)
  } else {
    rownames(opt$corr) <- colnames(opt$corr) <- c(Z1.name, delta.name, yearZ.name, sigma.name)
    Z.var <- deltamethod_MSM3(Z = data.pred[[4]], delta = opt$par[(nspec+1):(nspec+ncp)],
                              M = tmb.dat$LH[, 3], cov = covariance[1:(nspec+ncp), 1:(nspec+ncp)])
    results.rownames <- c(Z.name, delta.name, yearZ.name, sigma.name)
  }
  est <- c(as.numeric(t(data.pred[[4]])), opt$par[(nspec+1):length(opt$par)])
  Z.allvar <- cbind(diag(covariance)[1:nspec], Z.var)
  var <- c(as.numeric(t(Z.allvar)), diag(covariance)[(nspec+1):length(opt$par)])
  results.matrix <- matrix(c(est, sqrt(var)), ncol = 2)
  rownames(results.matrix) <- results.rownames
  colnames(results.matrix) <- c("Estimate", "Std. Error")

  return(list(opt = opt, results.matrix = results.matrix))
}

deltamethod_MSM2 <- function(Z, delta, epsilon, M, cov) {
  ncp <- ncol(Z) - 1
  nspec <- nrow(Z)

  VarZ <- matrix(NA, nrow = nspec, ncol = ncp)
  cov.temp1 <- matrix(0, ncol = 4, nrow = 4)
  cov.temp2 <- matrix(0, ncol = 5, nrow = 5)

  for(i in 1:ncp) {
    ii <- nspec + i
    for(n in 1:nspec) {
      nn <- nspec + ncp + n - 1
      if(n == 1) {
        cov.temp <- matrix(0, ncol = 3, nrow = 3)
        cov.temp[1:2, 1:2] <- cov[c(n, ii), c(n, ii)]
        VarZ[n, i] <- deltamethod(~ x2 * (x1 - x3) + x3, mean = c(Z[n, i], delta[i], M[n]),
                                  cov = cov.temp, ses = FALSE)
        if(i < ncp) {
          cov.temp1[1:3, 1:3] <- cov[c(n, ii, ii+1), c(n, ii, ii+1)]
          zz <- deltamethod(list(~ x2 * (x1 - x4) + x4, ~ x3),
                            mean = c(Z[n, i], delta[i:(i+1)], M[n]),
                            cov = cov.temp1, ses = FALSE)
          cov[n,ii+1] <- cov[ii+1,n] <- zz[1,2]
        }
      }
      if(n > 1) {
        cov.temp <- matrix(0, ncol = 4, nrow = 4)
        cov.temp[1:3, 1:3] <- cov[c(n, ii, nn), c(n, ii, nn)]
        VarZ[n, i] <- deltamethod(~ x2 * x3 * (x1 - x4) + x4,
                                  mean = c(Z[n, i], delta[i], epsilon[n-1], M[n]),
                                  cov = cov.temp, ses = FALSE)
        if(i < ncp) {
          cov.temp2[1:4, 1:4] <- cov[c(n, ii, ii+1, nn), c(n, ii, ii+1, nn)]
          zz <- deltamethod(list(~ x2 * x4 * (x1 - x5) + x5, ~ x3, ~ x4),
                            mean = c(Z[n, i], delta[i:(i+1)], epsilon[n-1], M[n]),
                            cov = cov.temp2, ses = FALSE)
          cov[n,ii+1] <- cov[ii+1,n] <- zz[1,2]
          cov[n,nn] <- cov[nn,n] <- zz[1,3]
          cov[ii+1,nn] <- cov[ii+1,nn] <- zz[2,3]
        }
      }
    }
    if(i < ncp) diag(cov[1:nspec, 1:nspec]) <- VarZ[, i]
  }
  return(VarZ)
}

deltamethod_MSM3 <- function(Z, delta, M, cov) {
  ncp <- ncol(Z) - 1
  nspec <- nrow(Z)

  VarZ <- matrix(NA, nrow = nspec, ncol = ncp)
  cov.temp <- matrix(0, ncol = 3, nrow = 3)
  cov.temp2 <- matrix(0, ncol = 4, nrow = 4)

  for(i in 1:ncp) {
    ii <- nspec + i
    for(n in 1:nspec) {
      cov.temp[1:2, 1:2] <- cov[c(n, ii), c(n, ii)]
      VarZ[n, i] <- deltamethod(~ x2 * (x1 - x3) + x3, mean = c(Z[n, i], delta[i], M[n]),
                                cov = cov.temp, ses = FALSE)
      if(i < ncp) {
        cov.temp2[1:3, 1:3] <- cov[c(n, ii, ii+1), c(n, ii, ii+1)]
        zz <- deltamethod(list(~ x2 * (x1 - x4) + x4, ~ x3), mean = c(Z[n, i], delta[i:(i+1)], M[n]),
                          cov = cov.temp2, ses = FALSE)
        cov[n,ii+1] <- cov[ii+1,n] <- zz[1,2]
      }
    }
    if(i < ncp) diag(cov[1:nspec, 1:nspec]) <- VarZ[, i]
  }
  return(VarZ)
}


.SSM_MSM1_profile <- function(tmb.dat, stZ, parallel, ygrid, spCont) {
  nspec <- tmb.dat$nspec
  ncp <- tmb.dat$nbreaks
  startZ <- rep(stZ, nspec*(ncp + 1))

  if(!parallel) {
    loglikesp <- matrix(0, nrow = nrow(ygrid), ncol = nspec)
    nll.vec <- numeric(nrow(ygrid))

    for(i in 1:nrow(ygrid)) {
      opt <- try(optim(startZ, SSM_MSM1_profile, year = as.numeric(ygrid[i, ]),
                       Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                       LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                       spCont = spCont,
                       method = "BFGS", control = list(maxit = 1e7)))
      if(inherits(opt, "try-error")) {
        nll.vec[i] <- NA
        loglikesp[i, ] <- rep(NA, nspec)
      } else {
        nll.vec[i] <- opt$value
        data.pred <- SSM_MSM1_pred(c(opt$par, as.numeric(ygrid[i, ])),
                                   Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                                   LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                                   isMSM1 = 1L, spCont = spCont)
        loglikesp[i, ] <- -1 * data.pred[[3]]
      }
    }
  }

  if(parallel) {
    start.yearZ <- split(ygrid, 1:nrow(ygrid))
    cl <- makeCluster(detectCores())
    clusterExport(cl, ls(), envir = environment())
    start.yearZ <- parLapply(cl, start.yearZ, as.numeric)
    opt <- parLapply(cl, start.yearZ, function(x) try(optim(startZ, SSM_MSM1_profile, year = x,
                                                            Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                                                            LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                                                            spCont = spCont, method = "BFGS", control = list(maxit = 1e7))))
    data.pred <- clusterMap(cl, function(x, y) if(!inherits(x, "try-error")) SSM_MSM1_pred(c(x$par, y),
                                                                                           Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                                                                                           LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                                                                                           isMSM1 = 1L, spCont),
                            x = opt, y = start.yearZ)
    stopCluster(cl)
    nll.vec <- vapply(opt, returnNAoptvalue, numeric(1))
    loglikesp <- t(vapply(data.pred, function(x) if(!is.null(x)) -x[[3]] else as.numeric(rep(NA, nspec)), numeric(nspec)))
  }

  return(list(negLL = nll.vec, loglikesp = loglikesp))
}


.MSM23_profile <- function(tmb.dat, stZ, parallel, ygrid, model, spCont) {
  nspec <- tmb.dat$nspec
  ncp <- tmb.dat$nbreaks
  startZ1 <- rep(stZ, nspec)
  delta <- rep(1, ncp)

  if("MSM2" %in% model) {
    MSM3.switch <- 0L
    epsilon <- rep(1, nspec - 1)
    stpar <- c(startZ1, delta, epsilon)
  }
  if("MSM3" %in% model){
    MSM3.switch <- 1L
    stpar <- c(startZ1, delta)
  }

  if(!parallel) {
    loglikesp <- matrix(0, nrow = nrow(ygrid), ncol = nspec)
    nll.vec <- numeric(nrow(ygrid))

    for(i in 1:nrow(ygrid)) {
      opt <- try(optim(stpar, MSM23_profile, year = as.numeric(ygrid[i, ]),
                       Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                       LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks, isMSM3 = MSM3.switch,
                       spCont = spCont, method = "BFGS", control = list(maxit = 1e7)))
      if(inherits(opt, "try-error")) {
        nll.vec[i] <- NA
        loglikesp[i, ] <- rep(NA, nspec)
      } else {
        nll.vec[i] <- opt$value
        data.pred <- MSM23_pred(c(opt$par, as.numeric(ygrid[i, ])),
                                   Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                                   LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                                   isMSM3 = MSM3.switch, spCont = spCont)
        loglikesp[i, ] <- -1 * data.pred[[3]]
      }
    }
  }

  if(parallel) {
    start.yearZ <- split(ygrid, 1:nrow(ygrid))
    cl <- makeCluster(detectCores())
    clusterExport(cl, ls(), envir = environment())
    start.yearZ <- parLapply(cl, start.yearZ, as.numeric)
    opt <- parLapply(cl, start.yearZ, function(x) try(optim(stpar, MSM23_profile, year = x,
                                                            Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                                                            LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                                                            isMSM3 = MSM3.switch, spCont = spCont,
                                                            method = "BFGS", control = list(maxit = 1e7))))
    data.pred <- clusterMap(cl, function(x, y) if(!inherits(x, "try-error")) MSM23_pred(c(x$par, y),
                                                                                        Lbar = tmb.dat$Lbar, ss = tmb.dat$ss,
                                                                                        LH = tmb.dat$LH, Lc = tmb.dat$Lc, nbreaks = tmb.dat$nbreaks,
                                                                                        isMSM3 = MSM3.switch, spCont = spCont),
                            x = opt, y = start.yearZ)
    stopCluster(cl)
    nll.vec <- vapply(opt, returnNAoptvalue, numeric(1))
    loglikesp <- t(vapply(data.pred, function(x) if(!is.null(x)) -x[[3]] else as.numeric(rep(NA, nspec)), numeric(nspec)))
  }

  return(list(negLL = nll.vec, loglikesp = loglikesp))
}





plot_mlcr_2ncp_contour <- function(output, value.var) {
  z.matrix <- acast(output, Year1 ~ Year2, value.var = value.var)
  output.vec <- getElement(output, value.var)
  min.index <- which.min(output.vec)
  x.vec <- as.numeric(rownames(z.matrix))
  y.vec <- as.numeric(colnames(z.matrix))

  contour(x = x.vec, y = y.vec, z = z.matrix,
          xlab = "First Change Point", ylab = "Second Change Point", labcex = 1, las = 1)
  lines(x.vec, y.vec)
  points(output[min.index, 1], output[min.index, 2], pch = 16, cex = 1.5)
  segments(x0 = x.vec[1], x1 = output[min.index, 1], y0 = output[min.index, 2], lty = 2, lwd = 2)
  segments(x0 = output[min.index, 1], y0 = y.vec[1], y1 = output[min.index, 2], lty = 2, lwd = 2)

  if(value.var == "ML") title("Negative log-likelihood surface for Mean Length")
  if(value.var == "CR") title("Negative log-likelihood surface for CPUE")
  if(value.var == "negLL") title("Total negative log-likelihood surface")

  invisible()

}


plot_mlcr_2ncp_filled.contour <- function(output) {
  new.output <- output[3:5]
  mins <- apply(new.output, 2, min, na.rm = TRUE)
  new.output <- t(t(new.output) - mins)
  z.range <- range(new.output, finite = TRUE)

  output.tmp <- cbind(output[, 1:2], ML = new.output[, 2])
  z.matrix <- acast(output.tmp, Year1 ~ Year2, value.var = "ML")

  output.vec <- getElement(output.tmp, "ML")
  min.index <- which.min(output.vec)
  x.vec <- as.numeric(rownames(z.matrix))
  y.vec <- as.numeric(colnames(z.matrix))

  layout(matrix(c(1,2,3,3,3,3,4,4), nrow = 2), widths = c(1, 1, 0.8, 0.25))

  filled.contour3(x = x.vec, y = y.vec, z = z.matrix, zlim = z.range,
                  xlab = "First Change Point", ylab = "Second Change Point",
                  plot.axes = {
                    axis(1);
                    axis(2);
                    lines(x.vec, y.vec);
                    points(output.tmp[min.index, 1], output.tmp[min.index, 2], pch = 16, cex = 1.5);
                    segments(x0 = x.vec[1], x1 = output.tmp[min.index, 1], y0 = output.tmp[min.index, 2], lty = 2, lwd = 2);
                    segments(x0 = output.tmp[min.index, 1], y0 = y.vec[1], y1 = output.tmp[min.index, 2], lty = 2, lwd = 2)
                  }
  )
  title("Negative log-likelihood for Mean Length")

  output.tmp <- cbind(output[, 1:2], CR = new.output[, 3])
  z.matrix <- acast(output.tmp, Year1 ~ Year2, value.var = "CR")

  output.vec <- getElement(output.tmp, "CR")
  min.index <- which.min(output.vec)
  x.vec <- as.numeric(rownames(z.matrix))
  y.vec <- as.numeric(colnames(z.matrix))

  filled.contour3(x = x.vec, y = y.vec, z = z.matrix,
                  xlab = "First Change Point", ylab = "Second Change Point",
                  plot.axes = {
                    axis(1);
                    axis(2);
                    lines(x.vec, y.vec);
                    points(output.tmp[min.index, 1], output.tmp[min.index, 2], pch = 16, cex = 1.5);
                    segments(x0 = x.vec[1], x1 = output.tmp[min.index, 1], y0 = output.tmp[min.index, 2], lty = 2, lwd = 2);
                    segments(x0 = output.tmp[min.index, 1], y0 = y.vec[1], y1 = output.tmp[min.index, 2], lty = 2, lwd = 2)
                  }
  )
  title("Negative log-likelihood for CPUE")


  output.tmp <- cbind(output[, 1:2], negLL = new.output[, 1])
  z.matrix <- acast(output.tmp, Year1 ~ Year2, value.var = "negLL")

  output.vec <- getElement(output.tmp, "negLL")
  min.index <- which.min(output.vec)
  x.vec <- as.numeric(rownames(z.matrix))
  y.vec <- as.numeric(colnames(z.matrix))

  filled.contour3(x = x.vec, y = y.vec, z = z.matrix, zlim = z.range,
                  xlab = "First Change Point", ylab = "Second Change Point",
                  plot.axes = {
                    axis(1);
                    axis(2);
                    lines(x.vec, y.vec);
                    points(output.tmp[min.index, 1], output.tmp[min.index, 2], pch = 16, cex = 1.5);
                    segments(x0 = x.vec[1], x1 = output.tmp[min.index, 1], y0 = output.tmp[min.index, 2], lty = 2, lwd = 2);
                    segments(x0 = output.tmp[min.index, 1], y0 = y.vec[1], y1 = output.tmp[min.index, 2], lty = 2, lwd = 2)
                  }
  )
  title("Total negative log-likelihood surface")

  filled.legend(z = z.matrix, zlim = z.range)


  invisible()

}


# Following code obtained from: http://wiki.cbr.washington.edu/qerm/index.php/R/Contour_Plots
# Thanks QERM!

# modification by Ian Taylor of the filled.contour function
# to remove the key and facilitate overplotting with contour()
# further modified by Carey McGilliard and Bridget Ferris
# to allow multiple plots on one page

filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
            col = color.palette(length(levels) - 1), plot.title, plot.axes,
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
            axes = TRUE, frame.plot = axes,mar, ...)
  {

    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")
    par(las = las)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) stop("no proper 'z' matrix specified")
    if (!is.double(z)) storage.mode(z) <- "double"

    .filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col)

    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot)
      box()
    if (missing(plot.title))
      title(...)
    else plot.title
    invisible()
  }

filled.legend <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
            col = color.palette(length(levels) - 1), plot.title, plot.axes,
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
            axes = TRUE, frame.plot = axes, ...)
  {
    # modification of filled.contour by Carey McGilliard and Bridget Ferris
    # designed to just plot the legend
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")

    mar.orig <- par()$mar
    on.exit(par(mar = mar.orig))

    mar <- mar.orig
    mar[2] <- 0
    par(las = las, mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
      if (axes)
        axis(4)
    }
    else key.axes
    box()
  }
