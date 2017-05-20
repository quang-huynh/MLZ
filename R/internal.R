
plot_mlcr_2ncp_contour <- function(output, min.time, value.var) {
  z.matrix <- acast(output, Year1 ~ Year2, value.var = value.var)
  output.vec <- getElement(output, value.var)
  min.index <- which.min(output.vec)
  x.vec <- as.numeric(rownames(z.matrix))
  y.vec <- as.numeric(colnames(z.matrix))

  contour(x = x.vec, y = y.vec, z = z.matrix,
          xlab = "First Change Point", ylab = "Second Change Point", labcex = 1, las = 1)
  lines(x.vec, x.vec + min.time)
  points(output[min.index, 1], output[min.index, 2], pch = 16, cex = 1.5)
  segments(x0 = x.vec[1], x1 = output[min.index, 1], y0 = output[min.index, 2], lty = 2, lwd = 2)
  segments(x0 = output[min.index, 1], y0 = y.vec[1], y1 = output[min.index, 2], lty = 2, lwd = 2)

  if(value.var == "ML") title("Negative log-likelihood surface for Mean Length")
  if(value.var == "CR") title("Negative log-likelihood surface for CPUE")
  if(value.var == "negLL") title("Total negative log-likelihood surface")

  invisible()

}


plot_mlcr_2ncp_filled.contour <- function(output, min.time) {
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
                    lines(x.vec, x.vec + min.time);
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
                    lines(x.vec, x.vec + min.time);
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
                    lines(x.vec, x.vec + min.time);
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
