plotPcaFd <- function(pcafd, nx = 128, pointplot = TRUE, harm = 0,
                        expand = 0, cycle = FALSE, ...)
{
#
#  Plots the harmonics produced by PCA.FD.
#
#   If pointplot=T, then the harmonics are plotted as + and -
#    otherwise lines are used.  Another thing that needs doing is an
#     arrowplot option.
#
# If harm = 0 (the default) then all the computed harmonics are plotted.
#   Otherwise those in jharm are plotted.
# If expand =0 then effect of +/- 2 standard deviations of each pc are given
#   otherwise the factor expand is used.
# If cycle=T and there are 2 variables then a cycle plot will be drawn
#  If the number of variables is anything else, cycle will be ignored.
#

#  Last modified 27 June 2001

  if (!(inherits(pcafd, "pca.fd"))) stop('Argument PCAFD is not a pca.fd object.')

  harmfd  <- pcafd[[1]]
  basisfd <- getbasis(harmfd)
  rangex  <- basisfd$rangeval
  x <- seq(rangex[1], rangex[2], length = nx)
  fdmat   <- eval.fd(x, harmfd)
  meanmat <- eval.fd(x, pcafd$meanfd)
  dimfd   <- dim(fdmat)
  nharm   <- dimfd[2]
  harm <- as.vector(harm)
  if(harm[1] == 0) harm <- (1:nharm)
  if(length(dimfd) == 2) {
    for(iharm in harm) {
      if(expand == 0) fac <- sqrt(pcafd$values[iharm]) else fac <- expand
      vecharm <- fdmat[, iharm]
      pcmat <- cbind(meanmat + fac * vecharm, meanmat - fac * vecharm)
      if (pointplot) plottype <- "p" else plottype <- "l"
      percentvar <- round(100 * pcafd$varprop[iharm], 1)
      matplot(x, pcmat, lty = 2:3, pch = "+-",
              sub = paste("PCA function", iharm,
                          "(Percentage of variability", percentvar, ")"),
              col = 2:3, type = plottype, ...)
      lines(x, meanmat)
      mtext("Click to advance to next plot", side = 3, line = -3, outer = TRUE)
      text(locator(1), "")
    }
  } else {
    if(cycle && dimfd[3] == 2) {
      meanmat <- drop(meanmat)
      for(iharm in harm) {
        if(expand == 0) {
          fac <- 2 * sqrt(pcafd$values[iharm])
        } else {
          fac <- expand
        }
        matharm <- fdmat[, iharm,  ]
        mat1 <- meanmat + fac * matharm
        mat2 <- meanmat - fac * matharm
        if (pointplot) plottype <- "p" else plottype <- "l"
        percentvar <- round(100 * pcafd$varprop[iharm],1)
        matplot(cbind(mat1[, 1], mat2[, 1]),
                cbind(mat1[, 2], mat2[, 2]), lty = 2:3, pch = "+-",
                sub = paste("PCA function", iharm,
                            "(Percentage of variability", percentvar, ")"),
                            col = 2:3, type = plottype, ...)
        lines(meanmat)
        mtext("Click to advance to next plot",
              side = 3, line = -3, outer = TRUE)
        text(locator(1), "")
      }
    } else {
      for (iharm in harm) {
        if (expand == 0) fac <- sqrt(pcafd$values[iharm]) else fac <- expand
        meanmat <- drop(meanmat)
        matharm <- fdmat[, iharm,  ]
        nvar    <- dim(matharm)[2]
        for (jvar in 1:nvar) {
          pcmat <- cbind(meanmat[, jvar] + fac * matharm[, jvar],
                         meanmat[, jvar] - fac * matharm[, jvar])
          if (pointplot) plottype <- "p" else plottype <- "l"
          percentvar <- round(100 * pcafd$varprop[iharm], 1)
          matplot(x, pcmat, lty = 2:3, pch = "+-",
                  sub = paste("PCA function", iharm,
                              "(Percentage of variability", percentvar,")"),
                  main = dimnames(fdmat)[[3]][jvar],
                  col = 2:3, type = plottype, ...)
          lines(x, meanmat[, jvar])
        }
        mtext("Click to advance to next set of plots",
              side = 3, line = -3, outer = TRUE)
        text(locator(1), "")
      }
    }
  }
  invisible()
}
