plotCcaFd <- function(ccafd, overplt = FALSE, jcan = 0, flip = FALSE, ...)
{
#  Plot a functional canonical correlation analysis object cca.fd
#
#  If overplt=T  then each pair of weight functions is plotted in
#     a single plot.  The line types and colours of the
#     "x" and "y" curves respectively are specified as in plotFd   .
#  If overplt=F  then the weight functions are plotted in separate
#     plots, side by side if a command like par(mfrow=c(2,2)) is
#       used.
#
#  If jcan=0, then all the pairs of variates are plotted.  Otherwise
#     only the variates jcan are plotted (eg if jcan=1, only the leading
#     variate is plotted, if jcan=c(1,3) only the first and third.)
#
#  If flip[j] is T then the jth pair of weight functions are multiplied
#     by -1.  If flip is a scalar it is replicated to the necessary length.
#
#  Other arguments are passed to plotFd
#

#  Last modified 6 Feb 2001

  if (!(inherits(ccafd, "cca.fd"))) stop("First argument not of CCA.FD class.")

  wtfd   <- ccafd[[1]]
  wtcoef <- getcoef(wtfd)
  if (jcan[1] != 0) wtcoef <- wtcoef[, jcan,  , drop = FALSE]
  wtcoef <- aperm((-1)^flip * aperm(wtcoef, c(3, 2, 1)), c(3, 2, 1))
  if (overplt) {
    wtcoef <- aperm(wtcoef, c(1, 3, 2))
    wtfd[[1]] <- wtcoef
    templabs <- wtfd$fdnames[[3]]
    wtfd$fdnames[[3]] <- wtfd$fdnames[[2]]
    wtfd$fdnames[[2]] <- templabs
    plot(wtfd, ylab = "Weight function", ...)
  } else {
    ncan <- dim(wtcoef)[2]
    for (jj in (1:ncan)) {
      wtfdtemp <- wtfd
      wtfdtemp[[1]] <- wtfdtemp[[1]][, jj,  , drop = FALSE]
      plot(wtfdtemp, ylab = "Weight function",
           sub = dimnames(wtfdtemp[[1]])[[2]], ...)
    }
  }
  invisible()
}
