cycleplot.fd <- function(fd, matplt = TRUE, nx = 128, ...)
{

#  Performs a cycle plot of a functional data object FD,
#   assuming that FD is a bivariate function...the first component
#   of which is the x-coordinate and the second the y-coordinate
#
#  If MATPLT is T, matplot is used to plot all curves in
#     a single plot.
#  Otherwise, each curve is plotted separately, and the
#     next curve is plotted when the mouse is clicked.
#  NX is the number of sampling points to use (default 128)


 #  Last modified 6 Feb 2001

  if (!(inherits(fd, "fd"))) stop(
     "First argument is not a functional data object.")

  coef  <- getcoef(fd)
  coefd <- dim(coef)
  ndim  <- length(coefd)
  if(ndim < 3) stop("Univariate functions cannot be cycle plotted")
  nbasis <- coefd[1]
  nrep <- coefd[2]
  nvar <- coefd[3]
  basisfd  <- getbasis(fd)
  if(nvar > 2)  warning("Only first two functions used")
  crvnames <- dimnames(coef)[[2]]
  varnames <- dimnames(coef)[[3]][1:2]
  rangex   <- basisfd$rangeval
  x        <- seq(rangex[1], rangex[2], length = nx)
  fdmat    <- eval.fd(x, fd)
  fdnames  <- getnames(fd)
  crvnames <- fdnames[[2]]
  varnames <- fdnames[[3]]
  if(matplt) {
    matplot(fdmat[,  , 1], fdmat[,  , 2], type = "l", lty = 1,
            xlab=varnames[1], ylab=varnames[2], ...)
  }
  if(!matplt) {
    for (irep in 1:nrep) {
      plot(fdmat[, irep, 1], fdmat[, irep, 2], type = "l",
         lty = 1, xlab=varnames[1], ylab=varnames[2],
         main = paste("Curve", irep, crvnames[irep]), ...)
      mtext("Click to advance to next plot", side = 3,
                  line = -3, outer = TRUE)
      text(locator(1), "")
    }
  }
  invisible()
}
