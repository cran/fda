smooth.fd <- function(fd, lambda = 0, Lfd = NULL, rebase = TRUE)
{
#  Smooths a functional data object.
#  Arguments for this function:
#
#  FD      ... A functional data object.
#
#  LAMBDA  ... The smoothing parameter determining the weight to be
#              placed on the size of the derivative in smoothing.  This
#              is 0 by default, but this will produce a warning
#              message that no smoothing has been carried out.
#
#  LFD     ... The order of derivative or a linear differential
#              operator to be penalized in the smoothing phase.
#              By default Lfd is set in function GETBASISPENALTY
#
#  If rebase=TRUE and the basis type is "polyg" then the basis
#    is changed to a cubic bspline  basis and before smoothing
#
#  Returns a functional data object containing a smoothed version
#    of the input functional data object
#

#  Last modified 6 Feb 2001
 
#
# Rebase to default B spline basis if rebase is T and basistype is
#    polygonal.  Then test to see if any smoothing is actually required.
#

  if (!(inherits(fd, "fd"))) stop("Argument FD not a functional data object.")

  basisfd <- getbasis(fd)
  if(rebase == TRUE && basisfd$type == "polyg") {
    fd <- data2fd(getcoef(fd), basisfd$params, fdnames = fd$fdnames)
    basisfd <- getbasis(fd)
  }
  if(lambda <= 0) {
    warning("LAMBDA was not positive. No smoothing carried out.")
    return(fd)
  }
#
#  Main smoothing step
#
  coef  <- NA * getcoef(fd)
  coefd <- dim(coef)
  ndim  <- length(coefd)
  Bmat  <- inprod(basisfd, basisfd)
#
#  set up coefficient matrix for normal equations
#
  penmat <- getbasispenalty(basisfd, Lfd)
  Cmat   <- Bmat + lambda * penmat
#
#  solve normal equations for each observation
#
  if(ndim < 3) {
    Dmat <- inprod(basisfd, fd)
    coef <- symsolve(Cmat, Dmat)
  }
  else {
    for(ivar in (1:coefd[3])) {
      Dmat <- inprod(basisfd, fd[,ivar])
      coef[,,ivar] <- symsolve(Cmat, Dmat)
    }
  }
#
#  replace coefficient matrix in fd, leaving other properties alone
#
  fd[[1]] <- coef
  return(fd)
}
