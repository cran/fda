derivFd <- function(fd, Lfd=1)
{
  #  Applies linear differential operator LFD to functional data object FD
  #    and returns the result as functional data object DERIVFD.

  #  Last modified 6 Feb 2001

  if (!(inherits(fd, "fd"))) stop("Argument  FD not a functional data object.")

  basisfd  <- getbasis(fd)
  nbasis   <- basisfd$nbasis
  rangeval <- basisfd$rangeval

  evalarg  <- seq(rangeval[1], rangeval[2], len=10*nbasis+1)
  Lfdmat   <- eval.fd(evalarg, fd, Lfd)

  Lfdcoef  <- project.basis(Lfdmat, evalarg, basisfd)

  Dfdnames <- fd$fdnames
  Dfdnames[[3]] <- paste("D",Dfdnames[[3]])

  derivfd <- create.fd(Lfdcoef, basisfd, Dfdnames)

  return(derivfd)
}
