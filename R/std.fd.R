std.fd <- function(fd)
{
  #  Compute the standard deviation functions for functional observations
  #  Argument:
  #  FD    ... a functional data object
  #  Return:
  #  STDFD ... a functional data for the standard deviation functions

  #  Last modified 6 Feb 2001

  if (!(inherits(fd, "fd"))) stop("Argument  FD not a functional data object.")

  coef     <- getcoef(fd)
  coefd    <- dim(coef)
  ndim     <- length(coefd)
  if (coefd[1] == 1) stop("Only one replication found.")

  basisfd  <- getbasis(fd)
  nbasis   <- basisfd$nbasis
  rangeval <- basisfd$rangeval
  fdnames  <- getnames(fd)

  varbifd  <- var.fd(fd)

  neval    <- 10*nbasis + 1
  evalarg  <- seq(rangeval[1], rangeval[2], length=neval)
  vararray <- eval.bifd(evalarg, evalarg, varbifd)
  nvdim    <- length(dim(vararray))

  if (ndim == 2) {
    stdmat  <- matrix(sqrt(diag(vararray)), neval, 1)
  } else {
    nvar <- coefd[3]
    stdmat <- matrix(0, neval, nvar)
    m <- 0
    for (j in 1:nvar) {
      m <- m + j
      if (nvdim == 3) {
        stdmat[,j] <- sqrt(diag(varray[,,1,m]))
      } else {
        stdmat[,j] <- sqrt(diag(varray[,,m]))
      }
    }
  }
  stdcoef <- project.basis(stdmat, evalarg, basisfd)
  names(fdnames)[2] <- "Std. Dev."
  names(fdnames)[3] <- paste("Std. Dev.",names(fdnames)[3])
  stdfd <- create.fd(stdcoef, basisfd, fdnames)
  return(stdfd)
}
