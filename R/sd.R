sd.fd <- function(fdobj)
{
  #  Compute the standard deviation functions for functional observations
  #  Argument:
  #  fdobj    ... a functional data object
  #  Return:
  #  STDFD ... a functional data for the standard deviation functions

  #  Last modified 26 October 2005

  if (!(inherits(fdobj, "fd"))) stop(
		"Argument  fdobj not a functional data object.")

  coef     <- fdobj$coefs
  coefd    <- dim(coef)
  ndim     <- length(coefd)
  if (coefd[1] == 1) stop("Only one replication found.")

  basisobj <- fdobj$basis
  fdnames  <- fdobj$fdnames
  nbasis   <- basisobj$nbasis
  rangeval <- basisobj$rangeval

  varbifd  <- var.fd(fdobj)

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
  stdcoef <- project.basis(stdmat, evalarg, basisobj)
  names(fdnames)[2] <- "Std. Dev."
  names(fdnames)[3] <- paste("Std. Dev.",names(fdnames)[3])
  stdfd <- fd(stdcoef, basisobj, fdnames)
  return(stdfd)
}

std.fd <- function(fdobj)sd.fd(fdobj)
stdev.fd <- function(fdobj)sd.fd(fdobj)
stddev.fd <- function(fdobj)sd.fd(fdobj)
