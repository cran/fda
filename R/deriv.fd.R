#. added by Andrea Gilardi 29 August 24

deriv <- function(expr, ...) 

  UseMethod('deriv')

deriv.fd <- function(expr, Lfdobj=int2Lfd(1), ...)
{
  #  Applies linear differential operator LFD to functional data object FD
  #    and returns the result as functional data object DERIVFD.

  #  Last modified 6 January 2020 by Jim Ramsay

  fdobj <- expr
  if (!inherits(fdobj, "fd")) stop(
		"Argument  FD not a functional data object.")

  Lfdobj   <- int2Lfd(Lfdobj)

  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  rangeval <- basisobj$rangeval

  evalarg  <- seq(rangeval[1], rangeval[2], len=10*nbasis+1)
  Lfdmat   <- eval.fd(evalarg, fdobj, Lfdobj)

  Lfdcoef  <- project.basis(Lfdmat, evalarg, basisobj)

  Dfdnames <- fdobj$fdnames
  Dfdnames[[3]] <- paste("D",Dfdnames[[3]])

  Dfdobj <- fd(Lfdcoef, basisobj, Dfdnames)

  return(Dfdobj)
}
