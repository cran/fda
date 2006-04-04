lines.fd <- function(x, Lfdobj=int2Lfd(0), ...)
{
  fdobj <- x
  #  Plot a functional data object FD using lines in a pre-existing plot.
  #  If there are multiple variables, each curve will appear in the same plot.
  #  The remaining optional arguments are the same as those available
  #     in the regular "lines" function.

  #  Last modified 1 October 2005

  if (!inherits(fdobj,  "fd"))  stop(
		"First argument is not a functional data object.")
  if (!inherits(Lfdobj, "Lfd")) stop(
      "Second argument is not a linear differential operator.")

  coef   <- fdobj$coefs
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  nbasis <- coefd[1]
  nrep   <- coefd[2]
  if (ndim > 3) nvar <- coefd[3] else nvar <- 1
  crvnames <- fdobj$fdnames[[2]]
  varnames <- fdobj$fdnames[[3]]

  basisobj <- fdobj$basis
  rangex   <- basisobj$rangeval
  x        <- seq(rangex[1],rangex[2],length=101)
  fdmat    <- eval.fd(x,fdobj,Lfdobj)

  if (length(dim(coef)) < 2) {
    lines (x,fdmat,...)
  }
  if (length(dim(coef)) ==2 ) {
    matlines (x,fdmat,...)
  }
  if (length(dim(coef)) == 3) {
    for (ivar in 1:nvar) {
      matlines (x,fdmat[,,ivar],type="l",lty=1,
                main=varnames[ivar],...)
    }
  }
  invisible()
}
