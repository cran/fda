eval.monfd <- function(evalarg, Wfd, Lfd=0) {
  #  Evaluates a monotone functional data observation, or the value of a linear
  #  differential operator LFD applied to the object,
  #  at argument values in an array EVALARGS.
  #  Functional data object LFD, if an integer, defines NDERIV, the
  #  order of derivative to be evaluated.
  #  Functional data object LFD, if a fd object, defines weight
  #  functions for computing the value of a linear differential operator
  #  applied to the functions that are evaluated.

  #  A monotone functional data object h  is in the form

  #           h(x) = [D^{-1} exp Wfd](x)

  #  where  D^{-1} means taking the indefinite integral.
  #  The interval over which the integration takes places is defined in
  #  the basisfd object in WFD.

  coef  <- getcoef(Wfd)
  coefd <- dim(coef)
  ndim  <- length(coefd)
  if (ndim > 2) stop("WFD is not a univariate function")
  if (ndim == 2) ncurve <- coefd[2] else ncurve <- 1

  if (is.numeric(Lfd)) {
    if (length(Lfd) == 1) {
      nderiv <- Lfd
      if (nderiv != as.integer(nderiv)) {
        stop("Order of derivative must be an integer")
      }
      if (nderiv < 0) {
        stop("Order of derivative must be 0 or positive")
      }
    } else {
      stop("Order of derivative must be a single number")
    }
  } else {
    stop("General linear differential operators not implemented yet.")
  }

  n       <- length(evalarg)

  hmat <- matrix(0, n, ncurve)

  for (icurve in 1:ncurve) {
	
  if (nderiv == 0) {
    hval <- monfn(evalarg, Wfd[icurve])
    hmat[,icurve] <- hval
  }

  if (nderiv == 1) {
    Dhval <- exp(eval.fd(evalarg, Wfd[icurve]))
    hmat[,icurve] <- Dhval
  }

  if (nderiv == 2) {
    basisfd <- getbasis(Wfd)
    Dwmat   <- getbasismatrix(evalarg, basisfd, 1)
    D2hval  <- (Dwmat %*% coef) * exp(eval.fd(evalarg, Wfd[icurve]))
    hmat[,icurve] <- D2hval
  }

  if (nderiv == 3) {
    basisfd <- getbasis(Wfd)
    Dwmat   <- eval.fd(evalarg, basisfd, 1)
    D2wmat  <- eval.fd(evalarg, basisfd, 2)
    D3hval  <- ((D2wmat %*% coef) + (Dwmat %*% coef)^2) *
                 exp(eval.fd(evalarg, Wfd[icurve]))
    hmat[,icurve] <- D3hval
  }

  if (nderiv > 3) stop ("Derivatives higher than 3 not implemented.")

  }
  
  return(hmat)

}
