eval.bifd <- function(sevalarg, tevalarg, bifd, sLfd = 0, tLfd = 0) {

  #  Evaluates a bi-functional data object BIFD at argument values in arrays
  #  SEVALARG and TEVALARG.  Differential operators SLFD and TLFD are
  #     applied to BIFD if present.

  #  Last modified 6 Feb 2001

  if (!is.vector(sevalarg)) stop(
     "First argument is not a vector.")
  if (!is.vector(tevalarg)) stop(
     "Second argument is not a vector.")

  ns   <- length(sevalarg)
  nt   <- length(tevalarg)

  if (!(inherits(bifd, "bifd"))) stop("Third argument is not a bifd object")

  sbasisfd <- bifd[[2]]
  snbasis  <- sbasisfd$nbasis
  rangeval <- sbasisfd$rangeval
  if (min(sevalarg) < rangeval[1] || max(sevalarg) > rangeval[2]) stop(
    "Values of the first argument are outside of permitted range.")

  tbasisfd <- bifd[[3]]
  tnbasis  <- tbasisfd$nbasis
  rangeval <- tbasisfd$rangeval
  if (min(tevalarg) < rangeval[1] || max(tevalarg) > rangeval[2]) stop(
    "Values of the second argument are outside of permitted range.")

  coef  <- bifd[[1]]
  coefd <- dim(coef)
  ndim  <- length(coefd)

  if (is.numeric(sLfd)) {
    if (length(sLfd) == 1) {
      snderiv <- sLfd
      if (snderiv != as.integer(snderiv)) {
        stop("Order of derivative must be an integer")
      }
      if (snderiv < 0) {
        stop("Order of derivative must be 0 or positive")
      }
    } else {
      stop("Order of derivative must be a single number")
    }
    sLfd <- NULL
    if (snderiv < 0) stop ("Order of derivative cannot be negative")
  } else if (inherits(sLfd, "fd")) {
    sderivcoef <- getcoef(sLfd)
    snderiv <- ncol(sderivcoef)
  } else {
    stop("Third argument must be an integer or a functional data object")
  }

  sbasismat <- getbasismatrix(sevalarg, sbasisfd, snderiv)
  if (snderiv > 0 && !is.null(sLfd)) {
    sLfdmat <- eval.fd(sevalarg, sLfd)
    onerow <- rep(1,snbasis)
    for (j in 1:snderiv) {
      if (any(abs(sLfdmat[,j])) > 1e-7) {
        sbasismat <- sbasismat + outer(sLfdmat[,j],onerow)*
                         getbasismatrix(sevalarg, sbasisfd, j-1)
      }
    }
  }

  if (is.numeric(tLfd)) {
    if (length(tLfd) == 1) {
      tnderiv <- tLfd
      if (tnderiv != as.integer(tnderiv)) {
        stop("Order of derivative must be an integer")
      }
      if (tnderiv < 0) {
        stop("Order of derivative must be 0 or positive")
      }
    } else {
      stop("Order of derivative must be a single number")
    }
    tLfd <- NULL
    if (tnderiv < 0) stop ("Order of derivative cannot be negative")
  } else if (inherits(tLfd, "fd")) {
    tderivcoef <- getcoef(tLfd)
    tnderiv <- ncol(tderivcoef)
  } else {
    stop("Third argument must be an integer or a functional data object")
  }

  tbasismat <- getbasismatrix(tevalarg, tbasisfd, tnderiv)
  if (tnderiv > 0 && !is.null(tLfd)) {
    tLfdmat <- eval.fd(tevalarg, tLfd)
    onerow <- rep(1,tnbasis)
    for (j in 1:tnderiv) {
      if (any(abs(tLfdmat[,j])) > 1e-7) {
        tbasismat <- tbasismat + outer(tLfdmat[,j],onerow)*
                         getbasismatrix(tevalarg, tbasisfd, j-1)
      }
    }
  }

  if (ndim == 2) {
    evalbifd <- sbasismat %*% coef %*% t(tbasismat)
  }
  if (ndim == 3) {
    nrep  <- coefd[3]
    evalbifd <- array(0,c(ns,nt,nrep))
    for (i in 1:nrep) {
      evalbifd[,,i] <- sbasismat %*% coef[,,i] %*% t(tbasismat)
    }
    dimnames(evalbifd) <- list(NULL,NULL,dimnames(coef)[[3]])
  }
  if (ndim > 3) {
    nrep  <- coefd[3]
    nvar  <- coefd[4]
    evalbifd <- array(0,c(ns,nt,nrep,nvar))
    for (i in 1:nrep) for (j in 1:nvar) {
      evalbifd[,,i,j] <-
        sbasismat %*% coef[,,i,j] %*% t(tbasismat)
    }
    dimnames(evalbifd) <-
        list(NULL,NULL,dimnames(coef)[[3]],dimnames(coef)[[4]])
  }
  return(evalbifd)
}
