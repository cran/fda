eval.bifd <- function(sevalarg, tevalarg, bifd, sLfdobj = 0, tLfdobj = 0) {

  #  Evaluates a bi-functional data object BIFD at argument values in arrays
  #  SEVALARG and TEVALARG.  Differential operators SLFD and TLFD are
  #     applied to BIFD if present.

  #  Last modified 26 October 2005

  if (!is.vector(sevalarg)) stop(
     "First argument is not a vector.")
  if (!is.vector(tevalarg)) stop(
     "Second argument is not a vector.")

  ns   <- length(sevalarg)
  nt   <- length(tevalarg)

  if (!(inherits(bifd, "bifd"))) stop("Third argument is not a bifd object")

  sbasisobj <- bifd$sbasis
  snbasis   <- sbasisobj$nbasis
  rangeval  <- sbasisobj$rangeval
  if (min(sevalarg) < rangeval[1] || max(sevalarg) > rangeval[2]) stop(
    "Values of the first argument are outside of permitted range.")

  tbasisobj <- bifd$tbasis
  tnbasis   <- tbasisobj$nbasis
  rangeval  <- tbasisobj$rangeval
  if (min(tevalarg) < rangeval[1] || max(tevalarg) > rangeval[2]) stop(
    "Values of the second argument are outside of permitted range.")

  coef  <- bifd$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)

  sLfdobj <- int2Lfd(sLfdobj)
  tLfdobj <- int2Lfd(tLfdobj)

  snderiv <- sLfdobj$nderiv
  tnderiv <- tLfdobj$nderiv

  sbasismat <- getbasismatrix(sevalarg, sbasisobj, snderiv)
  if (snderiv > 0 && !is.null(sLfd)) {
    sLfdmat <- eval.fd(sevalarg, sLfd)
    onerow  <- rep(1,snbasis)
    for (j in 1:snderiv) {
      if (any(abs(sLfdmat[,j])) > 1e-7) {
        sbasismat <- sbasismat + outer(sLfdmat[,j],onerow)*
                         getbasismatrix(sevalarg, sbasisobj, j-1)
      }
    }
  }

  tbasismat <- getbasismatrix(tevalarg, tbasisobj, tnderiv)
  if (tnderiv > 0 && !is.null(tLfd)) {
    tLfdmat <- eval.fd(tevalarg, tLfd)
    onerow <- rep(1,tnbasis)
    for (j in 1:tnderiv) {
      if (any(abs(tLfdmat[,j])) > 1e-7) {
        tbasismat <- tbasismat + outer(tLfdmat[,j],onerow)*
                         getbasismatrix(tevalarg, tbasisobj, j-1)
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
