lambda2df <- function (argvals, basisfd, wtvec=rep(1,n), Lfd=NULL, lambda=0)
{
  #  Computes the the degrees of freedom associated with a regularized
  #    basis smooth by calculating the trace of the smoothing matrix.

  #  Arguments for this function:
  #
  #  ARGVALS  ... A set of argument values.
  #  BASISFD  ... A basis.fd object created by function create.basis.fd.
  #  WTVEC    ... A vector of N weights, set to one by default, that can
  #               be used to differentially weight observations in the
  #               smoothing phase
  #  LFD      ... The order of derivative or a linear differential
  #               operator to be penalized in the smoothing phase.
  #               By default Lfd is set in function GETBASISPENALTY
  #  LAMBDA   ... The smoothing parameter determining the weight to be
  #               placed on the size of the derivative in smoothing.  This
  #               is 0 by default.
  #  Returns:
  #  DF    ...  a degrees of freedom measure

  #  Last modified:  17 May 2000

  n        <- length(argvals)
  nbasis   <- basisfd$nbasis
  if (lambda == 0) {
    df <- nbasis
    return( df )
  }
  if (length(wtvec) != n) stop('WTVEC of wrong length')
  if (min(wtvec) <= 0)    stop('All values of WTVEC must be positive.')
  basismat <- getbasismatrix(argvals, basisfd)
  basisw   <- basismat*outer(wtvec,rep(1,nbasis))
  Bmat     <- crossprod(basisw,basismat)
  penmat   <- getbasispenalty(basisfd, Lfd)
  Bnorm    <- sqrt(sum(c(Bmat)^2))
  pennorm  <- sqrt(sum(c(penmat)^2))
  condno   <- pennorm/Bnorm
  if (lambda*condno > 1e12) {
    lambda <- 1e12/condno
    warning(paste("lambda reduced to",lambda,"to prevent overflow"))
  }
  Cmat     <- Bmat + lambda*penmat
  Cmat     <- (Cmat + t(Cmat))/2
  if (is.diag(Cmat)) {
      Cmatinv <- diag(1/diag(Cmat))
  } else {
      Lmat    <- chol(Cmat)
      Lmatinv <- solve(Lmat)
      Cmatinv <- crossprod(t(Lmatinv))
  }
  hatmat <- Cmatinv %*% Bmat
  df <- sum(diag(hatmat))
  return( df )
}
