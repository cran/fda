cca.2fun.fd <- function(fd1, fd2=fd1, ncan = 2,
                         lambda1 = 0.00025, Lfd1 = 2,
                         lambda2 = lambda1, Lfd2 = Lfd1)
{
#  carry out a functional CCA with regularization using two different
#    functional data samples.  These may have different bases, and
#    different penalty functionals.  It is assumed that both are
#    univariate.
#  Arguments:
#  FD1       ... Functional data object.
#  FD2       ... Functional data object.
#  NCAN     ...  Number of pairs of canonical variates to be found
#  LAMBDA1  ...  Smoothing or regularization parameter for FD1.
#  LFD1     ... The order derivative of FD1 to be penalized if an integer, or
#           a linear differential operator if a functional data object.
#  LAMBDA2  ...  Smoothing or regularization parameter for FD2.
#  LFD2     ... The order derivative of FD2 to be penalized if an integer, or
#           a linear differential operator if a functional data object.
#
#  Returns:  A list with the entries
#  weight functions  ... A functional data object for the canonical
#                  variate weight functions
#  correlations    ... The corresponding set of canonical correlations
#  variates        ... An array of the values taken by the canonical variates
#                  (ie the scores on the variates.)  This is a 3-way array
#                       with first dimension corresponding to replicates,
#                       second to the different variates (dimension NCAN)
#                       and third (dimension 2) to the "x" and "y" scores.
#


 #  Last modified 6 Feb 2001

  if (!(inherits(fd1, "fd"))) stop("Argument FD1 not a functional data object.")
  if (!(inherits(fd2, "fd"))) stop("Argument FD2 not a functional data object.")

  ctrfd1 <- center.fd(fd1)
  coef1  <- getcoef(ctrfd1)
  coefd1 <- dim(coef1)
  ndim1  <- length(coefd1)
  nrep1      <- coefd1[2]

  ctrfd2 <- center.fd(fd2)
  coef2  <- getcoef(ctrfd2)
  coefd2 <- dim(coef2)
  ndim2  <- length(coefd2)
  nrep2      <- coefd2[2]

  if (nrep1 != nrep2) stop("Numbers of replications are not equal.")
  if (nrep1 < 2) stop("CCA not possible without replications.")
  nrep <- nrep1

  if (ndim1 > 2 || ndim2 > 2) stop("One or both functions are not univariate.")

  basisfd1   <- getbasis(fd1)
  nbasis1    <- basisfd1$nbasis
  type1      <- getbasistype(basisfd1)

  basisfd2   <- getbasis(fd2)
  nbasis2    <- basisfd2$nbasis
  type2      <- getbasistype(basisfd2)

#
#   Set up cross product matrices
#
  Jmat1 <- getbasispenalty(basisfd1, 0)
  Jmat2 <- getbasispenalty(basisfd2, 0)
  Jx    <- t(Jmat1 %*% coef1)
  Jy    <- t(Jmat2 %*% coef2)
  PVxx  <- crossprod(Jx)/nrep
  PVyy  <- crossprod(Jy)/nrep
  if (inherits(Lfd1, "fd") || (is.numeric(Lfd1) && Lfd1 >= 0)) {
    Kmat1 <- getbasispenalty(basisfd1, Lfd1)
    PVxx  <- PVxx + lambda1 * Kmat1
  }
  if (inherits(Lfd2, "fd") || (is.numeric(Lfd2) && Lfd2 >= 0)) {
    Kmat2 <- getbasispenalty(basisfd2, Lfd2)
    PVyy  <- PVyy + lambda2 * Kmat2
  }
  Vxy   <- crossprod(Jx,Jy)/nrep
  #  do eigenanalysis
  result <- geigen(Vxy, PVxx, PVyy)
  #  set up canonical correlations and coefficients for weight functions
  canwtcoef      <- array(0,c(nbasis,ncan,2))
  canwtcoef1 <- result$Lmat[,1:ncan]
  canwtcoef2 <- result$Mmat[,1:ncan]
  corrs          <- result$values
  #   Normalize the coefficients for weight functions
  for (j in 1:ncan) {
      temp <- canwtcoef1[,j]
      temp <- temp/sqrt(sum(temp^2))
      canwtcoef1[,j] <- temp
      temp <- canwtcoef2[,j]
      temp <- temp/sqrt(sum(temp^2))
      canwtcoef2[,j] <- temp
  }
  #  set up the canonical weight functions
  canwtfdnames      <- getnames(fd)
  canwtfdnames[[2]] <- paste("Can. Fn.",as.character(1:ncan))
  names(canwtfdnames)[2] <- "Canonical functions"
  names(canwtfdnames)[3] <-
            paste("CCA wt. fns. for",names(canwtfdnames)[3])
  canwtfd1 <- create.fd(canwtcoef1, basisfd1, canwtfdnames)
  canwtfd2 <- create.fd(canwtcoef2, basisfd2, canwtfdnames)
  #  set up canonical variable values
  canvarvalues      <- array(0, c(nrep, ncan, 2))
  canvarvalues[,,1] <- Jx %*% canwtcoef1
  canvarvalues[,,2] <- Jy %*% canwtcoef2
  #  set up return list
  cancorlist        <- list(canwtfd1, canwtfd2, corrs, canvarvalues)
  setOldClass("cca.2fun.fd")
  oldClass(cancorlist) <- "cca.2fun.fd"
  names(cancorlist) <- c("weight function 1", "weight function 2",
                         "correlations", "variates")
  return(cancorlist)
}
