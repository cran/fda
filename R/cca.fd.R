cca.fd <- function(fd, ncan = 2, lambda = 0.00025, Lfd = 2)
{
#  carry out a functional CCA with regularization.
#  Arguments:
#  FD        ... Functional data object.  It is assumed that there are
#                  two functions for each replication.
#  NCAN     ...  Number of pairs of canonical variates to be found
#  LAMBDA    ... Smoothing or regularization parameter.  The value 1 is
#                  arbitrary.  If lambda is a 2-vector then the first
#                  component will be applied to the "x" and the second to
#                  the "y" functions.
#  LFD  ... The order derivative to be penalized if an integer, or
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

  #  last modified 15 November 2001

  #  Check arguments

  if (!(inherits(fd, "fd"))) stop("Argument FD not a functional data object.")

  #  compute mean function and center if required

  ctrfd <- center.fd(fd)
  coef  <- getcoef(ctrfd)
  coefd <- dim(coef)
  ndim  <- length(coefd)

  if (ndim < 3 || coefd[3] == 1) stop(
     "CCA only possible with bivariate functions")
  if (coefd[3] > 2) warning("Only first two of multiple functions used")
  lambda <- c(lambda, lambda)
  if (lambda[1] <= 0 || lambda[2] <= 0) stop(
     "Smoothing parameters must be strictly positive")
  basisfd   <- getbasis(ctrfd)

  nbasis    <- basisfd$nbasis
  type      <- getbasistype(basisfd)
  nrep      <- coefd[2]
  if (nrep < 2) stop("CCA not possible without replications.")
#
#   Set up cross product matrices
#
  Jmat  <- getbasispenalty(basisfd, 0)
  Jx    <- t(Jmat %*% coef[,  , 1])
  Jy    <- t(Jmat %*% coef[,  , 2])
  PVxx  <- crossprod(Jx)/nrep
  PVyy  <- crossprod(Jy)/nrep
  if (inherits(Lfd, "fd") || (is.numeric(Lfd) && Lfd >= 0)) {
    Kmat  <- getbasispenalty(basisfd, Lfd)
    if (lambda[1] > 0) PVxx  <- PVxx + lambda[1] * Kmat
    if (lambda[2] > 0) PVyy  <- PVyy + lambda[2] * Kmat
  }
  Vxy   <- crossprod(Jx,Jy)/nrep
  #  do eigenanalysis
  result <- geigen(Vxy, PVxx, PVyy)
  #  set up canonical correlations and coefficients for weight functions
  canwtcoef      <- array(0,c(nbasis,ncan,2))
  canwtcoef[,,1] <- result$Lmat[,1:ncan]
  canwtcoef[,,2] <- result$Mmat[,1:ncan]
  corrs          <- result$values
  #   Normalize the coefficients for weight functions
  for (j in 1:ncan) for(k in 1:2) {
      temp <- canwtcoef[,j,k]
      temp <- temp/sqrt(sum(temp^2))
      canwtcoef[,j,k] <- temp
  }
  #  set up the canonical weight functions
  canwtfdnames      <- getnames(fd)
  canwtfdnames[[2]] <- paste('Can. Fn.',as.character(1:ncan))
  names(canwtfdnames)[2] <- 'Canonical functions'
  names(canwtfdnames)[3] <-
            paste('CCA wt. fns. for',names(canwtfdnames)[3])
  canwtfd           <- create.fd(canwtcoef, basisfd, canwtfdnames)
  #  set up canonical variable values
  canvarvalues      <- array(0, c(nrep, ncan, 2))
  canvarvalues[,,1] <- Jx %*% canwtcoef[,,1]
  canvarvalues[,,2] <- Jy %*% canwtcoef[,,2]
  #  set up return list
  cancorlist        <- list(canwtfd, corrs, canvarvalues)
  setOldClass("cca.fd")  
  oldClass(cancorlist) <- "cca.fd"
  names(cancorlist) <- c("weight functions", "correlations", "variates")

  return(cancorlist)
}
