cca.fd <- function(fdobj1, fdobj2=fdobj1, ncan = 2,
                   ccafdParobj1=fdPar(basisobj1, 2, 1e-10),
                   ccafdParobj2=ccafdParobj1,
                   centerfns=TRUE)
{
#  Carry out a functional CCA with regularization using two different
#    functional data samples.  These may have different bases, and
#    different penalty functionals.  It is assumed that both are
#    univariate.
#  Arguments:
#  FDOBJ1        ... Functional data object.
#  FDOBJ2        ... Functional data object.
#  NCAN          ... Number of pairs of canonical variates to be found
#  CCAFDPAROBJ1  ... A functional parameter object for the first set of
#                    canonical variables.
#  CCAFDPAROBJ2  ... A functional parameter object for the second set of
#                    canonical variables.
#  CENTERFNS     ... A logical variable indicating whether or not to
#                    center the functions before analysis.  Default is TRUE.
#
#  Returns:  An object of the CCA.FD class containing:
#  CCWTFD1       ... A functional data object for the first set of
#                    canonical variables
#  CCWTFD2       ... A functional data object for the second set of
#                    canonical variables
#  CANCORR       ... A vector of canonical correlations
#  CCAVAR1       ... A matrix of scores on the first canonical variable.
#  CCAVAR2       ... A matrix of scores on the second canonical variable.
#

# Last modified 2007.09.07 by Spencer Graves   
#  previously modified 27 March 2006

  if (!(inherits(fdobj1, "fd"))) stop(
		"Argument FDOBJ1 not a functional data object.")
  if (!(inherits(fdobj2, "fd"))) stop(
		"Argument FDOBJ2 not a functional data object.")
		
  if (centerfns) {
      fdobj1 <- center.fd(fdobj1)
      fdobj2 <- center.fd(fdobj2)
  }
  coef1  <- fdobj1$coefs
  coefd1 <- dim(coef1)
  ndim1  <- length(coefd1)
  nrep1  <- coefd1[2]

  coef2  <- fdobj2$coefs
  coefd2 <- dim(coef2)
  ndim2  <- length(coefd2)
  nrep2  <- coefd2[2]

  if (nrep1 != nrep2) stop("Numbers of replications are not equal.")
  if (nrep1 < 2) stop("CCA not possible without replications.")

  nrep <- nrep1

  if (ndim1 > 2 || ndim2 > 2) stop(
		"One or both functions are not univariate.")

  basisobj1  <- fdobj1$basis
  nbasis1    <- basisobj1$nbasis
  type1      <- basisobj1$type

  basisobj2  <- fdobj2$basis
  nbasis2    <- basisobj2$nbasis
  type2      <- basisobj2$type

  #   Set up cross product matrices

  Jmat1 <- eval.penalty(basisobj1, 0)
  Jmat2 <- eval.penalty(basisobj2, 0)
  Jx    <- t(Jmat1 %*% coef1)
  Jy    <- t(Jmat2 %*% coef2)
  PVxx  <- crossprod(Jx)/nrep
  PVyy  <- crossprod(Jy)/nrep
  Vxy   <- crossprod(Jx,Jy)/nrep

  #  add roughness penalty

  if (inherits(ccafdParobj1, "fd") || inherits(ccafdParobj1, "basisfd"))
      ccafdParobj1 <- fdPar(ccafdParobj1)
  if (inherits(ccafdParobj2, "fd") || inherits(ccafdParobj2, "basisfd"))
      ccafdParobj2 <- fdPar(ccafdParobj2)

  if (!inherits(ccafdParobj1, "fdPar")) stop(
		"ccafdParobj1 is not a fdPar object.")
  if (!inherits(ccafdParobj2, "fdPar")) stop(
		"ccafdParobj2 is not a fdPar object.")
		
  Lfdobj1 <- ccafdParobj1$Lfd
  Lfdobj2 <- ccafdParobj2$Lfd
  lambda1 <- ccafdParobj1$lambda
  lambda2 <- ccafdParobj2$lambda

  if (lambda1 > 0) {
    Kmat1 <- eval.penalty(basisobj1, Lfdobj1)
    PVxx  <- PVxx + lambda1 * Kmat1
  }
  if (lambda2 > 0) {
    Kmat2 <- eval.penalty(basisobj2, Lfdobj2)
    PVyy  <- PVyy + lambda2 * Kmat2
  }

  #  do eigenanalysis

  result <- geigen(Vxy, PVxx, PVyy)

  #  set up canonical correlations and coefficients for weight functions

  canwtcoef1 <- result$Lmat[,1:ncan]
  canwtcoef2 <- result$Mmat[,1:ncan]
  corrs      <- result$values

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

  canwtfdnames      <- fdobj1$fdnames
  canwtfdnames[[2]] <- "Canonical Variable"
  names(canwtfdnames)[2] <- "Canonical functions"
  names(canwtfdnames)[3] <-
            paste("CCA wt. fns. for",names(canwtfdnames)[3])
  canwtfd1 <- fd(canwtcoef1, basisobj1, canwtfdnames)
  canwtfd2 <- fd(canwtcoef2, basisobj2, canwtfdnames)
  canwtfd1$fdnames <- fdobj1$fdnames
  canwtfd2$fdnames <- fdobj2$fdnames

  #  set up canonical variable values

  canvarvalues1 <- Jx %*% canwtcoef1
  canvarvalues2 <- Jy %*% canwtcoef2

  #  set up return list

  ccafd        <- list(canwtfd1, canwtfd2, corrs,
                       canvarvalues1, canvarvalues2)
  names(ccafd) <- c("ccawtfd1", "ccawtfd2", "ccacorr",
                    "ccavar1",  "ccavar2")

#  setOldClass("cca.fd")
#  oldClass(ccafd) <- "cca.fd"
  class(ccafd) <- "cca.fd"

  return(ccafd)
}
