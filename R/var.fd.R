var.fd <- function(fdx, fdy = fdx)
{
  #  compute the variance and covariance functions for functional observations

  #  Last modified 6 Feb 2001

  if (!(inherits(fdx, "fd"))) stop("Argument FDX not a functional data object.")
  if (!(inherits(fdy, "fd"))) stop("Argument FDY not a functional data object.")

  coefx   <- getcoef(fdx)
  coefdx  <- dim(coefx)
  basisx  <- getbasis(fdx)
  nbasisx <- basisx$nbasis
  coefy   <- getcoef(fdy)
  coefdy  <- dim(coefy)
  basisy  <- getbasis(fdy)
  nbasisy <- basisy$nbasis
  if (coefdx[2] != coefdy[2]) stop(
    "Number of replications are not equal.")
  if (length(coefdx) == 2) {
    if(length(coefdy) == 2) {
      coefvar <- var(t(coefx),t(coefy))
      coefnames <- list(dimnames(coefx)[[1]], dimnames(coefy)[[1]])
      varbifd <- create.bifd(coefvar, basisx, basisy, coefnames)
    } else {
      stop("Both arguments must be univariate.")
    }
  } else {
    nvar <- coefdx[3]
    npair <- nvar*(nvar+1)/2
    coefnames <- list(dimnames(coefx[[1]]), dimnames(coefx[[1]]),
                      "Covariance", rep(" ",npair))
    coefvar <- array(0,c(nbasisx,nbasisx,1,npair))
    m <- 0
    for (i in 1:nvar) for (j in 1:i) {
      m <- m + 1
      coefvar[,,1,m] <- var(t(coefx[,,i]),t(coefx[,,j]))
      dimnames(coefx)[[4]][m]  <- paste(dimnames(coefx)[[3]][i],"vs",
                                        dimnames(coefx)[[3]][j])
    }
    varbifd <- create.bifd(coefvar, basisx, basisx, coefnames)
  }
  return(varbifd)
}
