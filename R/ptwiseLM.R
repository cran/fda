
#  ------------------------------------------------------------------------

ptwiseLM <- function (xfd, yfd, wbasis=xbasis, n=5*nbasisw,
                      estimate=rep(TRUE,ncoef),  constant=TRUE,
                      lambda=rep(0,ncoef), wfd0=rep(0,ncoef))
{
#  A function to compute the basis function expansion of the
#    estimate of the coefficient functions
#    for a pointwise linear model predicting the function(S)
#    in functional data object YFD from the p functions in
#    functional data object XFD.
#  The coefficient functions are expanded in terms of the
#    basis functions specified in wbasis.

#  Arguments:
#  XFD       ...  functional data object for P independent variable
#                   functions.  These must be univariate.
#  YFD       ...  functional data object for dependent variable
#                   functions
#  WBASIS    ...  basis object for regression functions.
#  N         ...  number of sampling points for numerical integration
#  ESTIMATE  ...  logical array of length P, if a value is T, the
#                 corresponding coefficient function is estimated, otherwise
#                 the target value is used.
#  CONSTANT  ...  If true, a constant function is added to the fit.
#  LAMBDA    ...  penalty parameter for penalizing the departure of the
#                 estimated weight functions from those defined in WFD0
#  WFD0      ...  A specification of a functional data object that is used for
#                 those weight functions not estimated, or as target functions
#                 toward which the estimated weight functions are smoothed. WFD0
#                 can either be a vector of NCOEF constants, or a functional
#                 data object with the same structure as WFN that is returned
#                 by this function.

#  Returns:
#  WFN       ...  estimated weight functional data object.  It has P + 1
#                 replicates if CONSTANT is T, and P otherwise

#  Last modified 6 Feb 2001

  if (!(inherits(xfd, "fd"))) stop(
       "Argument XFD not a functional data object.")
  if (!(inherits(yfd, "fd"))) stop(
       "Argument YFD not a functional data object.")

  coefx  <- as.matrix(getcoef(xfd))
  coefdx <- dim(coefx)
  ndimx  <- length(coefdx)
  ncurve <- coefdx[2]

  coefy  <- as.matrix(getcoef(yfd))
  coefdy <- dim(coefy)
  ndimy  <- length(coefdy)
  if (ndimy > 2) nvar <- coefdy[3] else nvar <- 1

  if (coefdy[2] != ncurve) stop(
      "Number of replications for XFD and YFD are not the same.")

  xbasis  <- getbasis(xfd)
  nbasisx <- xbasis$nbasis
 #nbasisx <- xbasis@nbasis
  ybasis  <- getbasis(yfd)
  nbasisy <- ybasis$nbasis
 #nbasisy <- ybasis@nbasis
  if (ndimx == 2) coefx <- array(coefx,c(nbasisx,ncurve,1))
  if (ndimx > 2) p  <- coefdx[3] else p <- 1

  typew   <- wbasis$type
  nbasisw <- wbasis$nbasis
  rangew  <- wbasis$rangeval
 #typew   <- wbasis@type
 #nbasisw <- wbasis@nbasis
 #rangew  <- wbasis@rangeval

  if (any(rangew != xbasis$rangeval)) stop(
    "Weight function range not equal to range in XFD")
 #if (any(rangew != xbasis@rangeval)) stop(
 #  "Weight function range not equal to range in XFD")

  if (typew == "bspline") {
    nbreaksw <- length(wbasis$params)
    norderw  <- nbasisw - nbreaksw
  }

  delta <- (rangew[2]-rangew[1])/(n-1)
  tfine <- seq(rangew[1],rangew[2],delta)

  yarray <- eval.fd(tfine, yfd)
  estimate <- as.logical(estimate)
  if (constant) {
    ncoef  <- length((1:(p+1))[estimate]) 
    xarray <- array(1,c(n,ncurve,p+1))
    xarray[,,2:(p+1)] <- eval.fd(tfine, xfd)
  } else {
    ncoef  <- length((1:p)[estimate])
    xarray <- eval.fd(tfine, xfd)
  }


  basismat <- getbasismatrix(tfine, wbasis)

  if (ncurve == 1) {
    DV <- -delta*yarray
    IV <- matrix(0,n,ncoef*nbasisw)
    mi <- 0
    for (i in 1:ncoef)
    {
      if(estimate[i]) {
        mi <- mi + 1
        index <- (1 + (mi-1)*nbasisw):(mi*nbasisw)
        IV[,index] <- delta*outer(xarray[,,i],rep(1,nbasisw))*basismat
      }
    }
    result <- lsfit(IV,DV,int=FALSE)
    dvec   <- result$coef
  } else {
    mi   <- 0
    mij  <- 0
    Swgt <- matrix(0,n,ncoef)
    Rwgt <- matrix(0,n,ncoef*(ncoef+1)/2)
    for (i in 1:ncoef)
    {
      if(estimate[i]) {
        mi <- mi + 1
        index <- (1 + (mi-1)*nbasisw):(mi*nbasisw)
        Swgt[,mi] <- apply(xarray[,,i]*yarray,1,mean)
        mj <- 0
        for (j in 1:i) {
          if(estimate[j]) {
	        mj <- mj + 1
            mij <- mij + 1
            Rwgt[,mij] <- apply(xarray[,,mi]*xarray[,,mj],1,mean)
          }
        }
      }
    }

    result <- SRsetup1(ncoef, nbasisw, Swgt, Rwgt, basismat)

    Cmat <- result[[2]]
    Dmat <- result[[1]]
    if (any(lambda > 0)) {
      if (!(inherits(wfd0, "fd")) && is.numeric(wfd0)) {
        if (length(wfd0) != ncoef) stop(
          "WFD0 is a vector of incorrect length")
        wcoef0 <- matrix(wfd0,nbasisw,ncoef)
        wfn0 <- create.fd(wcoef0, wbasis)
      } else {
        stop("WFN0 is neither a vector nor a FD object")
      }
      Hmat   <- getbasispenalty(wbasis)
      for (i in 1:ncoef) {
        index <- (1 + (i-1)*nbasisw):(i*nbasisw)
        if (lambda[i] > 0) {
          Cmat[index,index] <- Cmat[index,index] - lambda[i]*Hmat
          Dmat[index,1] <- Dmat[index,1] + lambda[i]*inprod(wbasis,wfn0[i])
        }
      }
    }
    dvec   <- solve(Cmat,Dmat)
  }

  dmat <- matrix(0,nbasisw,ncoef)
  mi  <- 0
  for (i in 1:ncoef)
  {
    if(estimate[i]) {
      mi <- mi + 1
      index <- (1 + (mi-1)*nbasisw):(mi*nbasisw)
      dmat[,i] <- dvec[index]
    }
  }

  wfnfd <- create.fd(dmat, wbasis)
 #wfnfd <- new("fd", coefs=dmat, basis=wbasis)
 
  return( wfnfd )
}

#  ------------------------------------------------------------------------

SRsetup1 <- function(nwgt, nbasis, Swgt, Rwgt, basismat)
{
  #  sets up coefficient matrices for basis expansion of weight functions

  Smat <- matrix(0, nwgt*nbasis, 1)
  Rmat <- matrix(0, nwgt*nbasis, nwgt*nbasis)
  n  <- nrow(Swgt)
  m1 <- ncol(basismat)
  m  <- 0
  one <- rep(1, nrow(basismat))
  for (i in 1:nwgt){
    indexi <- (1:nbasis) + (i-1)*nbasis
    temp     <- basismat * outer(Swgt[,i], rep(1,m1))
    temp[1,] <- temp[1,]/2
    temp[n,] <- temp[n,]/2
    Smat[indexi] <- crossprod(temp, one)
    for (j in 1:i) {
      m <- m + 1
      indexj <- (1:nbasis) + (j-1)*nbasis
      temp     <- basismat * outer(Rwgt[,m],rep(1,m1))
      temp[1,] <- temp[1,]/2
      temp[n,] <- temp[n,]/2
      Rmat[indexi,indexj] <- crossprod(temp, basismat)
      if (i != j) Rmat[indexj,indexi] <- Rmat[indexi,indexj]
    }
  }
  return (list(Smat, Rmat) )
}

