linmod <- function(xfdobj, yfdobj, wtvec=rep(1,nrep),
                   xLfdobj=int2Lfd(2), yLfdobj=int2Lfd(2),
                   xlambda=0, ylambda=0)
{

  #  This function fits an unrestricted functional linear model for a 
  #  a functional dependent variable and a single function covariate.
  #  Smoothing is controlled by two parameters XLAMBDA and YLAMBDA,
  #  corresponding to the independent and dependent functional
  #  variables, respectively.

  #  Argument:
  #  XFDOBJ  ... If the independent variable is multivariate, a design matrix.
  #              If the independent variable is functional, a "fd" object.
  #  YFDOBJ  ... If the dependent variable is multivariate, a design matrix.
  #              If the dependent variable is functional, a "fd" object.
  #  WTVEC   ... a vector of weights
  #  XLFDOBJ ... A linear differential operator object for the independent
  #              variable.
  #  YLFDOBJ ... A linear differential operator object for the dependent
  #              variable.
  #  XLAMBDA ... a smoothing parameter for the independent variable
  #  YLAMBDA ... a smoothing parameter for the   dependent variable
  #  ZMATRNK ... actual rank of independent variable matrix for the
  #              functional DV/multivariate IV case

  #  Returns:  a list containing
  #  ALPHA  ... a vector of intercept values
  #  REGFD  ... a functional data object for the regression function

  #  Last modified:   4 December 2005

  #  check XLfdobj and YLfdobj

  xLfdobj <- int2Lfd(xLfdobj)
  yLfdobj <- int2Lfd(yLfdobj)

  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #              The multivariate IV and functional DV case
  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (inherits(yfdobj, "fd") && !inherits(xfdobj, "fd"))  {
      stop ("Use function fRegress for this model.")
    }


  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #             The functional IV and multivariate DV case
  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (inherits(xfdobj, "fd") && !(inherits(yfdobj, "fd"))) {

      stop ("Use function fRegress for this model.")
    }


  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #             The functional IV and functional DV case
  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if (inherits(xfdobj, "fd") && inherits(yfdobj, "fd")) {

    coefx   <- xfdobj$coefs
    coefy   <- yfdobj$coefs
    coefdx  <- dim(coefx)
    coefdy  <- dim(coefy)
    ndimx   <- length(coefdx)
    ndimy   <- length(coefdy)
    if (ndimx < 2) stop(
             "Linear modeling impossible with 1 replication")
    if (ndimx == 3) stop(
             "This version cannot accommodate multiple functional IVs")
    nrep <- coefdx[2]
    if (coefdy[2] != nrep) stop (
      "Numbers of observations in the first two arguments do not match.")

    rangewt <- range(wtvec)
    if (rangewt[1] < 0) stop("WTVEC must not contain negative values.")
    if (rangewt[1] == rangewt[2]) wtvar <- FALSE else wtvar <- TRUE

    xbasisobj <- xfdobj$basis
    nbasisx   <- xbasisobj$nbasis
    typex     <- xbasisobj$type
    rangevx   <- xbasisobj$rangeval

    ybasisobj <- yfdobj$basis
    nbasisy   <- ybasisobj$nbasis
    typey     <- ybasisobj$type
    rangevy   <- ybasisobj$rangeval

    if (length(wtvec) != nrep) stop("WTVEC of wrong length")
    if (min(wtvec) <= 0)    stop("All values of WTVEC must be positive.")
    if (xlambda < 0) warning (
              "Value of LAMBDA was negative, and 0 used instead.")

    jmatx   <- inprod(xbasisobj, xbasisobj)
    penmatx <- inprod(xbasisobj, xbasisobj, xLfdobj, xLfdobj)
    if (ndimx == 2) {
      zmatx   <- t(rbind(matrix(1,1,nrep),jmatx %*% coefx))
    } else {
      zmatx   <- t(rbind(matrix(1,1,nrep),jmatx %*% coefx[,,1]))
    }

    jmaty   <- inprod(ybasisobj, ybasisobj)
    penmaty <- inprod(ybasisobj, ybasisobj, yLfdobj, yLfdobj)

    alpha  <- rep(0,nbasisy)
    index  <- 2:(nbasisx+1)
    kmatx <- matrix(0,nbasisx+1,nbasisx+1)
    kmatx[index,index] <- penmatx
    if (wtvar) {
      zmatw <- sweep(zmatx,2,wtmat,"*")
      tempx <- solve(crossprod(zmatw,zmatx) + xlambda*kmatx)
      tempy <- solve(jmaty + ylambda*penmaty)
      gmat  <- tempx %*% crossprod(zmatw,t(coefy[,,1])) %*% jmaty %*% tempy
    } else {
      tempx <- solve(crossprod(zmatx) + xlambda*kmatx)
      tempy <- solve(jmaty + ylambda*penmaty)
      if (ndimy == 2) {
        gmat  <- tempx %*% crossprod(zmatx,t(coefy)) %*% jmaty %*% tempy
      } else {
        gmat  <- tempx %*% crossprod(zmatx,t(coefy[,,1])) %*% jmaty %*% tempy
      }
    }
    yhatcoef <- t(zmatx %*% gmat)
    bcoef <- matrix(0,nbasisx,nbasisy)
    bcoef <- gmat[index,]
    alpha <- as.matrix(gmat[1,])

    fdnames <- yfdobj$fdnames

    alphafdnames      <- fdnames
    alphafdnames[[2]] <- "Intercept"
    alphafd <- fd(alpha, ybasisobj, alphafdnames)

    regfd   <- bifd(bcoef, xbasisobj, ybasisobj)
    yhatfd  <- fd(yhatcoef,ybasisobj, fdnames)

    linmodlist <- list(alphafd=alphafd, regfd=regfd, yhatfd=yhatfd)

    return( linmodlist )
  }

}
