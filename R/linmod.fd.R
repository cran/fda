linmod.fd <- function(xfd, yfd, wtvec=rep(1,nrep),
                      xLfd=2, yLfd=2, xlambda=0, ylambda=0, zmatrnk=p)
{

  #  This function one of three types of linear model,
  #  each model consisting of a constant term
  #  plus either a conventional independent variable matrix or a single
  #  functional independent variable.  The modeling problem may
  #  be functional either in terms of the independent variable or in
  #  terms of the dependent variable, but at least one variable must
  #  be functional.
  #  Smoothing is controlled by two parameters XLAMBDA and YLAMBDA,
  #  corresponding to the independent and dependent functional
  #  variables, respectively.

  #  Argument:
  #  XFD     ... If the independent variable is multivariate, a design matrix.
  #              If the independent variable is functional, a 'fd' object.
  #  YFD     ... If the dependent variable is multivariate, a design matrix.
  #              If the dependent variable is functional, a 'fd' object.
  #  WTVEC   ... a vector of weights
  #  XLFD    ... For the independent variable, the order derivative to be
  #              penalized if an integer, or
  #              a linear differential operator if a functional data object.
  #  YLFD    ... For the dependent variable, the order derivative to be
  #              penalized if an integer, or
  #              a linear differential operator if a functional data object.
  #  XLAMBDA ... a smoothing parameter for the independent variable
  #  YLAMBDA ... a smoothing parameter for the   dependent variable
  #  ZMATRNK ... actual rank of independent variable matrix for the
  #              functional DV/multivariate IV case

  #  Returns:  a list containing
  #  ALPHA  ... a vector of intercept values
  #  REGFD  ... a functional data object for the regression function

  #  Last modified:   4 July 2001

  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #              The multivariate IV and functional DV case
  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (inherits(yfd, "fd") && !(inherits(xfd, "fd")))  {
    if (is.matrix(xfd)) {
      zmat <- as.matrix(xfd)
      fd   <- yfd
    }  else {
      stop ("First argument not a vector or a matrix.")
    }

    coef    <- getcoef(fd)
    coefd   <- dim(coef)
    ndim    <- length(coefd)
    if (ndim < 2) stop(
             "Linear modeling impossible with 1 replication")

    nrep <- nrow(zmat)
    if(ndim == 3) nvar <- coefd[3] else nvar <- 1
    basisfd <- getbasis(fd)
    nbasis  <- basisfd$nbasis

    rangewt <- range(wtvec)
    if (rangewt[1] < 0) stop("WTVEC must not contain negative values.")
    if (rangewt[1] == rangewt[2]) wtvar <- FALSE else wtvar <- TRUE

    if (length(wtvec) != nrep) stop("WTVEC of wrong length")
    if (min(wtvec) <= 0)    stop("All values of WTVEC must be positive.")
    if (ylambda < 0) warning (
       "Value of LAMBDA was negative, and 0 used instead.")
    if (ylambda > 0 && yLfd < 0) stop(
       "Order of derivative must be nonnegative.")

    if (dim(coef)[2] != nrep) stop(
       "Number of rows of ZMAT must equal number of replications")
    p    <- ncol(zmat)

    if (nvar > 1)
    {
      bcoef <- array(0,c(nbasis,p,nvar))
    } else {
      bcoef <- matrix(0,nbasis,p)
    }
    if (zmatrnk > p) stop(
       "Specified rank of ZMAT must not be greater than no. columns.")
    if (zmatrnk < p) {
      rootw   <- sqrt(wtvec)
      zmatsvd <- svd(sweep(zmat,2,rootw,"*"))
      zmatd   <- zmatsvd$d
      if (zmatd[zmatrnk] <= 0) stop("ZMAT is not of specified column rank")
      index  <- 1:zmatrnk
      zginvt <- sweep(zmatsvd$u[,index],1,zmatd[index],"/") %*%
                t(zmatsvd$v[,index])
      if (nvar == 1)
      {
        bcoef <- sweep(coef,1,rootw,"*") %*% zginvt
      } else {
        for (j in 1:nvar) {
          bcoef[,,j] <- sweep(coef[,,j],1,rootw,"*") %*% zginvt
        }
      }
    } else {
      if (nvar == 1)
      {
        bcoef <- t(lsfit(zmat,t(coef),wtvec,int=FALSE)$coef)
      } else {
        for (j in 1:nvar)
        {
          bcoef[,,j] <- t(lsfit(zmat,t(coef[,,j]),wtvec,int=FALSE)$coef)
        }
      }
    }
    yhatcoef <- bcoef %*% t(zmat)

    if (nvar > 1)
    {
      dimnames(bcoef) <- list(NULL,dimnames(zmat)[[2]],
                                 dimnames(yfd[[1]])[[3]])
    } else {
      dimnames(bcoef) <- list(NULL,dimnames(zmat)[[2]])
    }

    fdnames <- getnames(fd)

    regfdnames      <- fdnames
    regfdnames[[2]] <- paste('Reg. Coef.',1:p)
    regfdnames[[3]] <- 'Reg. Coef.'
    names(regfdnames)[2] <- 'Reg. Coefficients'
    regfd  <- create.fd (bcoef, basisfd, regfdnames)

    yhatfd <- create.fd (yhatcoef, basisfd, fdnames)

    linmodlist <- list(0, regfd, yhatfd)
    names(linmodlist) <- c('alpha', 'regfd', 'yhatfd')
    return( linmodlist )
  }

  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #             The functional IV and multivariate DV case
  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (inherits(xfd, "fd") && !(inherits(yfd, "fd"))) {

    if (is.numeric(yfd)) {
      ymat <- as.matrix(yfd)
      fd   <- xfd
    }  else {
      stop (paste("Second argument not a functional data object",
                  " and is not numeric",
                  "when first argument is functional."))
    }

    coef   <- getcoef(fd)
    coefd  <- dim(coef)
    ndim   <- length(coefd)
    if (ndim < 2) stop(
             "Linear modeling impossible with 1 replication")
    if (length(dim(coef)) == 3) stop(
             "This version cannot accommodate multiple functional IVs")
    nrep <- coefd[2]
    nvar <- ncol(ymat)

    rangewt <- range(wtvec)
    if (rangewt[1] < 0) stop("WTVEC must not contain negative values.")
    if (rangewt[1] == rangewt[2]) wtvar <- FALSE else wtvar <- TRUE

    basisfd <- getbasis(fd)
    nbasis  <- basisfd$nbasis
    type    <- getbasistype(basisfd)
    if (nrow(ymat) != nrep) stop(
      "Number of rows of YMAT must equal number of replications")
    one  <- rep(1,nrep)

    if (length(wtvec) != nrep) stop("WTVEC of wrong length")
    if (xlambda < 0) warning (
              "Value of XLAMBDA was negative, and 0 used instead.")

    jmat <- inprod(basisfd, basisfd)
    zmat <- t(rbind(one,jmat %*% coef))

    bcoef <- matrix(0,nbasis,nvar)
    alpha <- rep(0,nvar)
    index <- 2:(nbasis+1)

    if (xlambda <= 0) {
      #  no smoothing required, do ordinary least squares
      if (ncol(zmat) > nrow(zmat)) stop(paste(
         "For XLAMBDA = 0, no. of basis functions exceeds",
         "number of replications. No fit possible.")
)
      if (wtvar)
      {
        temp <- sweep(zmat,2,wtvec,"*")
        Cmat <- crossprod(temp,zmat)
        Dmat <- crossprod(temp,ymat)
      } else {
        Cmat <- crossprod(zmat)
        Dmat <- crossprod(zmat,ymat)
      }
    } else {
      #  smoothing required
      kmat <- matrix(0,nbasis+1,nbasis+1)
      kmat[index,index] <- inprod(basisfd, basisfd, xLfd, xLfd)
      if (wtvar)
      {
        temp <- sweep(zmat,2,wtvec,"*")
        Cmat <- crossprod(temp,zmat) + xlambda*kmat
        Dmat <- crossprod(temp,ymat)
      } else {
        Cmat <- crossprod(zmat) + xlambda*kmat
        Dmat <- crossprod(zmat,ymat)
      }
    }
    temp  <- symsolve( Cmat, Dmat )
    yhat  <- zmat %*% temp
    bcoef <- as.matrix(temp[index,])

    alpha <- temp[1,]

    fdnames <- getnames(xfd)

    regfdnames      <- fdnames
    regfdnames[[2]] <- paste('Reg. Coef.',1:nvar)
    regfdnames[[3]] <- 'Reg. Coef.'
    names(regfdnames)[2] <- 'Reg. Coefficients'
    regfd  <- create.fd (bcoef, basisfd, regfdnames)

    linmodlist <- list(alpha, regfd, yhat)
    names(linmodlist) <- c('alpha', 'regfd', 'yhat')
    return( linmodlist )
  }

  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #             The functional IV and functional DV case
  #  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if (inherits(xfd, "fd") && inherits(yfd, "fd")) {

    coefx   <- getcoef(xfd)
    coefy   <- getcoef(yfd)
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

    basisxfd <- getbasis(xfd)
    nbasisx  <- basisxfd$nbasis
    typex    <- basisxfd$type
    rangevx  <- basisxfd$rangeval

    basisyfd <- getbasis(yfd)
    nbasisy  <- basisyfd$nbasis
    typey    <- basisyfd$type
    rangevy  <- basisyfd$rangeval

    if (length(wtvec) != nrep) stop("WTVEC of wrong length")
    if (min(wtvec) <= 0)    stop("All values of WTVEC must be positive.")
    if (xlambda < 0) warning (
              "Value of LAMBDA was negative, and 0 used instead.")

    jmatx   <- inprod(basisxfd, basisxfd)
    penmatx <- inprod(basisxfd, basisxfd, xLfd, xLfd)
    if (ndimx == 2) {
      zmatx   <- t(rbind(matrix(1,1,nrep),jmatx %*% coefx))
    } else {
      zmatx   <- t(rbind(matrix(1,1,nrep),jmatx %*% coefx[,,1]))
    }

    jmaty   <- inprod(basisyfd, basisyfd)
    penmaty <- inprod(basisyfd, basisyfd, yLfd, yLfd)

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

    fdnames <- getnames(yfd)

    alphafdnames      <- fdnames
    alphafdnames[[2]] <- 'Intercept'
    names(alphafdnames)[2] <- 'Intercept'
    alphafd <- create.fd(alpha, basisyfd, alphafdnames)

    regfd   <- create.bifd(bcoef, basisxfd, basisyfd)
    yhatfd  <- create.fd(yhatcoef,basisyfd, fdnames)

    linmodlist <- list(alphafd, regfd, yhatfd)
    names(linmodlist) <- c('alphafd', 'regfd', 'yhatfd')
    return( linmodlist )
  }

}
