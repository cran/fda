smooth.basis <- function(argvals, y, fdParobj,
                         wtvec=NULL, fdnames=NULL)
{
#  check argvals
  if (!is.numeric(argvals)) stop("'argvals' is not numeric.")
#
  dima <- dim(argvals)
  nda <- length(dima)
  if(nda<2)
    return(smooth.basis1(argvals, y=y, fdParobj=fdParobj,
                         wtvec=wtvec, fdnames=fdnames) )
#
  dimy <- dim(y)
  ndy <- length(dimy)
  if(ndy<nda)
    stop('argvals has ', nda, ' dimensions;  y has only ', ndy)
#
  if(all(dima==dimy[1:nda])){
    if(nda<3){
      sb2 <- smooth.basis2(argvals, y=y, fdParobj=fdParobj,
                           wtvec=wtvec, fdnames=fdnames)
      return(sb2)
    }
#    stop('smooth.basis not programmed for ',
#         'length(dim(argvals)) > 2')
    if(nda<4)
      return(smooth.basis3(argvals, y=y, fdParobj=fdParobj,
                           wtvec=wtvec, fdnames=fdnames) )
    cat('dim(argvals =', paste(dima, collapse=', '), '\n')
    stop('smooth.basis not programmed for ',
         'length(dim(argvalse)) > 3')
  } else {
    cat('dim(argvals) =', paste(dima, collapse=', '), '\n')
    cat('dim(y) = ', paste(dimy, collapse=''), '\n')
    stop('Dimensions of argvals do not match those of y')
  }
}

smooth.basis3 <- function(argvals, y, fdParobj,
                          wtvec=NULL, fdnames=NULL){
##
## 1.  dim(argvals) == dim(y)?
##
  dima <- dim(argvals)
  nda <- length(dima)
  if(nda<3)stop('length(dim(argvals)) must be 3;  is ', nda)
#
  dimy <- dim(y)
  ndy <- length(dimy)
  if(ndy<3)stop('length(dim(y)) must be 3;  is ', ndy)
#
  if(any(dima != dimy)){
    stop('dim(argvals) = ', paste(dima, collapse=', '),
         ' != dim(y) = ', paste(dimy, collapse=', '))
  }
##
## 2.  Call smooth.basis2 repeatedly
##
#  2.1.  argvals[, , 1]
  sb1 <- smooth.basis2(argvals[, , 1], y=y[, , 1], fdParobj=fdParobj,
                       wtvec=wtvec, fdname=fdnames)
#  2.2.  set up output object
  coef1 <- sb1$fd$coefs
  dimc1 <- dim(coef1)
  dimc <- c(dimc1[1], dimy[-1])
  coefs <- array(NA, dim=dimc)
  argNames <- dimnames(argvals)
  yNames <- dimnames(y)
  c1Names <- dimnames(coef1)
  cNames <- vector('list', 3)
  if(!is.null(c1Names[[1]])) cNames[[1]] <- c1Names[[1]]
  if(!is.null(yNames[[2]])) cNames[[2]] <- yNames[[2]]
  if(!is.null(yNames[[3]])) cNames[[3]] <- yNames[[3]]
#
  dimnames(coefs) <- cNames
#
  for(i in seq(2, length=dimy[3]-1)){
    sbi <- smooth.basis2(argvals[,,i], y=y[,,i], fdParobj=fdParobj,
                         wtvec=wtvec, fdnames=fdnames)
    coefs[,,i] <- sbi$fd$coefs
  }
  if(is.null(fdnames)){
    fdnames <- list(time=NULL, reps=NULL, values=NULL)
    if(!is.null(yNames[[1]])){
      fdnames[[1]] <- yNames[[1]]
    } else {
      if(!is.null(argNames[[1]]))
        fdnames[[1]] <- argNames[[1]]
    }
    if(!is.null(yNames[[2]])){
      fdnames[[2]] <- yNames[[2]]
    } else {
      if(!is.null(argNames[[2]]))
        fdnames[[2]] <- argNames[[2]]
    }
    if(!is.null(yNames[[3]])){
      fdnames[[3]] <- yNames[[3]]
    } else {
      if(!is.null(argNames[[3]]))
        fdnames[[3]] <- argNames[[3]]
    }
  }
##
## 3.  done
##
  sb <- sb1
  sb$fd$coefs <- coefs
  sb$fd$fdnames <- fdnames
  sb
}

smooth.basis2 <- function(argvals, y, fdParobj,
                          wtvec=NULL, fdnames=NULL){
##
## 1.  number of  dimensions of y = 2 or 3?
##
  dimy <- dim(y)
  ndy <- length(dimy)
  ynames <- dimnames(y)
  argNames <- dimnames(argvals)
##
## 2.  ndy = 2
##
  if(ndy<3){
#  2.1.  argvals[, 1]
    sb1 <- smooth.basis1(argvals[, 1], y=y[, 1], fdParobj=fdParobj,
                         wtvec=wtvec)
#  2.2.  set up output object
    dimc1 <- dim(sb1$fd$coefs)
    dimc <- c(dimc1[1], dimy[-1])
    coefs <- array(NA, dim=dimc)
    c1names <- dimnames(sb1$fd$coefs)
    cNames <- vector('list', 2)
    if(!is.null(c1names[[1]]))cNames[[1]] <- c1names[[1]]
    if(!is.null(ynames[[2]]))cNames[[2]] <- ynames[[2]]
    dimnames(coefs) <- cNames
    coefs[, 1] <- sb1$fd$coefs
#    dimnames(coefs) <- argNames
#
    for(i in seq(2, length=dimy[2]-1)){
      sbi <- smooth.basis1(argvals[, i], y=y[, i], fdParobj=fdParobj,
                           wtvec=wtvec)
      coefs[, i] <- sbi$fd$coefs
    }
    if(is.null(fdnames)){
      fdnames <- sb1$fdnames
      if(is.null(fdnames))
        fdnames <- list(time=NULL, reps=NULL, values='value')
      valueChk <- ((length(fdnames$values)==1)
                   && (fdnames$values=='value')
                   && (length(fdnames$reps)==1)
                   && (!is.null(ynames[[2]])) )
      if(valueChk)fdnames$values <- fdnames$reps
      if(!is.null(ynames[[2]]))
        fdnames[[2]] <- ynames[[2]]
    }
  } else{
##
## 3.  ndy = 3
##
#  3.1.  argvals[, 1]
    sb1 <- smooth.basis1(argvals[, 1], y=y[, 1, ], fdParobj=fdParobj,
                         wtvec=wtvec)
#  3.2.  set up output object
    coef1 <- sb1$fd$coefs
    dimc1 <- dim(coef1)
    dimc <- c(dimc1[1], dimy[-1])
    coefs <- array(NA, dim=dimc)
    yNames <- dimnames(y)
    c1Names <- dimnames(coef1)
    cNames <- vector('list', 3)
    if(!is.null(c1Names[[1]])) cNames[[1]] <- c1Names[[1]]
    if(!is.null(yNames[[2]])) cNames[[2]] <- yNames[[2]]
    if(is.null(c1Names[[2]])){
      if(!is.null(yNames[[3]]))
        cNames[[3]] <- yNames[[3]]
    } else cNames[[3]] <- c1Names[[2]]
# ** the middle dimension was dropped, in calling smooth.basis1
# ** so the second dimension of its output
# ** becomes the third of the desired output.
#
#    if(is.null(argNames)){
#      argNames <- vector('list', 3)
#      if(!is.null(yNames[[2]]))
#        argNames[[2]] <- yNames[[2]]
#    }
#    if(!is.null(yNames[[3]]))
#       argNames[[3]] <- yNames[[3]]
#
    dimnames(coefs) <- cNames
    coefs[, 1, ] <- coef1
#
    for(i in seq(2, length=dimy[2]-1)){
      sbi <- smooth.basis1(argvals[, i], y=y[, i, ], fdParobj=fdParobj,
                           wtvec=wtvec, fdnames=fdnames)
      coefs[, i, ] <- sbi$fd$coefs
    }
    if(is.null(fdnames)){
      fdnames <- sb1$fdnames
      if(is.null(fdnames)){
#        fdnames <- vector('list', 3)
        fdnames <- list(time=NULL, reps=NULL, values=NULL)
        if(!is.null(argNames[[1]])){
          fdnames[[1]] <- argNames[[1]]
        } else {
          fdnames[[1]] <- ynames[[1]]
        }
        if(!is.null(ynames[[2]]))fdnames[[2]] <- ynames[[2]]
        if(!is.null(ynames[[3]]))fdnames[[3]] <- ynames[[3]]
      }
    }
  }
##
## 4.  done
##
  sb <- sb1
  sb$fd$coefs <- coefs
  sb$fd$fdnames <- fdnames
  sb
}

smooth.basis1 <- function (argvals, y, fdParobj,
                          wtvec=NULL, fdnames=NULL)
{
# ARGVALS ... A set of argument values, set by default to equally spaced
#             on the unit interval (0,1).
# Y       ... an array containing values of curves
#             If the array is a matrix, rows must correspond to argument
#             values and columns to replications, and it will be assumed
#             that there is only one variable per observation.
#             If Y is a three-dimensional array, the first dimension
#             corresponds to argument values, the second to replications,
#             and the third to variables within replications.
#             If Y is a vector, only one replicate and variable are assumed.
# FDPAROBJ... A functional parameter or fdPar object.  This object
#             contains the specifications for the functional data
#             object to be estimated by smoothing the data.  See
#             comment lines in function fdPar for details.
#             This argument may also be either a FD object, or a
#             BASIS object.  In this case, the smoothing parameter
#             LAMBDA is set to 0.
# WTVEC   ... A vector of N weights, set to one by default, that can
#             be used to differentially weight observations.
# DFFACTOR... A multiplier of df in GCV, set to one by default
# FDNAMES ... A cell of length 3 with names for
#             1. argument domain, such as 'Time'
#             2. replications or cases
#             3. the function.
# Returns a list containing:
#   FDOBJ ...  an object of class fd containing coefficients.
#   DF    ...  a degrees of freedom measure.
#   GCV   ...  a measure of lack of fit discounted for df.
#              If the function is univariate, GCV is a vector
#              containing the error  sum of squares for each
#              function, and if the function is multivariate,
#              GCV is a NVAR by NCURVES matrix.
#   COEF  ...  the coefficient matrix for the basis function
#                expansion of the smoothing function
#   SSE   ...  the error sums of squares.
#              SSE is a vector or matrix of the same size as
#              GCV.
#   PENMAT...  the penalty matrix.
#   Y2CMAP...  the matrix mapping the data to the coefficients.

# last modified 2008.07.01 by Spencer Graves
#  previously modified:  2007.09.29 and 1 March 2007

#  ---------------------------------------------------------------------
#                      Check argments
#  ---------------------------------------------------------------------

  #  check ARGVALS

  if (!is.numeric(argvals)) stop("'argvals' is not numeric.")

  argvals <- as.vector(argvals)
  n      <- length(argvals)

  #  check Y

  if (is.vector(y)) y <- as.matrix(y)

  if(!inherits(y, "matrix") && !inherits(y, "array"))
    stop("'y' is not of class matrix or class array.")

  ydim <- dim(y);

  if (ydim[1] != n)
    stop("'y' is not the same length as 'argvals'.")
  if(is.null(wtvec))wtvec <- rep(1, n)

  #  check fdParobj

  if (!inherits(fdParobj, "fdPar")) {
    if (inherits(fdParobj, "fd") || inherits(fdParobj, "basisfd"))
        fdParobj <- fdPar(fdParobj)
    else
        stop(paste("'fdParobj' is not a functional parameter object,",
               "not a functional data object, and",
               "not a basis object."))
  }

  #  check WTVEC

  if (!is.vector(wtvec))  stop("'wtvec' is not a vector.")
  if (length(wtvec) != n) stop("'wtvec' of wrong length")
  if (min(wtvec) <= 0)    stop("All values of 'wtvec' must be positive.")

  #  extract information from fdParobj

  Lfdobj   <- fdParobj$Lfd
  nderiv   <- Lfdobj$nderiv
  fdobj    <- fdParobj$fd
  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  onebasis <- rep(1,nbasis)
  lambda   <- fdParobj$lambda

  #  check LAMBDA

  if (lambda < 0) {
    warning ("Value of 'lambda' was negative;  0 used instead.")
    lambda <- 0
  }

#  ---------------------------------------------------------------------
#                      Set up analysis
#  ---------------------------------------------------------------------

  #  set number of curves and number of variables

  ndim  <- length(ydim)

  tnames <- dimnames(y)[[1]]
  bnames <- fdParobj$fd$basis$names
#
  if(is.null(tnames))tnames <- 1:n
  if (ndim == 1) {
    nrep <- 1
    nvar <- 1
    coef <- rep(0,nbasis)
#   This should be unnecessary as when ndim<3, coef is overwritten below
#   names(coef) <- bnames
    y <- matrix(y,n,1, dimnames=list(tnames, NULL))
    ynames <- 'rep'
    vnames <- 'value'
  }
  if (ndim == 2)  {
    nrep <- ncol(y)
    nvar <- 1
    coef <- matrix(0,nbasis,nrep)
#   This should be unnecessary as when ndim<3, coef is overwritten below
#   dimnames(coef) <- list(bnames, ynames)
    ynames <- dimnames(y)[[2]]
    vnames <- 'value'
  }
  if (ndim == 3)  {
    nrep <- dim(y)[2]
    nvar <- dim(y)[3]
    coef <- array(0,c(nbasis,nrep,nvar))
    ynames <- dimnames(y)[[2]]
    vnames <- dimnames(y)[[3]]
    dimnames(coef) <- list(bnames, ynames, vnames)
  }

  #  set up matrix of basis function values

  basismat <- eval.basis(as.vector(argvals), basisobj)

#  ----------------------------------------------------------------
#                set up the linear equations for smoothing
#  ----------------------------------------------------------------

#  if (n >= nbasis || lambda > 0) {
  if (n > nbasis || lambda > 0) {

    #  The following code is for the coefficients completely determined

    basisw <- basismat*wtvec
    Bmat0  <- crossprod(basisw,basismat)

    #  set up right side of equations

    {
      if (ndim < 3) {
      	Dmat <- crossprod(basisw,y)
      } else {
        Dmat <- array(0, c(nbasis, nrep, nvar))
      	for (ivar in 1:nvar)
          Dmat[,,ivar] <- crossprod(basisw,y[,,ivar])
      }
    }
    if (lambda > 0) {
#  smoothing required, set up coefficient matrix for normal equations
      penmat    <- eval.penalty(basisobj, Lfdobj)
      Bnorm   <- sqrt(sum(c(Bmat0)^2))
      pennorm <- sqrt(sum(c(penmat)^2))
      condno  <- pennorm/Bnorm
      if (lambda*condno > 1e12) {
        lambda <- 1e12/condno
        warning(paste("lambda reduced to",lambda,
                      "to prevent overflow"))
      }
      Bmat <- Bmat0 + lambda*penmat
    } else {
      penmat <- matrix(0,nbasis,nbasis)
      Bmat   <- Bmat0
    }

    #  compute inverse of Bmat

    Bmat    <- (Bmat+t(Bmat))/2
    Lmat    <- try(chol(Bmat), silent=TRUE)
    {
      if(class(Lmat)=="try-error"){
        Beig <- eigen(Bmat, symmetric=TRUE)
        BgoodEig <- (Beig$values>0)
        Brank <- sum(BgoodEig)
        if(Brank<dim(Bmat)[1])
          warning("Matrix of basis function values has rank ",
                  Brank, " < dim(fdobj$basis)[2] = ",
                  length(BgoodEig), ";  ignoring null space")
        goodVec <- Beig$vectors[, BgoodEig]
        Bmatinv <- (goodVec %*% (Beig$values[BgoodEig] * t(goodVec)))
      }
      else {
        Lmatinv <- solve(Lmat)
        Bmatinv <- Lmatinv %*% t(Lmatinv)
      }
    }

#  ----------------------------------------------------------------
#       Compute the coefficients defining the smooth and
#            summary properties of the smooth
#  ----------------------------------------------------------------

    #  compute map from y to c

    y2cMap = Bmatinv %*% t(basisw)

    #  compute degrees of freedom of smooth

    BiB0 <- (Bmatinv %*% Bmat0)

    df. <- sum(diag(BiB0))

    #  solve normal equations for each observation

    if (ndim < 3) coef <- Bmatinv %*% Dmat
    else for (ivar in 1:nvar)
			coef[,,ivar] <- Bmatinv %*% Dmat[,,ivar]

  } else {
#  if((n <= nbasis) && (lambda<=0))
    if(n == nbasis) {
      y2cMap <- solve(basismat[, 1:n])
      df. <- n
      {
        if(ndim==1)
          coef[1:n] <- (y2cMap %*% y)
        else if (ndim==2)
          coef[1:n, ] <-(y2cMap %*% y)
        else
          for(ivar in 1:var)
            coef[1:n, , ivar] <- (y2cMap %*% y[,,ivar])
      }
      penmat <- matrix(0,nbasis,nbasis)
    } else {
      warning("The number of basis functions = ", nbasis, " exceeds ",
              n, " = the number of points to be smoothed.  ",
              "With no smoothing (lambda = 0), this will produce ",
              "a perfect fit to data that typically has wild ",
              "excursions between data points.")
#  ------------------------------------------------------------------
#  The following code is for the underdetermined coefficients:
#  the number of basis functions exceeds the number of argument values.
#  ------------------------------------------------------------------
#      qrlist <- qr(t(basismat))
#      Qmat   <- qr.Q(qrlist, complete=TRUE)
#      Rmat   <- t(qr.R(qrlist))
#      Q1mat  <- Qmat[,1:n]
#      Q2mat  <- as.matrix(Qmat[,(n+1):nbasis])
#      Hmat   <- getbasispenalty(basisobj)
#      Q2tHmat   <- crossprod(Q2mat,Hmat)
#      Q2tHQ2mat <- Q2tHmat %*% Q2mat
#      Q2tHQ1mat <- Q2tHmat %*% Q1mat
      w.5 <- sqrt(wtvec)
      b.5 <- (basismat * w.5)
      svdB <- svd(b.5)
      dGood <- with(svdB, d > sqrt(.Machine$double.eps)*d[1])
#      vd <- with(svdB, v %*% diag(1/d))
      v.d <- with(svdB, v / rep(d, each=nbasis))
      vdu <- tcrossprod(v.d, svdB$u)
      y2cMap <- (vdu * rep(w.5, each=nbasis))
#
      if (ndim < 3) {
#        z1mat <- solve(Rmat,y)
#        z2mat <- solve(Q2tHQ2mat, Q2tHQ1mat %*% z1mat)
#        coef <- Q1mat %*% z1mat + Q2mat %*% z2mat
        coef <- y2cMap %*% y
      } else {
        for (ivar in 1:nvar) {
#          z1mat <- solve(Rmat,y[,,ivar])
#          z2mat <- solve(Q2tHQ2mat, Q2tHQ1mat %*% z1mat)
#          coef[,,ivar] <- Q1mat %*% z1mat + Q2mat %*% z2mat
          coef[,,ivar] <- (y2cMap %*% y[, , ivar])
        }
      }
      df. <- n
      penmat <- matrix(0,nbasis,nbasis)
#      basisinv <- solve(basismat)
#      basisiw <- basisinv / rep(wt, each=n)
#      Bmatinv <- tcrossprod(basisiw, basisinv)
#      y2cMap = Bmatinv %*% t(basisw)
    }
  }
#  ----------------------------------------------------------------
#            compute SSE, yhat, GCV and other fit summaries
#  ----------------------------------------------------------------

  #  compute error sum of squares

  if (ndim < 3) {
    yhat <- basismat %*% coef
    SSE <- sum((y - yhat)^2)
    if(is.null(ynames))ynames <- dimnames(yhat)[[2]]
  } else {
    SSE <- 0
    yhat <- array(0,c(n, nrep, nvar))
    dimnames(yhat) <- list(dimnames(basismat)[[1]],
                           dimnames(coef)[[2]],
                           dimnames(coef)[[3]])
    for (ivar in 1:nvar) {
      yhat[,,ivar] <- basismat %*% coef[,,ivar]
      SSE <- SSE + sum((y[,,ivar] - yhat[,,ivar])^2)
    }
    if(is.null(ynames))ynames <- dimnames(yhat)[[2]]
    if(is.null(vnames))vnames <- dimnames(yhat)[[2]]
  }
  if(is.null(ynames))ynames <- paste('rep', 1:nrep, sep='')
  if(is.null(vnames))vnames <- paste('value', 1:nvar, sep='')

  #  compute  GCV index
  if (df. < n) {
    if (ndim < 3) {
      gcv <- rep(0,nrep)
      for (i in 1:nrep) {
        SSEi <- sum((y[,i] - yhat[,i])^2)
        gcv[i] <- (SSEi/n)/((n - df.)/n)^2
      }
      if(ndim>1)names(gcv) <- ynames
    } else {
      gcv <- matrix(0,nrep,nvar)
      for (ivar in 1:nvar) {
        for (i in 1:nrep) {
          SSEi <- sum((y[,i,ivar] - yhat[,i,ivar])^2)
          gcv[i,ivar] <- (SSEi/n)/((n - df.)/n)^2
        }
      }
      dimnames(gcv) <- list(ynames, vnames)
    }
  } else {
#    gcv <- NA
    gcv <- Inf
  }

  #  ------------------------------------------------------------------
  #          Set up the functional data objects for the smooths
  #  ------------------------------------------------------------------

  #  set up default fdnames

#  if (ndim == 1) defaultnames <- list("time", "reps", "values")
#  if (ndim == 2) defaultnames <- list("time",
#                                    paste("reps",as.character(1:nrep)),
#                                    "values")
#  if (ndim == 3) defaultnames <- list("time",
#                                    paste("reps",as.character(1:nrep)),
#                                    paste("values",as.character(1:nvar)) )
  if(is.null(fdnames)){
    fdnames <- list(time=tnames, reps=ynames, values=vnames)
#    if (ndim == 1) fdnames <- list("time", "reps", "values")
#    if (ndim == 2)fdnames <- list("time",
#          paste("reps",as.character(1:nrep)),"values")
#    if (ndim == 3) fdnames <- list("time",
#          paste("reps",as.character(1:nrep)),
#          paste("values",as.character(1:nvar)) )
#    names(fdnames) <- c("args", "reps", "funs")
  }
  fdobj <- fd(coef, basisobj, fdnames)
#
#  smoothlist <- list(fd=fdobj, df=df., gcv=gcv, coef=coef,
#                     SSE=SSE, penmat=penmat, y2cMap=y2cMap )
  smoothlist <- list(fd=fdobj, df=df., gcv=gcv,
                     SSE=SSE, penmat=penmat, y2cMap=y2cMap,
                     argvals=argvals, y=y)

  class(smoothlist) <- 'fdSmooth'
  return(smoothlist)
}
