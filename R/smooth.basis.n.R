smooth.basis.n <- function (y, argvals, basisfd, wt=rep(1,n), penspecs=NULL)
{

  #  this version inputs the penalty specs as a vector of lists

  #  Arguments for this function:
  #
  #  Y        ... An array of discrete curve values to smooth.
  #  ARGVALS  ... A vector of argument  values
  #  BASISFD  ... A basis.fd object created by function create.basis.fd.
  #  WT       ... A vector of N weights, set to one by default, that can
  #               be used to differentially weight observations in the
  #               smoothing phase
  #  PENSPECS ... A list.  The members are lists, each containing:
  #    LFD      ... An order of derivative or a
  #                 linear differential operator to be penalized in the
  #                 smoothing phase.
  #                 2  by default
  #    LAMBDA   ... A smoothing parameter determining the weight
  #                 to be placed on the size of the derivative in smoothing.
  #                 This is 0 by default.
  #    PENRNG   ... A vector containing 1 or 2 elements
  #                 specifying the range over which the
  #                 integration is to be carried out.  If missing, the
  #                 full range is used.  If of length 1, the functional
  #                 is evaluated at that point rather than integrated.
  #    WTFD     ... A functional data object specifying a weight function
  #                 The default is NULL

  #  Returns an object of class fd containing coefficients
  #    for the expansion and BASISFD

  #  Last modified 16 May 1999

  n <- length(argvals)

  basismat <- getbasismatrix(argvals, basisfd)
  nbasis   <- dim(basismat)[2]
  onebasis <- rep(1,nbasis)

  if (length(wt) != n) stop("WT of wrong length")
  if (min(wt) <= 0)    stop("All values of WT must be positive.")

  data  <- as.array(y)
  datad <- dim(data)
  ndim  <- length(datad)
  if (ndim == 1) {
    nrep <- 1
    nvar <- 1
    coef <- rep(0,nbasis)
    data <- matrix(data,datad,1)
  }
  if (ndim == 2)  {
    nrep <- ncol(data)
    nvar <- 1
    coef <- matrix(0,nbasis,nrep)
  }
  if (ndim == 3)  {
    nrep <- dim(data)[2]
    nvar <- dim(data)[3]
    coef <- array(0,c(nbasis,nrep,nvar))
  }

  if (n >= nbasis) {

    #  The following code is for the coefficients completely determined
    basisw   <- basismat*outer(wt,rep(1,nbasis))
    Bmat     <- crossprod(basisw,basismat)

    if (is.null(penspecs)) {
      nspecs <- 0
    } else {
      nspecs <- length(penspecs)
    }
    Cmat <- Bmat
    if (nspecs >  0) {
      #  smoothing required, set up coefficient matrix for normal equations
      for (ispec in 1:nspecs) {
        penspec <- penspecs[[ispec]]
        penmat  <- getbasispenalty.n(basisfd, penspec)
        lambda  <- penspec$lambda
        Cmat    <- Cmat + lambda*penmat
      }
    }
    Cmat <- (Cmat + t(Cmat))/2

    #  compute inverse of Cmat

    if (is.diag(Cmat)) {
        Cmatinv <- 1/Cmat
    } else {
      Lmat    <- chol(Cmat)
      Lmatinv <- solve(Lmat)
      Cmatinv <- crossprod(t(Lmatinv))
    }

    #  compute degrees of freedom of smooth

    df <- sum(diag(Cmatinv %*% Bmat))

    #  solve normal equations for each observation

    if (ndim < 3) {
      Dmat <- crossprod(basisw,data)
      coef <- Cmatinv %*% Dmat
    } else {
      for (ivar in 1:nvar) {
        Dmat <- crossprod(basisw,data[,,ivar])
        coef[,,ivar] <- Cmatinv %*% Dmat
      }
    }

  } else {
    #  The following code is for the underdetermined coefficients:
    #     the number of basis functions exceeds the number of argument values.

    qrlist <- qr(t(basismat))
    Qmat   <- qr.Q(qrlist, complete=TRUE)
    Rmat   <- t(qr.R(qrlist))
    Q1mat  <- Qmat[,1:n]
    Q2mat  <- as.matrix(Qmat[,(n+1):nbasis])
    Hmat   <- getbasispenalty(basisfd)
    Q2tHmat   <- crossprod(Q2mat,Hmat)
    Q2tHQ2mat <- Q2tHmat %*% Q2mat
    Q2tHQ1mat <- Q2tHmat %*% Q1mat
    if (ndim < 3) {
      z1mat <- solve(Rmat,data)
      z2mat <- solve(Q2tHQ2mat, Q2tHQ1mat %*% z1mat)
      coef <- Q1mat %*% z1mat + Q2mat %*% z2mat
    } else {
      for (ivar in 1:nvar) {
        z1mat <- solve(Rmat,data[,,ivar])
        z2mat <- solve(Q2tHQ2mat, Q2tHQ1mat %*% z1mat)
        coef[,,ivar] <- Q1mat %*% z1mat + Q2mat %*% z2mat
      }
    }
    df <- n
  }

  #  compute  GCV index

  if (df < n) {
    if (ndim < 3) {
      datahat <- basismat %*% coef
      SSE <- sum((data - datahat)^2)
    } else {
      SSE <- 0
      for (ivar in 1:nvar) {
        datahat <- basismat %*% coef[,,ivar]
        SSE <- SSE + sum((data[,,ivar] - datahat)^2)
      }
    }
    gcv <- (SSE/n)/(nvar*(n - df)/n)^2
  } else {
    gcv <- NA
  }

  prefdnames <- dimnames(prefd[[1]])
  if (ndim == 1) dimnames(coef) <- NULL
  if (ndim == 2) dimnames(coef) <- list(NULL, prefdnames[[2]])
  if (ndim == 3) dimnames(coef) <- list(NULL, prefdnames[[2]], prefdnames[[3]])

  fd <- create.fd(coef, basisfd, df = df, gcv = gcv)

  return(fd)
}
