smooth.basis <- function (y, argvals, basisfd, wtvec=rep(1,n),
                          Lfd=NULL, lambda=0,
                          fdnames=list(NULL, dimnames(y)[2], NULL))
{

  #  Arguments for this function:
  #
  #  Y        ... an array containing values of curves
  #               If the array is a matrix, rows must correspond to argument
  #               values and columns to replications, and it will be assumed
  #               that there is only one variable per observation.
  #               If Y is a three-dimensional array, the first dimension
  #               corresponds to argument values, the second to replications,
  #               and the third to variables within replications.
  #               If Y is a vector, only one replicate and variable are assumed.
  #  ARGVALS  ... A vector of argument values.
  #  BASISFD  ... A basis.fd object created by function create.basis.fd.
  #  WTVEC    ... A vector of N weights, set to one by default, that can
  #               be used to differentially weight observations in the
  #               smoothing phase
  #  LFD      ... The order of derivative or a nonhomogeneous linear differential
  #               operator to be penalized in the smoothing phase.
  #               By default Lfd is set in function GETBASISPENALTY
  #  LAMBDA   ... The smoothing parameter determining the weight to be
  #               placed on the size of the derivative in smoothing.  This
  #               is 0 by default.
  #  FDNAMES  ... A list of length 3 with members containing
  #               1. a single name for the argument domain, such as "Time"
  #               2. a vector of names for the replications or cases
  #               3. a name for the function, or a vector of names if there
  #                  are multiple functions.
  #  Returns a list containing:
  #  FD    ...  an object of class fd containing coefficients
  #  DF    ...  a degrees of freedom measure
  #  GCV   ...  a measure of lack of fit discounted for df.

  #  Last modified:  5 May 2003

  n        <- length(argvals)
  nbasis   <- basisfd$nbasis
  onebasis <- rep(1,nbasis)

  #  check WTVEC

  if (!is.vector(wtvec)) stop("WTVEC is not a vector.")
  if (length(wtvec) != n) stop("WTVEC of wrong length")
  if (min(wtvec) <= 0)    stop("All values of WTVEC must be positive.")

  #  check LAMBDA

  if (lambda < 0) {
    warning ("Value of LAMBDA was negative, and 0 used instead.")
    lambda <- 0
  }

  #  check data array Y

  data  <- as.array(y)
  datad <- dim(data)
  ndim  <- length(datad)
  if (datad[1] != n) stop("First dimension of Y not compatible with ARGVALS.")

  #  set number of curves and number of variables

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

  #  set up matrix of basis function values

  basismat <- getbasismatrix(argvals, basisfd)

  #  set up the linear equations for smoothing

  if (n >= nbasis || lambda > 0) {

    #  The following code is for the coefficients completely determined

    basisw <- basismat*outer(wtvec,rep(1,nbasis))
    Bmat   <- crossprod(basisw,basismat)
    Bmat0  <- Bmat

    #  set up right side of equations

    if (ndim < 3) {
      	Dmat <- crossprod(basisw,data)
    } else {
	 	Dmat <- array(0, c(nbasis, ncurves, nvar))
      	for (ivar in 1:nvar) {
        	Dmat[,,ivar] <- crossprod(basisw,data[,,ivar])
      	}
    }
    if (lambda > 0) {
        #  smoothing required, set up coefficient matrix for normal equations
        penmat  <- getbasispenalty(basisfd, Lfd)
        Bnorm   <- sqrt(sum(c(Bmat)^2))
        pennorm <- sqrt(sum(c(penmat)^2))
        condno  <- pennorm/Bnorm
        if (lambda*condno > 1e12) {
          lambda <- 1e12/condno
          warning(paste("lambda reduced to",lambda,
                        "to prevent overflow"))
        }
        Bmat <- Bmat + lambda*penmat
    } else {
	  	penmat <- matrix(0,nbasis,nbasis)
       Bmat   <- Bmat
    }

    #  compute inverse of Bmat

    Bmat <- (Bmat+t(Bmat))/2
    if (is.diag(Bmat)) {
      	Bmatinv <- diag(1/diag(Bmat))
    } else {
      	Bmatinv <- solve(Bmat)
    }

    #  compute degrees of freedom of smooth

    df <- sum(diag(Bmatinv %*% Bmat0))

    #  solve normal equations for each observation

    if (ndim < 3) {
      	coef <- Bmatinv %*% Dmat
    } else {
      	for (ivar in 1:nvar) {
        	coef[,,ivar] <- Bmatinv %*% Dmat
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
  }

  #  compute error sum of squares

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

  #  compute  GCV index

  if (df < n) {
    	gcv <- (SSE/n)/(nvar*(n - df)/n)^2
  } else {
    	gcv <- NA
  }

  basislabels <- as.character(1:nbasis)
  if (ndim == 1) {
    coeflabs <- list(basislabels,NULL)
    names(coeflabs) <- c("basisfns", "reps")
  }
  if (ndim == 2) {
    coeflabs <- list(basislabels, fdnames[[2]])
    names(coeflabs) <- c("basisfns", "reps")
  }
  if (ndim == 3) {
    coeflabs <- list(basislabels, fdnames[[2]], fdnames[[3]])
    names(coeflabs) <- c("basisfns", "reps", "funs")
  }
  dimnames(coef) <- coeflabs

  fd <- create.fd(coef, basisfd, fdnames = fdnames)

  smoothlist <- list( fd, df, gcv, coef, SSE, penmat )
  names(smoothlist) <- c("fd", "df", "gcv", "coef", "SSE", "penmat")
  return( smoothlist )
}
