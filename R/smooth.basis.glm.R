smooth.basis.glm <- function(argvals, y, fdParobj, wtvec=NULL, fdnames=NULL, 
                             covariates=NULL, family="binomial", dfscale=1) {
  #SMOOTH.GLM  Smooths discrete curve represented by basis function
  #  expansions fit by penalized least squares.
  #
  #  Required arguments for this function are
  #
  #  ARGVALS   A set of argument values, set by default to equally spaced
  #               on the unit interval (0,1).
  #  Y         If the family is not binomial, y is a matrix or an array 
  #               containing values of curves.
  #               If y is a matrix, rows must correspond to argument
  #               values and columns to replications, and it will be assumed
  #               that there is only one variable per observation.
  #               If Y is a three-dimensional array, the first dimension
  #               corresponds to argument values, the second to replications,
  #               and the third to variables within replications.
  #               If Y is a vector, only one replicate and variable are 
  #               assumed.
  #            If the family is binomial local sample sizes M.i,
  #               Y is a list vector of length 2, the first of which cantains
  #               the matrix or array as above containing observed frequencies,
  #               and the second of which contains the corresponding local
  #               sample sizes.
  #  FDPAROBJ  A functional parameter or fdPar object.  This object 
  #               contains the specifications for the functional data
  #               object to be estimated by smoothing the data.  See
  #               comment lines in function fdPar for details.
  #               This argument may also be either a FD object, or a 
  #               BASIS object.  If this argument is a basis object, the 
  #               smoothing parameter LAMBDA is set to 0.
  #
  #  Optional arguments are input in pairs  the first element of the pair
  #     is a string specifying the property that the argument value defines,
  #     and the second element is the value of the argument
  #
  #     Valid property/value pairs include
  #
  #     Property        Value
  #     ----------------------------------------------------------------
  #     weight          vector of the same length as the data vector to be
  #                     smoothed, containing nonnegative weights to be 
  #                     applied to the data values
  #     fdnames         A cell array of length 3 with names for
  #                       1. argument domain, such as "Time"
  #                       2. replications or cases
  #                       3. the function.
  #     covariates      A N by Q matrix Z of covariate values used to augment
  #                     the smoothing function, where N is the number of
  #                     data values to be smoothed and Q is the number of
  #                     covariates.  The process of augmenting a smoothing 
  #                     function in this way is often called "semi-parametric 
  #                     regression".  The default is the empty object NULL.
  #     dfscale         A scalar value multiplying the degrees of freedom
  #                     in the definition of the generalized 
  #                     cross-validated or GCV criterion for selecting the
  #                     bandwidth parameter LAMBDA.  It was recomm}ed by
  #                     Chong Gu that this be a number slightly larger than
  #                     1.0, such as 1.2, to prevent under-smoothing,
  #                     The default is 1.0.
  #     family          a character string containing one of
  #                       "normal"  
  #                       "binomial"
  #                       "poisson"
  #                       "gamma"
  #                       "inverse gaussian"
  #                     the value determines which of the link functions in
  #                     the generalized linear model (GLM) family is to be
  #                     used.  The default is "normal".
  #      control        a struct object controlling iterations with members
  #                       epsilon  convergence criterion (default 1e-8)
  #                       maxit    max. iterations       (default 25)
  #                       trace    output iteration info (0)
  #      start          a vector containing starting values for coefficients
  #                      
  #
  #  Returned objects are
  #
  #  FDOBJ    an object of class fd containing coefficients.
  #  DF       a degrees of freedom measure.
  #  GCV      a measure of lack of fit discounted for df.
  #                 If the function is univariate, GCV is a vector 
  #                 containing the stop  sum of squares for each 
  #                 function, and if the function is multivariate, 
  #                 GCV is a NVAR by NCURVES matrix.
  #  SSE      the stop sums of squares.  
  #                 SSE is a vector or matrix of the same size as 
  #                 GCV.
  #  PENMAT   the penalty matrix, if computed, otherwise NULL.
  #  Y2CMAP   the matrix mapping the data to the coefficients.
  #  ARGVALS  the input set of argument values.
  #  Y        the input array containing values of curves
  
  #  Last modified 15 May 2018 by Jim Ramsay
  
  n <- length(argvals)
  
  #  check ARGVALS
  
  if (!is.numeric(argvals)) stop("ARGVALS is not numeric.")
  argvals <- as.vector(argvals)
  if (length(argvals) < 2)  stop("ARGVALS does not contain at least two values.")
  
  #  check Y
  
  if (is.vector(y)) y <- as.matrix(y)
  c  
  #  check FDPAROBJ and get FDOBJ and LAMBDA
  
  if (!inherits(fdParobj, "fdPar")) {
    if (inherits(fdParobj, "fd") || inherits(fdParobj, "basisfd")) {
      fdParobj <- fdPar(fdParobj)
    } else
      stop(paste("'fdParobj' is not a functional parameter object,",
                 "not a functional data object, and",
                 "not a basis object."))
  }
  fdobj    <- fdParobj$fd
  lambda   <- fdParobj$lambda
  Lfdobj   <- fdParobj$Lfd
  
  #  check LAMBDA
  
  if (lambda < 0) { 
    lambda <- 0  
  }
  
  #  get BASIS and NBASIS
  
  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis - length(basisobj$dropind)
  
  #  check WTVEC
  
  # wtList <- wtcheck(n, wtvec)
  # wtvec  <- wtList[[1]]
  # onewt  <- wtList[[2]]
  # 
  # if (onewt) {
  #   wtvec <- matrix(1,n,1)
  # }
  
  if (is.null(wtvec)) wtvec <- matrix(1,n,1)
  
  #  check FDNAMES
  
  if (!is.null(fdnames) && !is.list(fdnames)) {
    stop("Optional argument FDNAMES is not a list object.")
  }
  
  if (is.list(fdnames) && length(fdnames) != 3) {
    stop("Optional argument FDNAMES is not of length 3.")
  }
  
  #  check COVARIATES
  
  q <- 0
  if (!is.null(covariates)) {
    if (!is.numeric(covariates)) {
      stop("Optional argument COVARIATES is not numeric.")
    }
    if (dim(covariates)[1] != n) {
      stop("Optional argument COVARIATES has incorrect number of rows.")
    }
    q <- dim(covariates)[2]
  }
  
  #  ------------------------------------------------------------------
  #                set up the linear equations for smoothing
  #  ------------------------------------------------------------------
  
  #  set up matrix of basis function values
  
  basismat <- eval.basis(argvals, basisobj)
  
  if (n >= nbasis || lambda > 0) {
    
    #  The following code is for the coefficients completely determined
    
    #  set up additional rows of the least squares problem for the
    #  penalty term.
    
    basismat0 <- basismat
    y0        <- y
    
    if (lambda > 0) {
      penmat  <- eval.penalty(basisobj, Lfdobj)
      lamRmat <- lambda*penmat
    } else {
      lamRmat <- NULL
    }
    
    #  augment BASISMAT0 and BASISMAT by the covariate matrix 
    #  if (it is supplied
    
    if (!is.null(covariates)) {
      basismat0 <- matrix(cbind(basismat0, covariates))
      basismat  <- matrix(cbind(basismat,  covariates))
      if (!is.null(lamRmat)) {
        lamRmat <- rbind(cbind(lamRmat,            matrix(0,nbasis,q)),
                         cbind(matrix(0,q,nbasis), matrix(0,q)       ))
      }
    }
    
    #  ------------------------------------------------------------------
    #               compute solution using Matlab function glmfit
    #  ------------------------------------------------------------------
    
    dimy   <- dim(y)
    nbasis <- dimy[1]
    ncurve <- dimy[2]
    ndim   <- length(dimy)
    if (ndim < 3) {
      coef  <- matrix(0,nbasis,ncurve)
      dev   <- matrix(0,ncurve,1)
      glmList <- glm.fda(basismat, y, family, lamRmat, wtvec)
      coef <- glmList[[1]]
      dev  <- glmList[[2]]
    } else {
      nvar  <- dimy[3]
      coef  <- array(0,c(nbasis,ncurve,nvar))
      dev   <- matrix(0,ncurve,nvar)
      for (ivar in 1:nvar) {
        yi <- as.matrix(y[,,ivar])
        glmList <- glm.fda(basismat, yi, family, lamRmat, wtvec)
        coefi  <- glmList[[1]]
        devi   <- glmList[[2]]
        statsi <- glmList[[3]]
        coef[,,ivar]  <- coefi
        dev[,ivar]    <- devi
        stats[[ivar]] <- statsi
      }
    }
    
    #  compute basismat*R^{-1}
    
    if (is.null(lamRmat)) {
      M <- crossprod(basismat)
    } else {
      M <- crossprod(basismat) + lamRmat
    }
    
    #  compute map from y to c
    
    y2cMap <- solve(M,t(basismat))
    
    #  compute degrees of freedom of smooth
    
    df <- sum(diag(basismat %*% y2cMap))
    
  } else {
    stop(paste("The number of basis functions exceeds the number of ", 
               "points to be smoothed."))    
  }
  
  #  ------------------------------------------------------------------
  #            compute SSE, yhat, GCV and other fit summaries
  #  ------------------------------------------------------------------
  
  # #  compute stop sum of squares
  # 
  # if (ndim < 3) {
  #   yhat <- basismat0 %*% coef
  #   SSE  <- sum((y0 - yhat)^2)
  # } else {
  #   SSE <- matrix(0,nvar,ncurve)
  #   for (ivar in 1:nvar) {
  #     coefi <- coef[,,ivar]
  #     yhati <- basismat %*% coefi
  #     yi    <- y[,,ivar]
  #     SSE[ivar,] <- sum((yi - yhati)^2)
  #   }
  # }
  # 
  # #  compute  GCV index
  # 
  # if (df < n) {
  #   gcv <- (SSE/n)/((n - dfscale*df)/n)^2
  # } else {
  #   gcv <- NA
  # }
  
  #  set up the functional data object
  
  if (ndim < 3) {
    fdobj <- fd(coef[1:nbasis,],  basisobj, fdnames)
  } else {
    fdobj <- fd(coef[1:nbasis,,], basisobj, fdnames)
  }
  
  #  set up the regression coefficient matrix beta
  
  if (q > 0) {
    ind <- (nbasis+1):(nbasis+q)
    if (ndim < 3) {
      beta <- coef[ind,]
    } else {
      beta <- coef[ind,,]
    }
  } else {
    beta <- NULL
  }
  
  return(list(fdobj=fdobj, beta=beta))
         
}



