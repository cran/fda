smooth.surp <- function(argvals, Wbin, Bmat0, WfdPar, wtvec=NULL, conv=1e-4,
                       iterlim=50, dbglev=0) {
  #  Smooths the relationship of Y to ARGVALS using weights in WTVEC by fitting 
  #     surprisal functions to a set of surprisal transforms of choice 
  #     probabilities, where the surprisal transformation of each probability is 
  #                      W(p_m) = -log_M (p_m), m=1, ..., M,
  #     where  W  is a function defined over the same range as ARGVALS.
  #  The fitting criterion is penalized mean squared error:
  #    PENSSE(Wlambda) <- \sum w_i[y_i - f(x_i)]^2 +
  #                     \Wlambda * \int [L W(x)]^2 dx
  #  where L is a linear differential operator defined in argument Lfd,
  #  and w_i is a positive weight applied to the observation.
  #  The function W(x) is expanded by the basis in functional data object
  #    Wfd.
  #
  #  This version uses Numerical Recipes function lnsrch
  #
  #  Arguments:
  #  ARGVALS ...  Argument value array of length NBIN, the number of 
  #               surprisal values for each curve.  It is assumed that
  #               that these argument values are common to all observed
  #               curves.  
  #  WBIN    ...  A matrix containingg the values to be fit.
  #               This will be an NBIN by M matrix, where NBIN is the 
  #               number of bins containing choice probabilities and M is
  #               the number of options in a specific question or rating
  #               scale.
  #  BMAT0   ...  An initial K by M-1 matrix defining the surprisal curves
  #               via spline functions.  K is the number of basis functions
  #               in the spline expansions, and M is the number of choices
  #               for a particular question in a test or rating scale.
  #  WFDPAR  ...  A functional parameter or fdPar object.  This object
  #               contains the specifications for the functional data
  #               object to be estimated by smoothing the data.  
  #               Note:  WFDPAR is only a container for its 
  #               functional basis object WBASIS, the penalty matrix WPEN, 
  #               and the smoothing parameter Wlambda.  A coefficient
  #               matrix in WFDPAR defined by using a function data object
  #               is discarded, and overwritten by argument BMAT0.
  #  WTVEC   ...  a vector of weights, a vector of N one's by default.
  #  CONV    ...  convergence criterion, 0.0001 by default
  #  ITERLIM ...  maximum number of iterations, 50 by default.
  #  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
  #               no output, if 1, output at each iteration, if higher,
  #               output at each line search iteration. 1 by default.
  #               enabling this option.
  
  #  Returns are:
  #  WFD     ...  Functional data object for W.
  #               Its coefficient matrix an N by NCURVE (by NVAR) matrix
  #               (or array), depending on whether the functional
  #               observations are univariate or multivariate.
  #  FLIST ... A list object or a vector of list objects, one for
  #            each curve (and each variable if functions are multivariate).
  #            Each list object has slots:
  #                 f    ... The sum of squared errors
  #                 grad ... The gradient
  #                 norm ... The norm of the gradient
  #  When multiple curves and variables are analyzed, the lists containing
  #  FLIST objects are indexed linear with curves varying inside
  #  variables.
  
  #  Last modified 25 June 2021 by Jim Ramsay
  
  #  check ARGVALS, a vector of length n
  
  if (!is.numeric(argvals)) stop("ARGVALS is not numeric.")
  argvals <- as.vector(argvals)
  if (length(argvals) < 2) stop("ARGVALS does not contain at least two values.")
  n       <- length(argvals)
  onesobs <- matrix(1,n,1)
  
  #  Check Wbin, an n by M-1 matrix of surprisal values.  
  #  It may not contain negative values.
  
  Wbin <- as.matrix(Wbin)
  Wbindim <- dim(Wbin)
  M <- Wbindim[2]
  if (Wbindim[1] != n) 
      stop("The length of ARGVALS  and the number of rows of WBIN differ.")
  # if (min(Wbin) < 0) stop("WBIN contains negative values.")
  
  #  Check WfdPar and extract WBASIS, WNBASIS, Wlambda and WPENALTY.  
  #  Note that the coefficient matrix is not used.
  
  WfdPar   <- fdParcheck(WfdPar)
  Wbasis   <- WfdPar$fd$basis
  Wnbasis  <- Wbasis$nbasis
  Wlambda  <- WfdPar$lambda
  Wpenalty <- eval.penalty(Wbasis, WfdPar$Lfd)
  
  #  Check BMAT0, the WNBASIS by M-1 coefficient matrix
  
  if (is.null(Bmat0)) stop("BMAT0 is  NULL.")
  
  Bmatdim <- dim(Bmat0)
  if (Bmatdim[1] != Wnbasis) 
    stop("The first dimension of BMAT0 is not equal to WNBASIS.")
  if (Bmatdim[2] != M-1) 
    stop("The second dimension of BMAT0 is not equal to M - 1.")
  
  #  convert Bmat0 to a vector and NPAR to its length
  
  cvec <- as.vector(Bmat0)
  npar <- length(cvec)
  
  #  Set up the transformation from dimension M-1 to M
  #  where M-vectors sum to zero
  
  M <- dim(Bmat0)[2] + 1
  if (M == 2) {
    root2 <- sqrt(2)
    Zmat <- matrix(1/c(root2,-root2),2,1)
  } else {
    Zmat <- zerobasis(M)
  }
  
  #  Set up the matrix of basis function values 
  
  Phimat <- fda::eval.basis(argvals, Wbasis)
  
  #  check WTVEC
  
  if (is.null(wtvec)) wtvec<-rep(1,n)
  wtvec <- fda::wtcheck(n, wtvec)$wtvec
  
  #  initialize matrix Kmat defining penalty term
  
  if (Wlambda > 0) 
  {
    Kmat <- Wlambda*Wpenalty
  } else {
    Kmat <- matrix(0,Wnbasis,Wnbasis)
  }
  
  #  Set up list object for data required by PENSSEfun
  
  surpList <- list(argvals=argvals, Wbin=Wbin, wtvec=wtvec, Kmat=Kmat,
                   Zmat=Zmat, Phimat=Phimat, M=M)
  #  --------------------------------------------------------------------
  #              loop through variables and curves
  #  --------------------------------------------------------------------
  
  #  evaluate log likelihood
  #    and its derivatives with respect to these coefficients
  
  xold <- matrix(Bmat0, Wnbasis*(M-1),1)
  result    <- surp.fit(xold, surpList)
  PENSSE    <- result[[1]]
  DPENSSE   <- result[[2]]
  D2PENSSE  <- result[[3]]
  
  #  compute initial badness of fit measures
  
  fold <- PENSSE
  gvec <- DPENSSE
  hmat <- D2PENSSE
  
  # print("gvec0:")
  # print(round(t(gvec0),4))
  # print("hmat0[1:10,1:10]:")
  # print(round(hmat0[1:10,1:10],4))
  
  Flist <- list(f = fold, grad = gvec, norm = sqrt(mean(gvec^2)))
  
  Foldlist <- Flist
  
  #  evaluate the initial update vector for correcting the initial bmat
  
  pvec   <- -solve(hmat,gvec)
  # print("pvec:")
  # print(round(t(pvec),4))
  cosangle <- -sum(gvec*pvec)/sqrt(sum(gvec^2)*sum(pvec^2))
  
  #  initialize iteration status arrays
  
  iternum <- 0
  status <- c(iternum, Foldlist$f, Foldlist$norm)
  if (dbglev >= 1) {
    cat("\n")
    cat("\nIter.   PENSSE   Grad Length")
    cat("\n")
    cat(iternum)
    cat("        ")
    cat(round(status[2],4))
    cat("      ")
    cat(round(status[3],4))
  }
  
  #  -------  Begin iterations  -----------
  
  STEPMAX <- 10
  
  iternum <- 0
  for (iter in 1:iterlim) {
    iternum <- iternum + 1
    #  take optimal stepsize
    lnsrch_result <- 
      lnsrch(xold, fold, gvec, pvec, surp.fit, surpList, STEPMAX)
    x     <- lnsrch_result$x
    check <- lnsrch_result$check
    if (check) stop("lnsrch failure")
    Bmatnew <- matrix(x,Wnbasis,M-1)
    func_result <- surp.fit(Bmatnew, surpList)
    f     <- func_result[[1]]
    gvec  <- func_result[[2]]
    hmat  <- func_result[[3]]
    Flist$f    <- f
    Flist$grad <- gvec
    Flist$norm <- sqrt(mean(gvec^2))
    xold <- x
    fold <- f
    status <- c(iternum, Flist$f, Flist$norm)
    if (dbglev > 0) {
      cat("\n")
      cat(iternum)
      cat("        ")
      cat(round(status[2],4))
      cat("      ")
      cat(round(status[3],4))
    }
    #  test for convergence
    if (abs(Flist$f - Foldlist$f) < conv) {
      # cat("\n")
      break
    }
    if (Flist$f >= Foldlist$f) break
    #  evaluate the update vector
    pvec <- -solve(hmat,gvec)
    cosangle  <- -sum(gvec*pvec)/sqrt(sum(gvec^2)*sum(pvec^2))
    if (cosangle < 0) {
      if (dbglev > 1) print("cos(angle) negative")
      pvec <- -gvec
    }
    Foldlist <- Flist
    #  end of iteration loop
  }
  
  Bmat <- matrix(xold, Wnbasis, M-1)
  Wfd  <- fda::fd(Bmat, Wbasis)
  
  result <- list(Wfd=Wfd, Bmat=Bmat)
  
  return(result)
}

# ------------------------------------------------------------------
ycheck <- function(y, n) {
  
  #  check Y
  
  if (is.vector(y)) y <- as.matrix(y)
  
  if (!inherits(y, "matrix") && !inherits(y, "array"))
    stop("Y is not of class matrix or class array.")
  
  ydim <- dim(y)
  
  if (ydim[1] != n) stop("Y is not the same length as ARGVALS.")
  
  #  set number of curves and number of variables
  
  ndim  <- length(ydim)
  if (ndim == 2) {
    ncurve <- ydim[2]
    nvar   <- 1
  }
  if (ndim == 3) {
    ncurve <- ydim[2]
    nvar   <- ydim[3]
  }
  if (ndim > 3) stop("Second argument must not have more than 3 dimensions")
  
  
  return(list(y=y, ncurve=ncurve, nvar=nvar, ndim=ndim))
  
}

# ------------------------------------------------------------------
fdParcheck <- function (fdPar) {
  if (!inherits(fdPar, "fdPar")) {
    if (inherits(fdPar, "fd") || inherits(fdPar, "basisfd")) {
      fdPar <- fdPar(fdPar)
    } else
      stop(paste("'fdPar' is not a functional parameter object,",
                 "not a functional data object, and",
                 "not a basis object."))
  }
  
  return(fdPar)
  
}

# ------------------------------------------------------------------

zerobasis <- function(k) {
# ZEROBASIS constructes a K by K-1 matrix that maps an unrestricted matrix B with K - 1 rows by 
#  the linear transformation ZEROBASIS %*% B = C into the subspace of matrices with K rows having #  column sums equal to zero.  
#  The matrix has orthonormal columns, so that crossprod(ZEROBASIS) is the identity matrix
#  of order K - 1.

  tk <- 0:(k-1) + 0.5
  fbasis     <- fda::create.fourier.basis(k,k)
  fbasmat    <- fda::eval.basis(tk, fbasis)
  fbasmat    <- fbasmat[,2:k]
  fbasnorm   <- sqrt(apply(fbasmat^2,2,sum))
  zerobasmat <- fbasmat/outer(rep(1,k),fbasnorm)
  return(zerobasmat)
}
