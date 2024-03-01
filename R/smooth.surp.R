smooth.surp <- function(argvals, y, Bmat0, WfdPar, wtvec=NULL, conv=1e-4,
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
  #  Y       ...  A matrix containing the values to be fit.
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
  
  #  Last modified 14 August 2023 by Jim Ramsay
  
  #  check ARGVALS, a vector of length n
  
  if (!is.numeric(argvals)) stop("ARGVALS is not numeric.")
  argvals <- as.vector(argvals)
  if (length(argvals) < 2) stop("ARGVALS does not contain at least two values.")
  n       <- length(argvals)
  onesobs <- matrix(1,n,1)
  
  #  Check y, an n by M-1 matrix of surprisal values.  
  #  It may not contain negative values.
  
  y    <- as.matrix(y)
  ydim <- dim(y)
  M    <- ydim[2]
  if (ydim[1] != n) 
      stop("The length of ARGVALS  and the number of rows of Y differ.")
  
  #  Check WfdPar and extract WBASIS, WNBASIS, Wlambda and WPENALTY.  
  #  Note that the coefficient matrix is not used.
  
  WfdPar   <- fdParcheck(WfdPar,M)
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
  
  surpList <- list(argvals=argvals, y=y, wtvec=wtvec, Kmat=Kmat,
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
  
  Flist <- list(f = fold, grad = gvec, norm = sqrt(mean(gvec^2)))
  
  Foldlist <- Flist
  
  #  evaluate the initial update vector for correcting the initial bmat
  
  pvec   <- -solve(hmat,gvec)
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
    #  set up new Bmat and evaluate function, gradient and hessian
    Bmatnew <- matrix(x,Wnbasis,M-1)
    func_result <- surp.fit(Bmatnew, surpList)
    f     <- func_result[[1]]
    gvec  <- func_result[[2]]
    hmat  <- func_result[[3]]
    SSE   <- func_result[[4]]
    DSSE  <- func_result[[5]]
    D2SSE <- func_result[[6]]
    #  set up list object for current fit
    Flist$f    <- f
    Flist$grad <- gvec
    Flist$norm <- sqrt(mean(gvec^2))
    xold <- x
    fold <- f
    #  display results at this iteration if dbglev > 0
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
      break
    }
    #  also terminate iterations if new fit is worse than old
    if (Flist$f >= Foldlist$f) break
    #  set up objects for new iteration
    #  evaluate the update vector using Newton Raphson direction
    pvec <- -solve(hmat,gvec)
    cosangle  <- -sum(gvec*pvec)/sqrt(sum(gvec^2)*sum(pvec^2))
    if (cosangle < 0) {
      if (dbglev > 1) print("cos(angle) negative")
      pvec <- -gvec
    }
    Foldlist <- Flist
  }
  
  #  end of iteration loop, output results
  
  Bmat <- matrix(xold, Wnbasis, M-1)
  Wfd  <- fda::fd(Bmat, Wbasis)
  surpResult <- surp.fit(Bmat, surpList)
  
  PENSSE   <- surpResult$PENSSE
  DPENSSE  <- surpResult$DPENSSE 
  D2PENSSE <- surpResult$D2PENSSE
  SSE      <- surpResult$SSE
  DSSE     <- surpResult$DSSE
  D2SSE    <- surpResult$D2SSE
  DvecSmatDvecB <- surpResult$DvecSmatDvecB

  surpFd <- list(Wfd=Wfd, Bmat=Bmat, f=f, gvec=gvec, hmat=hmat,
                 PENSSE=PENSSE, DPENSSE=DPENSSE, D2PENSSE=D2PENSSE,
                 SSE=SSE, DSSE=DSSE, D2SSE=D2SSE,
                 DvecSmatDvecB=DvecSmatDvecB)
  class(surpFd) <- 'surpfd'
  
  return(surpFd)
}

# ------------------------------------------------------------------

surp.fit <- function(x, surpList) {
  
  #  This function is called within smooth.surp() to
  #  evaluate fit at each iteration
  
  #  extract objects from surpList
  
  argvals <- surpList$argvals 
  y       <- surpList$y 
  wtvec   <- surpList$wtvec 
  Kmat    <- surpList$Kmat
  Zmat    <- surpList$Zmat 
  Phimat  <- surpList$Phimat 
  
  # set up dimensions and Bmat
  
  n       <- length(argvals)
  M       <- surpList$M
  K       <- dim(Phimat)[2]
  Bmat    <- matrix(x, K, M-1)
  
  #  compute fit, gradient and hessian
  
  logM     <- log(M)
  onewrd   <- all(wtvec == 1)
  Xmat     <- Phimat %*% Bmat %*% t(Zmat)
  expXmat  <- M^Xmat
  sumexpXmat <- as.matrix(apply(expXmat,1,sum))
  Pmat     <- expXmat/(sumexpXmat %*% matrix(1,1,M))
  Smat     <- -Xmat + (log(sumexpXmat) %*% matrix(1,1,M))/logM
  Rmat     <- y - Smat
  vecBmat  <- matrix(Bmat,K*(M-1),1,byrow=TRUE)
  vecRmat  <- matrix(Rmat,n*M,    1,byrow=TRUE)
  vecKmat  <- kronecker(diag(rep(1,M-1)),Kmat)
  fitscale <- 1
  #  compute raw fit and its penalized version
  if (!onewrd) {
    vecwtmat <- diag(rep(wtvec,M))
    SSE <- t(vecRmat) %*% vecwtmat %*% vecRmat
  } else {
    SSE <- crossprod(vecRmat)
  }
  PENSSE   <- SSE/fitscale + t(vecBmat) %*% vecKmat %*% vecBmat
  #  compute raw gradient and its penalized version
  DvecXmatDvecB <- kronecker(Zmat,Phimat)
  DvecSmatDvecX <- matrix(0,n*M,n*M)
  m2 <- 0
  for (m in 1:M) {
    m1 <- m2 + 1
    m2 <- m2 + n
    m4 <- 0
    for (l in 1:M) {
      m3 <- m4 + 1
      m4 <- m4 + n
      diagPl <- diag(Pmat[,l])
      DvecSmatDvecX[m1:m2,m3:m4] <- diagPl
    }
  }
  DvecSmatDvecX <- DvecSmatDvecX - diag(rep(1,n*M))
  DvecSmatDvecB <- DvecSmatDvecX %*% DvecXmatDvecB
  if (!onewrd) {
    DSSE <- -2*t(DvecSmatDvecB) %*% vecwtmat %*% vecRmat
  } else {
    DSSE <- -2*t(DvecSmatDvecB) %*% vecRmat
  }
  DPENSSE  <- DSSE/fitscale + 2*vecKmat %*% vecBmat
  #  compute raw hessian and its penalized version
  if (!onewrd) {
    D2SSE <- 2*t(DvecSmatDvecB) %*% vecwtmat %*% DvecSmatDvecB
  } else {
    D2SSE <- 2*crossprod(DvecSmatDvecB)
  }
  D2PENSSE <- D2SSE/fitscale + 2*vecKmat
  
  #  return list object containing raw and penalized fit data
  
  return(list(
    PENSSE   = PENSSE, 
    DPENSSE  = DPENSSE, 
    D2PENSSE = D2PENSSE,
    SSE      = SSE,
    DSSE     = DSSE,
    D2SSE    = D2SSE,
    DvecSmatDvecB = DvecSmatDvecB)
  )
  
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
