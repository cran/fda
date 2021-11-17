smooth.morph2 <- function(x, y, ylim, WfdPar, wt=matrix(1,nobs,1), 
                          conv=1e-4, iterlim=20, dbglev=0) {
  #SMOOTH_MORPH smooths the relationship of Y to ARGVALS 
  #  by fitting a monotone fn.  f(argvals) <- b_0 + b_1 D^{-1} exp W(t)
  #     where  W  is a function defined over the same range as ARGVALS,
  #  W + ln b_1 <- log Df and w <- D W <- D^2f/Df.
  #  b_0 and b_1 are chosen so that values of f
  #  are within the interval [ylim[1],ylim[2]].
  #  The fitting criterion is penalized mean squared stop:
  #    PENSSE(lambda) <- \sum [y_i - f(t_i)]^2 +
  #                     \lambda * \int [L W]^2 
  #  W(argvals) is expanded by the basis in functional data object Wfdobj.
  #
  #  Arguments are 
  #  X         argument value array
  #  Y         data array containing the the values to be fit
  #  YLIM      Ordinate value limits. The values of the estimated function
  #               range between these limits.  The abscissa range is 
  #               defined in WfdPar.
  #  WFDPAR   A functional parameter or fdPar object.  This object 
  #               contains the specifications for the functional data
  #               object to be estimated by smoothing the data.  See
  #               comment lines in function fdPar for details.
  #               The functional data object WFD in FDPAROBJ is used
  #               to initialize the optimization process.
  #               It's coefficient array has a single column, and these 
  #               are the starting values for the iterative minimization 
  #               of mean squared stop.
  #               This argument may also be either a FD object, or a 
  #               BASIS object.  In this case, the smoothing parameter 
  #               LAMBDA is set to 0.
  #  WT        a vector of weights
  #  CONV      convergence criterion, 0.0001 by default
  #  ITERLIM   maximum number of iterations, 20 by default
  #  DBGLEV    Controls the level of output on each iteration.  if (0,
  #               no output, if (1, output at each iteration, if (higher, 
  #               output at each line search iteration. 1 by default.
  #
  #  Returns are:
  #  WFD       Functional data object for W.  It's coefficient vector
  #               contains the optimized coefficients.
  
  #  last modified 31 October 2021
  
  nobs <- length(x)        #  number of observations
  
  # check consistency of x and y
  
  if (length(y) != nobs) {
    stop('Arguments X and Y are not of same length')
  }
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  #  check WfdPar
  
  if (!is.fdPar(WfdPar)) {
    if (is.fd(WfdPar)) WfdPar <- fdPar(WfdPar)
    else stop(paste('WFDPAR is not a functional parameter object, ', 
                    'and not a functional data object.'))
  }
  
  #  extract information from WfdPar
  
  Wfdobj   <- WfdPar$fd
  Wbasis   <- Wfdobj$basis       #  basis for Wfdobj
  Wnbasis  <- Wbasis$nbasis      #  no. basis functions
  Wtype    <- Wbasis$type
  xlim     <- Wbasis$rangeval
  WLfdobj  <- int2Lfd(WfdPar$Lfd)
  Wlambda  <- WfdPar$lambda
  if (any(wt < 0)) { 
    stop('One or more weights are negative.') 
  }
  
  cvec <- Wfdobj$coefs   #  initial coefficients
  
  #  first coefficient is not active for bspline or fourer bases
  
  if (Wtype == 'bspline' || Wtype == 'fourier') {
    active <- 2:Wnbasis
  } else {
    active <- 1:Wnbasis
  }
  inactive <- rep(TRUE,Wnbasis)
  inactive[active] = FALSE
  
  #  check range of argvals
  
  if (x[1] < xlim[1] || x[nobs] > xlim[2]) {
    stop('Values in ARGVALS are out of bounds.')
  }
  
  #  initialize matrix Kmat defining penalty term
  
  if (Wlambda > 0) {
    Kmat <- Wlambda*eval.penalty(Wbasis, WLfdobj)
  } else {
    Kmat  <- matrix(0,Wnbasis,Wnbasis)
  }
  
  #  load data into morphList
  
  morphList <- list(x=x, y=y, xlim=xlim, ylim=ylim, wt=wt, inactive=inactive,
                    Kmat=Kmat, Wlambda=Wlambda, Wbasis=Wbasis)
  
  #  Compute initial function and gradient values
  
  result <- fngrad2(cvec, morphList)   
  f    <- result$PENSSE 
  grad <- result$DPENSSE
  hmat <- result$D2PENSSE
  norm <- sqrt(sum(grad^2))
  
  #  compute initial badness of fit measures
  
  fold <- f
  cvecold <- cvec
  
  #  evaluate the initial update vector 
  
  pvec <- -solve(hmat,grad)
  
  #  initialize iteration status arrays
  
  iternum <- 0
  status  <- c(iternum, fold, norm)
  if (dbglev >= 1) {
    cat("\nIter.   PENSSE   Grad Length")
    cat("\n")
    cat(iternum)
    cat("        ")
    cat(round(status[2],4))
    cat("      ")
    cat(round(status[3],4))
  }
  iterhist <- matrix(0,iterlim+1,length(status))
  iterhist[1,] <- status
  
  #  ---------------------  Begin main iterations  ---------------
  
  STEPMAX <- 10
  itermax <- 20
  TOLX    <- 1e-10
  fold    <- f
  
  #  ---------------  beginning of optimization loop  -----------
  
  for (iter in 1:iterlim) {
    iternum <- iternum + 1
    #  line search
    result <- lnsrch(cvecold, fold, grad, pvec, fngrad2, morphList, 
                     STEPMAX, itermax, TOLX, dbglev)
    cvec <- result$x
    result <- fngrad2(cvec, morphList)   
    f    <- result$PENSSE 
    grad <- result$DPENSSE
    hmat <- result$D2PENSSE
    norm <- sqrt(sum(grad^2))
    status <- c(iternum, f, norm)
    iterhist[iter+1,] <- status
    if (dbglev >= 1) {
      cat("\n")
      cat(iternum)
      cat("        ")
      cat(round(status[2],4))
      cat("      ")
      cat(round(status[3],4))
    }
    #  test for convergence
    if (abs(f-fold) < conv) {
      break
    }
    if (f >= fold) { 
      warning('Current function value does not decrease fit.')
      break 
    }
    #  evaluate the update vector
    pvec   <- -solve(hmat,grad)
    cosangle <- -t(grad) %*% pvec/sqrt(sum(grad^2)*sum(pvec^2))
    if (cosangle < 0) {
      if (dbglev > 1) { 
        print('cos(angle) negative') 
      }
      pvec <- -grad
    }
    fold <- f
    cvecold <- cvec
  }
  
  #  construct output objects
  
  Wfdobj  <- fd(cvec, Wbasis)
  ywidth  <- ylim[2] - ylim[1]
  hraw    <- monfn(  x, Wfdobj)
  hmax    <- monfn(  xlim[2], Wfdobj)
  hmin    <- monfn(  xlim[1], Wfdobj)
  hwidth  <- hmax - hmin
  ysmth   <- (xlim[1] - hmin) + ywidth*hraw/hwidth
  cat("\n")
  
  return(list(Wfdobj=Wfdobj, f=f, grad=grad, hmat=hmat, norm=norm, 
              ysmth=ysmth, iternum=iternum, iterhist=iterhist)) 
}
#  ----------------------------------------------------------------

fngrad2 <- function(cvec, morphList) {
  #  extract data from morphList
  x        <- as.matrix(morphList$x)
  y        <- as.matrix(morphList$y)
  xlim     <- as.matrix(morphList$xlim)
  ylim     <- as.matrix(morphList$ylim)
  wt       <- as.matrix(morphList$wt)
  inactive <- as.matrix(morphList$inactive)
  Kmat     <- morphList$Kmat
  Wlambda  <- morphList$Wlambda
  Wbasis   <- morphList$Wbasis
  
  nobs     <- length(x)
  ywidth   <- ylim[2] - ylim[1]
  Wfdobj   <- fd(cvec, Wbasis)
  Wnbasis  <- Wbasis$nbasis
  #  get unnormalized function and gradient values
  hraw  <- as.matrix(monfn(  x, Wfdobj))
  Dhraw <- mongrad(x, Wfdobj)
  #  adjust functions and derivatives for normalization
  hmax    <- monfn(  xlim[2], Wfdobj) 
  hmin    <- monfn(  xlim[1], Wfdobj) 
  hwidth  <- hmax - hmin
  h  <- (xlim[1] - hmin) + ywidth*hraw/hwidth
  Dhmax   <- mongrad(xlim[2], Wfdobj)
  Dhmin   <- mongrad(xlim[1], Wfdobj)
  Dhwidth <- Dhmax - Dhmin
  Dh <- ywidth*(Dhraw*hwidth - hraw %*% Dhwidth)/hwidth^2
  res    <- y - h
  f <- mean(res^2*wt)
  #  if (required, gradient is computed
  temp   <- Dh*(wt %*% matrix(1,1,Wnbasis))
  grad   <- -2*t(temp) %*% res/nobs
  if (Wlambda > 0) {
    grad <- grad +         2 * Kmat %*% cvec
    f    <- f    + t(cvec) %*% Kmat %*% cvec
  }
  if (!is.null(inactive)) { 
    grad[inactive] <- 0 
  }
  #  if (required, hessian matrix is computed
  wtroot  <- sqrt(wt)
  temp <- Dh*(wtroot %*% matrix(1,1,Wnbasis))
  hmat <- 2*t(temp) %*% temp/nobs
  #  adjust for penalty
  if (Wlambda > 0) { 
    hmat <- hmat + 2*Kmat 
  }
  #  adjust for inactive coefficients
  if (!is.null(inactive)) {
    hmat[inactive,         ] <- 0
    hmat[,         inactive] <- 0
    hmat[inactive, inactive] <- diag(rep(1,Wnbasis)[inactive])
  }
  
  PENSSE   <- f
  DPENSSE  <- grad
  D2PENSSE <- hmat
  
  return(list(PENSSE=PENSSE, DPENSSE=DPENSSE, D2PENSSE=D2PENSSE))
  
}