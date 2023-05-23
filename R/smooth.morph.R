smooth.morph <- function(x, y, ylim, WfdPar,   
                         conv=1e-4, iterlim=20, dbglev=0) {
  #  SMOOTH_MORPH smooths the relationship of Y to X 
  #  by fitting a monotone fn.  f(x) <- b_0 + b_1 D^{-1} exp W(t)
  #     where  W  is a function defined over the same range as X,
  #  W + ln b_1 <- log Df and w <- D W <- D^2f/Df.
  #  b_0 and b_1 are chosen so that values of f
  #  are within the interval [ylim[1],ylim[2]].
  #  The fitting criterion is penalized mean squared stop:
  #    PENSSE(lambda) <- \sum [y_i - f(t_i)]^2 +
  #                     \lambda * \int [L W]^2 
  #  W(x) is expanded by the basis in functional data object Wfdobj.
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
  #  DBGLEV    Controls the level of output on each iteration.  If 0,
  #               no output, if 1, output at each iteration, if higher, 
  #               output at each line search iteration. 1 by default.
  #
  #  Returns are:
  #  WFD       Functional data object for W.  It's coefficient vector
  #               contains the optimized coefficients.
  
  #  last modified 17 May 2023 by Jim Ramsay
  
  #  number of observations and weights on x-values
  
  nobs <- length(x)        
  wt   <- matrix(1,nobs,1)
  wt[nobs] <- 10
  
  #  -----------------------------------------------------
  #                  Check arguments
  #  -----------------------------------------------------
  
  # check consistency of x and y and convert to column matrices
  
  if (length(y) != nobs) {
    stop('Arguments X and Y are not of same length')
  }
  
  if (!is.matrix(x)) x <- matrix(x,nobs,1)
  if (!is.matrix(y)) y <- matrix(y,nobs,1)
  
  #  check ylim for not being numeric
  
  if (!is.numeric(ylim)) {
    print("The third argument ylim is not numeric.  This argument
          should be a vector of length 2 containing the boundaries of
          the target interval.")
    stop("This argument has been added to allow morphs between two unequal invervals.")
  }
  
  #  check ylim for not having two strictly increasing numbers
  
  if (length(ylim) != 2 || ylim[1] >= ylim[2])
    stop("Argument ylim does not containing two strictly increasing numbers.")
  
  #  check WfdPar
  
  if (!is.fdPar(WfdPar)) {
    if (is.fd(WfdPar)) {
      WfdPar <- fdPar(WfdPar)
    } else {
      stop(paste("WFDPAR is not a functional parameter object,", 
                 "and not a functional data object."))
    }
  }
  
  #  -----------------------------------------------------
  #                  Initialize optimization
  #  -----------------------------------------------------
  
  #      extract information from WfdPar
  
  Wfdobj   <- WfdPar$fd
  cvec     <- Wfdobj$coef   #  initial coefficients
  Wbasis   <- Wfdobj$basis      #  basis for Wfdobj
  Wnbasis  <- Wbasis$nbasis      #  no. basis functions
  Wrange   <- Wbasis$rangeval
  Wtype    <- Wbasis$type
  xlim     <- Wbasis$rangeval
  WLfdobj  <- int2Lfd(WfdPar$Lfd)
  Wlambda  <- WfdPar$lambda
  if (any(wt < 0)) { 
    stop("One or more weights are negative.") 
  }
  
  #  check that values in x are within limits in xlim
  
  if (abs(x[1] - xlim[1]) > 1e-7 || abs(x[nobs] - xlim[2]) > 1e-7) {
    print("Argument vector x is out of range:")
    stop("Values are out of bounds by more than 1e-7")
  }  else {
    x[1   ] <- xlim[1]
    x[nobs] <- xlim[2]
  }
  
  #  transform coefficients to zero column sum 
  
  Zmat <- fda::zerobasis(length(cvec))
  bvec <- t(Zmat) %*% cvec
  cvec <- Zmat %*% bvec
  
  #  initialize matrix Kmat defining penalty term
  
  if (Wlambda > 0) {
    Kmat <- Wlambda*eval.penalty(Wbasis, WLfdobj)
  } else {
    Kmat  <- matrix(0,Wnbasis,Wnbasis)
  }
  
  #  load objects into morphList, used in fngrad_morph
  
  morphList <- list(x=x, y=y, xlim=xlim, ylim=ylim, wt=wt, Kmat=Kmat, 
                    Wlambda=Wlambda, Wbasis=Wbasis)
  
  #  -----------------------------------------------------
  #       Compute initial function and gradient values
  #  -----------------------------------------------------
  
  fnList <- fngrad_morph(bvec, morphList, Zmat)   
  f    <- fnList$f
  grad <- fnList$grad 
  hmat <- fnList$hmat
  norm <- fnList$norm
  
  #  compute initial badness of fit measures
  
  fold    <- f
  cvecold <- cvec
  
  #  evaluate the initial update vector 
  
  pvec <- -solve(hmat,grad)
  
  #  -----------------------------------------------------
  #  initialize iteration status arrays
  #  -----------------------------------------------------
  
  iternum <- 0
  status <- c(iternum, fold, norm)
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
  iterhist[1,]  <- status
  
  STEPMAX <- 10
  itermax <- 20
  TOLX    <- 1e-10
  fold    <- f
  
  #  -----------------------------------------------------
  #  --------  beginning of optimization loop  -----------
  #  -----------------------------------------------------
  
  for (iter in 1:iterlim) {
    iternum <- iternum + 1
    #  line search
    bvecold <- t(Zmat) %*% cvecold
    lnsrchList <- lnsrch_morph(bvecold, fold, grad, pvec, fngrad_morph, 
                                morphList, Zmat, STEPMAX, itermax, TOLX, 
                                dbglev)
    bvec   <- lnsrchList$x
    check  <- lnsrchList$check
    fnList <- fngrad_morph(bvec, morphList, Zmat)   
    f    <- fnList$f
    grad <- fnList$grad 
    hmat <- fnList$hmat
    norm <- fnList$norm
    cvec <- Zmat %*% bvec
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
    if (abs(f-fold) < conv) {
      break
    }
    if (f >= fold) { 
      warning("Current function value does not decrease fit.")
      cat("\n")
      break 
    }
    #  evaluate the update vector
    pvec <- -solve(hmat,grad)
    cosangle <- -t(grad) %*% pvec/sqrt(sum(grad^2)*sum(pvec^2))
    if (cosangle < 0) {
      if (dbglev > 1) {
        print("cos(angle) negative") 
      }
      pvec <- -grad
    }
    fold <- f
    cvecold <- cvec
  }
  
  #  -----------------------------------------------------
  #    construct output objects and return in a list
  #  -----------------------------------------------------
  
  Wfdobj  <- fd(cvec, Wbasis)
  
  xfine   <- as.matrix(seq(xlim[1],xlim[2],len=101))
  ywidth  <- ylim[2] - ylim[1]
  hfine   <- matrix(monfn(xfine, Wfdobj), 101, 1)
  hmax    <- monfn(xlim[2], Wfdobj)
  hmin    <- monfn(xlim[1], Wfdobj)
  hwidth  <- hmax - hmin
  hfine   <- (ylim[1] - hmin) + hfine*(ywidth/hwidth)
  
  return(list(Wfdobj=Wfdobj, f=f, grad=grad, hmat=hmat, 
              norm=norm, hfine=hfine, 
              iternum=iternum, iterhist=iterhist))
  
}

#  ------------------------------------------------------------------------

fngrad_morph <- function(bvec, morphList, Zmat) {
  #  -----------------------------------------------------
  #  extract data from morphList
  #  -----------------------------------------------------
  cvec     <- Zmat %*% bvec
  x        <- morphList$x
  y        <- morphList$y
  xlim     <- morphList$xlim
  ylim     <- morphList$ylim
  wt       <- morphList$wt
  Kmat     <- morphList$Kmat
  Wlambda  <- morphList$Wlambda
  Wbasis   <- morphList$Wbasis
  Wnbasis  <- Wbasis$nbasis
  Wfdobj   <- fd(cvec, Wbasis)
  #  -----------------------------------------------------
  #             compute fitting criterion f
  #  -----------------------------------------------------
  #  compute unnormalized monotone function values hraw
  
  nobs  <- length(x)
  hraw  <- matrix(monfn(  x, Wfdobj),nobs,1)
  #  adjust functions and derivatives for normalization
  hmax    <- monfn(  xlim[2], Wfdobj) 
  hmin    <- monfn(  xlim[1], Wfdobj) 
  #  width of hraw
  hwidth  <- hmax - hmin
  #  width of target interval
  ywidth  <- ylim[2] - ylim[1]
  #  normalized h varying horizontally over base interval and 
  #  vertically over target interval
  h   <- (ylim[1] - hmin) + hraw*(ywidth/hwidth)
  #  compute least squares fitting criterion
  res <- y - h
  f   <- mean(res^2*wt)
  #  -----------------------------------------------------
  #             compute fitting gradient grad
  #  -----------------------------------------------------
  #  un-normalized partial derivative of un-normalized h Dh
  Dhraw   <- matrix(mongrad(x, Wfdobj),nobs,Wnbasis)
  #  range of un-normalized gradient
  Dhmax   <- matrix(mongrad(xlim[2], Wfdobj), 1, Wnbasis)
  Dhmin   <- matrix(mongrad(xlim[1], Wfdobj), 1, Wnbasis)
  Dhwidth <- Dhmax - Dhmin
  #  normalized gradient
  Dh    <- ywidth*(Dhraw*hwidth - hraw %*% Dhwidth)/hwidth^2
  #  gradient of fitting function is computed
  temp  <- Dh*(wt %*% matrix(1,1,Wnbasis))
  grad  <- -2*t(temp) %*% res/nobs
  #  apply regularization if needed
  if (Wlambda > 0) {
    grad <- grad +           2 * Kmat  %*%  cvec
    f    <- f    + t(cvec)  %*%  Kmat  %*%  cvec
  }
  #  map parameter space into fitting space
  grad <- t(Zmat) %*% grad
  norm <- sqrt(sum(grad^2))   #  gradient norm
  #  -----------------------------------------------------
  #        compute fitting Hessian hmat
  #  -----------------------------------------------------
  wtroot  <- sqrt(wt)
  temp <- Dh * (wtroot %*% matrix(1,1,Wnbasis))
  hmat <- 2*t(temp) %*% temp/nobs
  #  apply regularization if needed
  if (Wlambda > 0) { 
    hmat <- hmat + 2*Kmat 
  }
  #  map parameter hessian into fitting hessian
  hmat <- t(Zmat) %*% hmat %*% Zmat
  
  return(list(f=f, grad=grad, hmat=hmat, norm=norm, h=h, Dh=Dh))
  
}

#  ------------------------------------------------------------------------

lnsrch_morph <- function(xold, fold, g, p, func, dataList, 
                          Zmat, stpmax, itermax=20, TOLX=1e-10, dbglev=0) {
  #  LNSRCH computes an approximately optimal parameter vector X given
  #  an initial parameter vector XOLD, an initial function value FOLD
  #  and an initial gradient vector G.  
  
  #  Arguments:
  #  XOLD      Initial parameter value
  #  FOLD      Initial function value
  #  G         Initial gradient value
  #  P         Search direction vector
  #  FUNC      Function object computing a function value and gradient
  #              vector
  #  DATALIST  List object used in function object FUNC
  #  STPMAX    Maximum step size
  #  ITERMAX   Maximum number of iterations
  #  TOLX      Tolerance for stop
  #  DBGLEV    Debugging output value:  none if 0, function value if 1,
  #              if greater than one, current step value, slope and
  #              function value.
  
  #  Last modified 12 February 2022
  
  #  set initial constants
  n <- length(xold)
  check <- FALSE
  f2    <- 0
  alam2 <- 0
  ALF   <- 1e-4
  psum  <- sqrt(sum(p^2))
  #  scale if attempted step is too big
  if (psum > stpmax) {
    p <- p*(stpmax/psum)
  }
  #  compute slope
  slope <- sum(g*p)
  # if (dbglev > 1) {
  #   status <- c(0, 0, slope, fold)
  #   cat("#10.f #10.4f #10.4f #10.4f\n", status)
  # }
  # check that initial slope is negative
  if (slope >= 0) {
    stop("Initial slope not negative.")
  }
  # compute lambdamin
  test <- 0
  for (i in 1:n) {
    temp <- abs(p[i])/max(abs(xold[i]),1)
    if (temp > test) {
      test <- temp
    }
  }
  alamin <- TOLX/test
  #  always try full Newton step first with step size 1
  alam   <- 1
  #  start of iteration loop
  iter <- 0
  while (iter <= itermax) {
    iter <- iter + 1
    x <- xold + alam*p
    #  -------------  function evaluation  -----------
    funcList <- func(x, dataList, Zmat)
    f    <- funcList$f
    gtmp <- funcList$grad
    slp <- sum(gtmp*p)
    # if (dbglev > 1) {
    #   status <- [iter, alam, slp, f]
    #   cat("#10.f #10.4f #10.4f #10.4f\n", status)
    # }
    #  -----------------------------------------------
    #  convergence on x change.
    if (alam < alamin) {
      x <- xold
      check <- TRUE
      return(list(x=x, check=check))
    } else {
      #  sufficient function decrease
      if (f <= fold + ALF*alam*slope) {
        return(list(x=x, check=check))
      }
      #  backtrack
      if (alam == 1) {
        #  first time
        tmplam <- -slope/(2*(f-fold-slope))
      } else {
        #  subsequent backtracks
        rhs1 <- f  - fold - alam *slope
        rhs2 <- f2 - fold - alam2*slope
        a <- (rhs1/alam^2 - rhs2/alam2^2)/(alam-alam2)
        b <- (-alam2*rhs1/alam^2 + alam*rhs2/(alam*alam2))/(alam-alam2)
        if (a == 0) {
          tmplam <- -slope/(2*b)
        } else {
          disc <- b^2 - 3*a*slope
          if (disc < 0) {
            tmplam <- 0.5*alam
          } else {
            if (b <= 0) {
              tmplam <- (-b+sqrt(disc))/(3*a)
            } else {
              tmplam <- -slope/(b+sqrt(disc))
            }
          }
          if (tmplam > 0.5*alam) {
            # lambda <= 0.5 lambda1
            tmplam <- 0.5*alam
          }
        }
      }
      alam2 <- alam
      f2    <- f
      #  lambda > 0.1 lambda1
      alam <- max(tmplam, 0.1*alam)
    }
    #  try again
    
  }
  
  return(list(x=x, check=check))
  
}


