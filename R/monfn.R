#  --------------------------------------------------------------

monfn <- function(x, Wfd, basislist=vector("list",JMAX)) {
#  evaluates a monotone function of the form
#            h(x) = [D^{-1} exp Wfd](x)
#  where  D^{-1} means taking the indefinite integral.
#  The interval over which the integration takes places is defined in
#  the basis object in Wfd.
#  Arguments:
#  X      ... argument values at which function and derivatives are evaluated
#  WFD    ... a functional data object
#  BASISLIST ... a list containing values of basis functions
#  Returns:
#  HVAL   ... value of h at input argument array X in first column.

#  Last modified 26 October 2003

  JMAX <- 15
  JMIN <- 11
  EPS  <- 1E-5

  coef  <- Wfd$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)
  if (ndim > 1 && coefd[2] != 1) stop("Wfd is not a single function")

  basisobj  <- Wfd$basis
  rangeval <- basisobj$rangeval

  #  set up first iteration

  width <- rangeval[2] - rangeval[1]
  JMAXP <- JMAX + 1
  h <- rep(1,JMAXP)
  h[2] <- 0.25
  #  matrix SMAT contains the history of discrete approximations to the
  #    integral
  smat <- matrix(0,JMAXP)
  #  array TVAL contains the argument values used in the approximation
  #  array FVAL contains the integral values at these argument values,
  #     rows corresponding to argument values
  #  the first iteration uses just the endpoints
  tval <- rangeval
  j   <- 1
  if (is.null(basislist[[j]])) {
      bmat <- getbasismatrix(tval, basisobj)
      basislist[[j]] <- bmat
  } else {
      bmat <- basislist[[j]]
  }
  fx   <- exp(bmat %*% coef)
  fval <- fx
  smat[1,]  <- width*apply(fx,2,sum)/2
  tnm <- 0.5

  #  now iterate to convergence
  for (j in 2:JMAX) {
    tnm  <- tnm*2
    del  <- width/tnm
    tj   <- seq(rangeval[1]+del/2, rangeval[2]-del/2, del)
    tval <- c(tval, tj)
    if (is.null(basislist[[j]])) {
        bmat <- getbasismatrix(tj, basisobj)
        basislist[[j]] <- bmat
    } else {
        bmat <- basislist[[j]]
    }
    fx   <- exp(bmat %*% coef)
    fval <- c(fval,fx)
    smat[j] <- (smat[j-1] + width*apply(fx,2,sum)/tnm)/2
    if (j >= JMIN) {
      ind <- (j-4):j
      result <- polintmat(h[ind],smat[ind],0)
      ss  <- result[[1]]
      dss <- result[[2]]
      if (all(abs(dss) < EPS*max(abs(ss)))) {
        # successful convergence
        # sort argument values and corresponding function values
        ordind <- order(tval)
        tval   <- tval[ordind]
        fval   <- fval[ordind]
        nx     <- length(tval)
        del    <- tval[2] - tval[1]
        fval   <- del*(cumsum(fval) - 0.5*(fval[1] + fval))
        hval   <- approx(tval, fval, x)$y
        return(hval)
      }
    }
    smat[j+1] <- smat[j]
    h[j+1]    <- 0.25*h[j]
  }
  stop(paste("No convergence after",JMAX," steps in MONFN"))
}
