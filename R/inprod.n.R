inprod.n <- function(fd1, fd2, Lfd1=0, Lfd2=0, penrng=rangeval, wtfd=NULL,
                     JMAX=12, EPS=1e-5) {

#  computes matrix of inner products of functions by numerical integration
#    using Romberg integration

#  Arguments:
#  FD1 and FD2  ...  these may be either functional data or basis function
#                    objects.  In the latter case, a functional data object
#                    is created from a basis function object by using the
#                    identity matrix as the coefficient matrix.
#                    Both functional data objects must be univariate.
#                    If inner products for multivariate objects are needed,
#                    use a loop and call inprod(fd1[i],fd2[i]).
#  LFD1 and LFD2  ...  order of derivatives for inner product for
#                    FD1 and FD2, respectively, or functional data objects
#                    defining linear differential operators
#  PENRNG ...  vector of length 2 giving the interval over which the
#               integration is to take place
#  WTFD   ...  functional data object defining a weight function
#  JMAX ...  maximum number of allowable iterations
#  EPS  ...  convergence criterion for relative error

#  Return:
#  A matrix of NREP1 by NREP2 of inner products for each possible pair
#  of functions.


#  Last modified 6 Feb 2001

  #  check arguments, and convert basis objects to functional data objects
  fdclass <- TRUE
  if (inherits(fd1, "fd") || inherits(fd1, "basis.fd")) {
    if (inherits(fd1, "basis.fd")) {
      coef1 <- diag(rep(1,fd1$nbasis))
      fd1 <- create.fd(coef1, fd1)
    } else coef1 <- getcoef(fd1)
  } else fdclass <- FALSE
  if (inherits(fd2, "fd") || inherits(fd2, "basis.fd")) {
    if (inherits(fd2, "basis.fd")) {
      coef2 <- diag(rep(1,fd2$nbasis))
      fd2 <- create.fd(coef2, fd2)
    } else coef2 <- getcoef(fd2)
  } else fdclass <- FALSE
  if (!fdclass) stop ("The two first arguments must be functional data objects")

  if (!is.null(wtfd)) {
    wtcoef <- getcoef(wtfd)
    if (dim(wtcoef)[[2]] != 1) stop("More than one weight function found")
  }

  #  determine NREP1 and NREP2, and check for common range
  coefd1 <- dim(coef1)
  coefd2 <- dim(coef2)
  ndim1  <- length(coefd1)
  ndim2  <- length(coefd2)
  if (ndim1 > 2 || ndim2 > 2) stop(
    "Functional data objects must be univariate")
  if (ndim1 > 1) nrep1 <- coefd1[2] else nrep1 <- 1
  if (ndim2 > 1) nrep2 <- coefd2[2] else nrep2 <- 1

  rangeval <- fd1$basis$rangeval
  #  set up first iteration
  width <- penrng[2] - penrng[1]
  JMAXP <- JMAX + 1
  h <- rep(1,JMAXP)
  h[2] <- 0.25
  s <- array(0,c(JMAXP,nrep1,nrep2))
  #  the first iteration uses just the endpoints
  fx1 <- eval.fd(penrng, fd1, Lfd1)
  fx2 <- eval.fd(penrng, fd2, Lfd2)
  if (is.null(wtfd)) {
    s[1,,]  <- width*crossprod(fx1,fx2)/2
  } else {
    wt <- c(eval.fd(penrng, wtfd))
    s[1,,]  <- width*crossprod(fx1,wt*fx2)/2
  }
  tnm <- 0.5
  j <- 1
  #print(j)
  #print(round(s[j,,],2))
  cat('\nNumerical integration iterations:  .')

  #  now iterate to convergence
  for (j in 2:JMAX) {
    tnm <- tnm*2
    del <- width/tnm
    x   <- seq(penrng[1]+del/2, penrng[2]-del/2, del)
    fx1 <- eval.fd(x, fd1, Lfd1)
    fx2 <- eval.fd(x, fd2, Lfd2)
    if (is.null(wtfd)) {
      s[j,,] <- (s[j-1,,] + width*crossprod(fx1,fx2)/tnm)/2
    } else {
      wt <- c(eval.fd(x, wtfd))
      s[j,,] <- (s[j-1,,] + width*crossprod(fx1,wt*fx2)/tnm)/2
    }
    cat('.')
    #print(j)
    #print(round(s[j,,],2))
    if (j >= 5) {
      ind <- (j-4):j
      result <- polintmat(h[ind],s[ind,,],0)
      ss  <- result[[1]]
      #print(round(log10(diag(ss)),2))
      #image(ss)
      #title(paste("J =",j))
      #text(locator(1),"")
      dss <- result[[2]]
      if (all(abs(dss) < EPS*max(abs(ss)))) return(ss)
    }
    s[j+1,,] <- s[j,,]
    h[j+1]   <- 0.25*h[j]
  }
  warning(paste("No convergence after",JMAX," steps in INPROD"))
}
