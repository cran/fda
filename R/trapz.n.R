trapz.n <- function(fd, Lfd=0, rng=rangeval, wtfd=NULL,
                     JMAX=12, EPS=1e-5) {

#  computes matrix of integrals of functions or values of a
#  differential operator applied to the functions by numerical integration
#    using Romberg integration

#  Arguments:
#  FD  ...  Either functional data or basis function
#           object.  In the latter case, a functional data object
#           is created from a basis function object by using the
#           identity matrix as the coefficient matrix.
#           The functional data objects must be univariate.
#  LFD ...  Order of derivative for integration for
#           FD, respectively, or functional data objects
#           defining linear differential operators
#  RNG ...  vector of length 2 giving the interval over which the
#           integration is to take place
#  WTFD ... functional data object defining a weight function
#  JMAX ... maximum number of allowable iterations
#  EPS  ... convergence criterion for relative error

#  Return:
#  A vector of length NREP of integralss for each function.

#  Last modified 6 Feb 2001

  #  check arguments, and convert basis objects to functional data objects
  fdclass <- TRUE
  if (inherits(fd, "fd") || inherits(fd, "basis.fd")) {
    if (inherits(fd, "basis.fd")) {
      coef <- diag(rep(1,fd$nbasis))
      fd <- create.fd(coef, fd)
    } else coef <- getcoef(fd)
  } else fdclass <- FALSE
  if (!fdclass) stop ("The first argument must be functional data objects")

  if (!is.null(wtfd)) {
    wtcoef <- getcoef(wtfd)
    if (dim(wtcoef)[[2]] != 1) stop("More than one weight function found")
  }

  #  determine NREP1 and NREP2, and check for common range
  coefd <- dim(coef)
  ndim  <- length(coefd)
  if (ndim > 2) stop("Functional data object must be univariate")
  if (ndim > 1) nrep <- coefd[2] else nrep <- 1

  rangeval <- fd$basis$rangeval
  #  set up first iteration
  width <- rng[2] - rng[1]
  JMAXP <- JMAX + 1
  h <- rep(1,JMAXP)
  h[2] <- 0.25
  s <- matrix(0,c(JMAXP,nrep))
  #  the first iteration uses just the endpoints
  fx <- eval.fd(rng, fd, Lfd)
  if (is.null(wtfd)) {
    s[1,]  <- width*fx/2
  } else {
    wt <- c(eval.fd(penrng, wtfd))
    s[1,]  <- width*(wt*fx)/2
  }
  tnm <- 0.5
  j <- 1
  #print(j)
  #print(round(s[j,],2))
  cat('\n.')

  #  now iterate to convergence
  for (j in 2:JMAX) {
    tnm <- tnm*2
    del <- width/tnm
    x   <- seq(rng[1]+del/2, rng[2]-del/2, del)
    fx <- eval.fd(x, fd, Lfd)
    if (is.null(wtfd)) {
      s[j,] <- (s[j-1,] + width*fx/tnm)/2
    } else {
      wt <- c(eval.fd(x, wtfd))
      s[j,] <- (s[j-1,] + width*(wt*fx)/tnm)/2
    }
    cat('.')
    #print(j)
    #print(round(s[j,,],2))
    if (j >= 5) {
      ind <- (j-4):j
      result <- polintmat(h[ind],s[ind,],0)
      ss  <- result[[1]]
      dss <- result[[2]]
      if (all(abs(dss) < EPS*max(abs(ss)))) return(list(ss, s[j,]))
    }
    s[j+1,] <- s[j,]
    h[j+1]   <- 0.25*h[j]
  }
  warning(paste("No convergence after",JMAX," steps in INPROD"))
}
