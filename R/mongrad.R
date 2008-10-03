#  --------------------------------------------------------------------------

mongrad <- function(x, Wfdobj, basislist=vector("list",JMAX)) {
#  Evaluates the gradient with respect to the coefficients in Wfdobj
#     of a monotone function of the form
#            h(x) = [D^{-1} exp Wfdobj](x)
#  where  D^{-1} means taking the indefinite integral.
#  The interval over which the integration takes places is defined in
#  the basisfd object in Wfdobj.
#  Arguments:
#  X      ... argument values at which function and derivatives are evaluated
#  WFDOBJ ... a functional data object
#  BASISLIST ... a list containing values of basis functions
#  Returns:
#  GVAL   ... value of gradient at input values in X.

#  Last modified 19 December 2007

  JMAX <- 15
  JMIN <- 11
  EPS  <- 1E-5

  coef  <- Wfdobj$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)
  if (ndim > 1 && coefd[2] != 1) stop("Wfdobj is not a single function")

  basisfd  <- Wfdobj$basis
  rangeval <- basisfd$rangeval
  nbasis   <- basisfd$nbasis
  onebas   <- rep(1,nbasis)
  width    <- rangeval[2] - rangeval[1]

  #  set up first iteration

  JMAXP <- JMAX + 1
  h <- rep(1,JMAXP)
  h[2] <- 0.25
  #  matrix SMAT contains the history of discrete approximations to the
  #    integral
  smat <- matrix(0,JMAXP,nbasis)
  #  array TVAL contains the argument values used in the approximation
  #  array FVAL contains the integral values at these argument values,
  #     rows corresponding to argument values
  #  the first iteration uses just the endpoints
  j   <- 1
  tval <- rangeval
  if (is.null(basislist[[j]])) {
      bmat <- getbasismatrix(tval, basisfd)
      basislist[[j]] <- bmat
  } else {
      bmat <- basislist[[j]]
  }
  fx   <- exp(bmat %*% coef)
  fval <- outer(c(fx),onebas)*bmat
  smat[1,]  <- width*apply(fval,2,sum)/2
  tnm <- 0.5

  #  now iterate to convergence
  for (j in 2:JMAX) {
    tnm  <- tnm*2
    del  <- width/tnm
    tj   <- seq(rangeval[1]+del/2, rangeval[2]-del/2, del)
    tval <- c(tval, tj)
    if (is.null(basislist[[j]])) {
        bmat <- getbasismatrix(tj, basisfd)
        basislist[[j]] <- bmat
    } else {
        bmat <- basislist[[j]]
    }
    fx   <- exp(bmat %*% coef)
    gval <- outer(c(fx),onebas)*bmat
    fval <- rbind(fval,gval)
    smat[j,] <- (smat[j-1,] + width*apply(fval,2,sum)/tnm)/2
    if (j >= max(c(5,JMIN))) {
      ind <- (j-4):j
      result <- polintmat(h[ind],smat[ind,],0)
      ss  <- result[[1]]
      dss <- result[[2]]
      if (all(abs(dss) < EPS*max(abs(ss))) || j == JMAX) {
        # successful convergence
        # sort argument values and corresponding function values
        ordind <- order(tval)
        tval   <- tval[ordind]
        fval   <- as.matrix(fval[ordind,])
        # set up partial integral values
        lval   <- outer(rep(1,length(tval)),fval[1,])
        del    <- tval[2] - tval[1]
        fval   <- del*(apply(fval,2,cumsum) - 0.5*(lval + fval))
        gval   <- matrix(0,length(x),nbasis)
        for (i in 1:nbasis) gval[,i] <- approx(tval, fval[,i], x)$y
        return(gval)
      }
    }
    smat[j+1,] <- smat[j,]
    h[j+1]     <- 0.25*h[j]
  }
  #stop(paste("No convergence after",JMAX," steps in MONGRAD"))
}
