smooth.monotone <- function(x, y, wt=rep(1,nobs), WfdParobj, zmat=matrix(1,nobs,1),
                            conv=.0001, iterlim=20, active=c(FALSE,rep(TRUE,ncvec-1)),
                            dbglev=1) {
#  Smooths the relationship of Y to X using weights in WT by fitting a
#     monotone function of the form
#                   f(x) = b_0 + b_1 D^{-1} exp W(x)
#     where  W  is a function defined over the same range as X,
#                 W + ln b_1 = log Df and w = D W = D^2f/Df.
#  The constant term b_0 in turn can be a linear combinations of covariates:
#                         b_0 = zmat * c.
#  The fitting criterion is penalized mean squared error:
#    PENSSE(lambda) = \sum w_i[y_i - f(x_i)]^2 +
#                     \lambda * \int [L W(x)]^2 dx
#  where L is a linear differential operator defined in argument Lfdobj,
#  and w_i is a positive weight applied to the observation.
#  The function W(x) is expanded by the basis in functional data object
#    Wfdobj.   The coefficients of this expansion are called "coefficients"
#    in the comments, while the b's are called "regression coefficients"

#  Arguments:
#  X       ...  vector of argument values
#  Y       ...  vector of function values to be fit
#  WT      ...  a vector of weights
#  WFDPAROBJ  ...  functional parameter object for W(x).  The coefficient array
#               for WFDPAROBJ$FD has a single column, and these are the
#               starting values for the iterative minimization of mean squared error.
#  ZMAT    ...  a matrix of covariate values for the constant term.
#               It defaults to a column of one's
#  CONV    ...  convergence criterion, 0.0001 by default
#  ITERLIM ...  maximum number of iterations, 20 by default
#  ACTIVE  ...  vector of 1's and 0's indicating which coefficients
#               are to be optimized (1) or remain fixed (0).  All values
#               are 1 by default, except that if a B-spline basis is used,
#               the first value is set to 0.
#  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
#               no output, if 1, output at each iteration, if higher, output
#               at each line search iteration. 1 by default.

#  Returns a list containing:
#  WFDOBJ    ...  functional data object for W(x).  Its coefficients are
#                   those that optimize fit.
#  BETA      ...  final regression coefficient values
#  FLIST     ...  A list containing the final function value, gradient and
#                 gradient norm.
#  ITERNUM   ...  number of iterations
#  ITERHIST  ...  ITERNUM+1 by 5 array containing iteration history

#  Last modified 26 October 2005

#  initialize some arrays

  x      <- as.vector(x)
  nobs   <- length(x)         #  number of observations

#  the starting values for the coefficients are in FD object WFDOBJ

  if (!(inherits(WfdParobj, "fdPar"))) stop(
		"Argument WFDPAROBJ is not a functional data object.")

  Wfdobj   <- WfdParobj$fd
  Lfdobj   <- WfdParobj$Lfd
  basisobj <- Wfdobj$basis     #  basis for W(x)
  nbasis   <- basisobj$nbasis  #  number of basis functions
  type     <- basisobj$type
  cvec     <- Wfdobj$coefs
  ncvec    <- length(cvec)

#  check some arguments

  if (any(wt < 0))  stop("One or more weights are negative.")
  if (all(wt == 0)) stop("All weights are zero.")
  zdim <- dim(zmat)
  if (zdim[1] != nobs) stop(
    "First dimension of ZMAT not correct.")

#  set up some variables

  ncov   <- zdim[2]   #  number of covariates
  ncovp1 <- ncov + 1  #  index for regression coef. for monotone fn.
  wtroot <- sqrt(wt)
  wtrtmt <- wtroot %*% matrix(1,1,ncovp1)
  yroot  <- y*wtroot
  climit <- c(-100*rep(1,nbasis), 100*rep(1,nbasis))
  inact  <- !active   #  indices of inactive coefficients

#  set up cell for storing basis function values

  JMAX <- 15
  basislist <- vector("list", JMAX)

#  initialize matrix Kmat defining penalty term

  if (WfdParobj$lambda > 0) Kmat <- WfdParobj$lambda*getbasispenalty(basisobj, Lfdobj)
  else            Kmat <- NULL

#  Compute initial function and gradient values

  result <- fngrad.smooth.monotone(y, x, zmat, wt, Wfdobj, WfdParobj$lambda,
                                   Kmat, inact, basislist)
  Flist  <- result[[1]]
  beta   <- result[[2]]
  Dyhat  <- result[[3]]

#  compute the initial expected Hessian

  hessmat <- hesscal.smooth.monotone(beta, Dyhat, wtroot, WfdParobj$lambda, Kmat, inact)

#  evaluate the initial update vector for correcting the initial cvec

  result   <- linesearch.smooth.monotone(Flist, hessmat, dbglev)
  deltac   <- result[[1]]
  cosangle <- result[[2]]

#  initialize iteration status arrays

  iternum <- 0
  status  <- c(iternum, Flist$f, Flist$norm, beta)
  if (dbglev >= 1) {
    cat("\nIter.   PENSSE   Grad Length Intercept   Slope")
    cat("\n")
    cat(iternum)
    cat("        ")
    cat(round(status[2],4))
    cat("      ")
    cat(round(status[3],4))
    cat("      ")
    cat(round(beta[1],4))
    cat("      ")
    cat(round(beta[ncovp1],4))
  }

  iterhist <- matrix(0,iterlim+1,length(status))
  iterhist[1,]  <- status
  if (iterlim == 0)
    return ( list( "Wfdobj" = Wfdobj, "beta" = beta, "Flist" = Flist,
                 "iternum" = iternum, "iterhist" = iterhist ) )

#  -------  Begin iterations  -----------

  MAXSTEPITER <- 10
  MAXSTEP     <- 100
  trial       <- 1
  reset       <- FALSE
  linemat     <- matrix(0,3,5)
  betaold     <- beta
  cvecold     <- cvec
  Foldlist    <- Flist
  dbgwrd      <- dbglev >= 2

  for (iter in 1:iterlim)
  {
     iternum <- iternum + 1
     #  initialize logical variables controlling line search
     dblwrd <- c(FALSE,FALSE)
     limwrd <- FALSE
     stpwrd <- FALSE
     ind    <- 0
     ips    <- 0
     #  compute slope at 0 for line search
     linemat[2,1] <- sum(deltac*Flist$grad)
     #  normalize search direction vector
      sdg     <- sqrt(sum(deltac^2))
      deltac  <- deltac/sdg
      dgsum   <- sum(deltac)
      linemat[2,1] <- linemat[2,1]/sdg
      # initialize line search vectors
      linemat[,1:4] <- outer(c(0, linemat[2,1], Flist$f),rep(1,4))
      stepiter <- 0
      if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,1]),6)))
      }
      #  return with error condition if initial slope is nonnegative
      if (linemat[2,1] >= 0) {
        if (dbgwrd >= 2) print("Initial slope nonnegative.")
        ind <- 3
        break
      }
      #  return successfully if initial slope is very small
      if (linemat[2,1] >= -1e-7) {
        if (dbglev >= 2) print("Initial slope too small")
        ind <- 0
        break
      }
      #  first step set to trial
      linemat[1,5]  <- trial
      #  Main iteration loop for linesearch
      cvecnew <- cvec
      Wfdnew  <- Wfdobj
      for (stepiter in 1:MAXSTEPITER)
      {
      #  ensure that step does not go beyond limits on parameters
        limflg  <- FALSE
        #  check the step size
        result <-
              stepchk(linemat[1,5], cvec, deltac, limwrd, ind,
                      climit, active, dbgwrd)
        linemat[1,5] <- result[[1]]
        ind          <- result[[2]]
        limwrd       <- result[[3]]
        if (linemat[1,5] <= 1e-7)
        {
          #  Current step size too small ... terminate
          if (dbglev >= 2) {
            print("Stepsize too small")
            print(avec[5])
          }
          if (limflg) ind <- 1 else ind <- 4
          break
        }
        #  compute new function value and gradient
        cvecnew <- cvec + linemat[1,5]*deltac
        Wfdnew$coefs <- as.matrix(cvecnew)
        result  <- fngrad.smooth.monotone(y, x, zmat, wt, Wfdnew, WfdParobj$lambda,
                                          Kmat, inact, basislist)
        Flist   <- result[[1]]
        beta    <- result[[2]]
        Dyhat   <- result[[3]]
        linemat[3,5] <- Flist$f
        #  compute new directional derivative
        linemat[2,5] <- sum(deltac*Flist$grad)
        if (dbglev >= 2) {
          cat("\n")
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,5]),6)))
        }
        #  compute next line search step, also test for convergence
        result  <- stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbglev)
        linemat <- result[[1]]
        ips     <- result[[2]]
        ind     <- result[[3]]
        dblwrd  <- result[[4]]
        trial   <- linemat[1,5]
        #  ind == 0  mean convergence
        if (ind == 0 | ind == 5) break
     }
     #  end iteration loop
     cvec    <- cvecnew
     Wfdobj  <- Wfdnew
     #  check that function value has not increased
     if (Flist$f > Foldlist$f) {
        # if it has, terminate iterations with a message
        if (dbglev >= 2) {
          cat("Criterion increased: ")
          cat(format(round(c(Foldlist$f, Flist$f),4)))
          cat("\n")
        }
        #  reset parameters and fit
        beta         <- betaold
        cvec         <- cvecold
        Wfdobj$coefs <- cvec
        Flist        <- Foldlist
        deltac       <- -Flist$grad
        if (reset) {
          # This is the second time in a row that
          #  this has happened ... quit
          if (dbglev >= 2) cat("Reset twice, terminating.\n")
          return ( list( "Wfdobj" = Wfdobj, "beta" = beta, "Flist" = Flist,
                         "iternum" = iternum, "iterhist" = iterhist ) )
        } else {
          reset <- TRUE
        }
     } else {
       if (abs(Foldlist$f - Flist$f) < conv) {
	       cat("\n")
	       break
       }
       cvecold  <- cvec
       betaold  <- beta
       Foldlist <- Flist
       hessmat  <- hesscal.smooth.monotone(beta, Dyhat, wtroot, WfdParobj$lambda, Kmat, inact)
       #  udate the line search direction
       result   <- linesearch.smooth.monotone(Flist, hessmat, dbglev)
       deltac   <- result[[1]]
       cosangle <- result[[2]]
       reset    <- FALSE
     }
     #  store iteration status
     status <- c(iternum, Flist$f, Flist$norm, beta)
     iterhist[iter+1,] <- status
     if (dbglev >= 1) {
        cat("\n")
        cat(iternum)
        cat("        ")
        cat(round(status[2],4))
        cat("      ")
        cat(round(status[3],4))
        cat("      ")
        cat(round(beta[1],4))
        cat("      ")
        cat(round(beta[ncovp1],4))
     }
  }
  return ( list( "Wfdobj" = Wfdobj, "beta" = beta, "Flist" = Flist,
                 "iternum" = iternum, "iterhist" = iterhist ) )
}

#  ----------------------------------------------------------------

linesearch.smooth.monotone <- function(Flist, hessmat, dbglev)
{
  deltac   <- -symsolve(hessmat,Flist$grad)
  cosangle <- -sum(Flist$grad*deltac)/sqrt(sum(Flist$grad^2)*sum(deltac^2))
  if (dbglev >= 2) {
    cat(paste("\nCos(angle) =",format(round(cosangle,4))))
    if (cosangle < 1e-7) {
      if (dbglev >=2)  cat("\nCosine of angle too small\n")
      deltac <- -Flist$grad
    }
  }
  return(list(deltac, cosangle))
}

#  ----------------------------------------------------------------

fngrad.smooth.monotone <- function(y, x, zmat, wt, Wfdobj, lambda,
                                   Kmat, inact, basislist)
{
  ncov   <- ncol(zmat)
  ncovp1 <- ncov + 1
  nobs   <- length(x)
  cvec   <- Wfdobj$coefs
  nbasis <- length(cvec)
  h      <- monfn(x, Wfdobj, basislist)
  Dyhat  <- mongrad(x, Wfdobj, basislist)
  xmat   <- cbind(zmat,h)
  Dxmat  <- array(0,c(nobs,ncovp1,nbasis))
  Dxmat[,ncovp1,] <- Dyhat
  wtroot <- sqrt(wt)
  wtrtmt <- wtroot %*% matrix(1,1,ncovp1)
  yroot  <- y*wtroot
  xroot  <- xmat*wtrtmt
  #  compute regression coefs.
  beta   <- lsfit(xmat, y, wt, int=FALSE)$coef
  #  update fitted values
  yhat   <- xmat %*% beta
  #  update residuals and function values
  res    <- y - yhat
  f      <- mean(res^2*wt)
  grad   <- matrix(0,nbasis,1)
  for (j in 1:nbasis) {
    Dxroot <- Dxmat[,,j]*wtrtmt
    yDx <- crossprod(yroot,Dxroot) %*% beta
    xDx <- crossprod(xroot,Dxroot)
    grad[j] <- crossprod(beta,(xDx+t(xDx))) %*% beta - 2*yDx
  }
  grad <- grad/nobs
  if (lambda > 0) {
    grad <- grad +         2 * Kmat %*% cvec
    f    <- f    + t(cvec) %*% Kmat %*% cvec
  }
  if (any(inact)) grad[inact] <- 0
  norm <- sqrt(sum(grad^2)) #  gradient norm
  Flist <- list("f"=f,"grad"=grad,"norm"=norm)
  return(list(Flist, beta, Dyhat))
}

#  ----------------------------------------------------------------

hesscal.smooth.monotone <- function(beta, Dyhat, wtroot, lambda,
                                    Kmat, inact)
{
  nbet    <- length(beta)
  Dydim   <- dim(Dyhat)
  nobs    <- Dydim[1]
  nbasis  <- Dydim[2]
  temp    <- beta[nbet]*Dyhat
  temp    <- temp*(wtroot %*% matrix(1,1,nbasis))
  hessmat <- 2*crossprod(temp)/nobs
  #  adjust for penalty
  if (lambda > 0) hessmat <- hessmat + 2*Kmat
  #  adjust for inactive coefficients
  if (any(inact)) {
    eyemat               <- diag(rep(1,nbasis))
    hessmat[inact,     ] <- 0
    hessmat[     ,inact] <- 0
    hessmat[inact,inact] <- eyemat[inact,inact]
  }
  return(hessmat)
}
