warpsmth <- function(x, y, wt=rep(1,nobs), Wfd,
                            Lfd=1, lambda=0, conv=.0001, iterlim=20,
                            active=rep(TRUE,nbasis), dbglev=0) {
#WARPSMTH smooths the relationship of Y to X using weights in WT by
#  fitting a monotone function of the form
#                   f(x) = b_0 + b_1 D^{-1} exp W(x)
#     where  W  is a function defined over the same range as X,
#                 W + ln b_1 = log Df and w = D W = D^2f/Df.
#  b_0 and b_1 are chosen so that f(x_1) = y_1 and f(x_n) = y_n.
#  The fitting criterion is penalized mean squared error:
#    PENSSE(lambda) = \sum [y_i - f(x_i)]^2 +
#                     \lambda * \int [L W(x)]^2 dx
#  where L is a linear differential operator defined in argument LFD.
#  The function W(x) is expanded by the basis in functional data object
#    WFD.

#  Arguments:
#  X       ...  vector of argument values
#  Y       ...  vector of function values to be fit
#  WT      ...  a vector of weights
#  WFD     ...  functional data object for W.  It's coefficient array
#               has a single column, and these are the starting values
#               for the iterative minimization of mean squared error.
#  LFD     ...  linear differential opr defining roughness penalty to
#               be applied to W.  This may be either a functional data
#               object defining a linear differential operator, or a
#               nonnegative integer.  If the latter, it specifies the
#               order of derivative to be penalized.
#               LFD = 1 by default, corresponding to L = D.
#  LAMBDA  ...  smoothing parameter determining the amount of penalty,
#               0 by default.
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
#  WFD       ...  functional data object for W.  Its coefficients are
#                   those that optimize fit.
#  BETA      ...  final regression coefficient values
#  FNEW      ...  final function value
#  MSG       ...  final gradient norm
#  ITERNUM   ...  number of iterations
#  ITERHIST  ...  ITERNUM+1 by 5 array containing iteration history

#  Last modified 4 July 2001
  if (!(inherits(Wfd, "fd"))) stop('Argument WFD is not a functional data object.')

#  initialize some arrays

  nobs   <- length(x)      #  number of observations
  basis  <- getbasis(Wfd)  #  basis for W
  nbasis <- basis$nbasis   #  number of basis functions

#  the starting values for the coefficients are in FD object WFD

  cvec   <- getcoef(Wfd)

#  check some arguments

  if (any(wt < 0))  stop("One or more weights are negative.")
  if (all(wt == 0)) stop("All weights are zero.")

#  set up some variables

  wtroot <- sqrt(wt)
  wtrtmt <- outer(wtroot,rep(1,ncovp1))
  yroot  <- y*wtroot
  climit <- matrix(100,nbasis,2) %*% matrix(c(-1,0,0,1),2,2)
  inact  <- !active   #  indices of inactive coefficients

#  initialize matrix Kmat defining penalty term

  if (lambda > 0) {
    Kmat <- lambda*getbasispenalty(basis, Lfd)
  } else {
    Kmat <- matrix(0,nbasis,nbasis)
  }

#  Compute initial function and gradient values

  result <- warpfngrad(y, x, wt, Wfd, lambda, Kmat, inact)
  Flist  <- result[[1]]
  Dyhat  <- result[[2]]

#  compute the initial expected Hessian

  hessmat <- warphesscal(Dyhat, wtroot, lambda, Kmat, inact)

#  evaluate the initial update vector for correcting the initial cvec

  result   <- linesearch(Flist, hessmat, dbglev)
  deltac   <- result[[1]]
  cosangle <- result[[2]]

#  initialize iteration status arrays

  iternum <- 0
  status  <- c(iternum, Flist$f, Flist$norm)
  if (dbglev >= 1) {
    cat ("Iter.   PENSSE   Grad Length Intercept\n")
    cat(paste(status[1],round(status[2],4),round(status[3],4),'\n'))
  }

  iterhist <- matrix(0,iterlim+1,length(status))
  iterhist[1,]  <- status
  if (iterlim == 0)
     return( list( Wfd, Flist, iternum, iterhist ) )

#  -------  Begin iterations  -----------

  MAXSTEPITER <- 10
  MAXSTEP <- 100
  trial   <- 1
  reset   <- FALSE
  linemat <- matrix(0,3,5)
  cvecold <- cvec
  Foldlist <- Flist
  dbgwrd  <- dbglev >= 2
  for (iter in 1:iterlim)
  {
     iternum <- iternum + 1
     #  initialize logical variables controlling line search
     dblwrd <- c(0,0)
     limwrd <- c(0,0)
     stpwrd <- 0
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
      if (dbglev >= 2)
                cat(c("                  ",
                 stepiter,format(round(c(linemat[,1]),6))),"\n")
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
      Wfdnew  <- Wfd
      for (stepiter in 1:MAXSTEPITER)
      {
      #  ensure that step does not go beyond limits on parameters
        limflg  <- FALSE
        #  check the step size
        result <-
              stepchk(linemat[1,5], cvec, deltac, limwrd, ind,
                      climit[,1], climit[,2], dbgwrd)
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
        Wfdnew[[1]] <- as.matrix(cvecnew)
        result  <- warpfngrad(y, x, wt, Wfdnew, lambda, Kmat, inact)
        Flist   <- result[[1]]
        Dyhat   <- result[[2]]
        linemat[3,5] <- Flist$f
        #  compute new directional derivative
        linemat[2,5] <- sum(deltac*Flist$grad)
        if (dbglev >= 2)
              cat(paste("                  ",
                 stepiter,format(round(c(linemat[,5]),6))),"\n")
        #  compute next line search step, also test for convergence
        result  <- stepitnew(linemat, ips, ind, dblwrd, MAXSTEP, dbglev)
        linemat <- result[[1]]
        ips     <- result[[2]]
        ind     <- result[[3]]
        dblwrd  <- result[[4]]
        trial   <- linemat[1,5]
        #  ind == 0  mean convergence
        if (ind == 0 | ind == 5) break
     }
     #  end iteration loop
     cvec <- cvecnew
     Wfd  <- Wfdnew
     #  check that function value has not increased
     if (Flist$f > Foldlist$f) {
        # if it has, terminate iterations with a message
        if (dbglev >= 2) {
          cat("Criterion increased: ")
          cat(format(round(c(Foldlist$f, Flist$f),4)))
          cat("\n")
        }
        #  reset parameters and fit
        cvec     <- cvecold
        Wfd[[1]] <- cvec
        Flist    <- Foldlist
        deltac   <- -Flist$grad
        if (reset) {
          # This is the second time in a row that
          #  this has happened ... quit
          if (dbglev >= 2) cat("Reset twice, terminating.\n")
          return ( list( Wfd, Flist, iternum, iterhist) )
        } else {
          reset <- TRUE
        }
     } else {
       if (abs(Foldlist$f - Flist$f) < conv) break
       cvecold  <- cvec
       Foldlist <- Flist
       hessmat  <- warphesscal(Dyhat, wtroot, lambda, Kmat, inact)
       #  udate the line search direction
       result   <- linesearch(Flist, hessmat, dbglev)
       deltac   <- result[[1]]
       cosangle <- result[[2]]
       reset    <- FALSE
     }
     #  store iteration status
     status <- c(iternum, Flist$f, Flist$norm)
     iterhist[iter+1,] <- status
     if (dbglev >= 1)
        cat(paste(status[1],round(status[2],4),round(status[3],4),"\n"))
  }
  return ( list( Wfd, Flist, iternum, iterhist ) )
}

#  ----------------------------------------------------------------

warpfngrad <- function(y, x, wt, Wfd, lambda, Kmat, inact)
{
  nobs   <- length(x)
  cvec   <- getcoef(Wfd)
  nbasis <- length(cvec)
  h      <-   monfn(x, Wfd, TRUE);
  Dyhat  <- Dcmonfn(x, Wfd, TRUE);
  #  update residuals and function values
  res    <- y - h
  f      <- mean(res^2*wt)
  temp   <- Dyhat*outer(wt,rep(1,nbasis))
  grad   <- -2*crossprod(temp,res)/nobs
  if (lambda > 0) {
    grad <- grad +         2 * Kmat %*% cvec
    f    <- f    + t(cvec) %*% Kmat %*% cvec
  }
  if (any(inact)) grad[inact] <- 0
  norm <- sqrt(sum(Flist$grad^2)) #  gradient norm
  Flist <- list(f=f,grad=grad,norm=norm)
  return(list(Flist, Dyhat))
}

#  ----------------------------------------------------------------

warphesscal <- function(Dyhat, wtroot, lambda, Kmat, inact)
{
  Dydim   <- dim(Dyhat)
  nobs    <- Dydim[1]
  nbasis  <- Dydim[2]
  temp    <- Dyhat*outer(wtroot,rep(1,nbasis))
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

