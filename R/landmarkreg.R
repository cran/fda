landmarkreg <- function(fd, fd0, ximarks, x0marks=xmeanmarks, wbasis = basis,
                        Lfd=0, sparval=1e-10, monwrd=FALSE)
{
#  Arguments:
#  FD      ... functional data object for curves to be registered
#  FD0     ... functional data object for the target curve
#  XIMARKS ... N by NL array of times of interior landmarks for
#                 each observed curve
#  XOMARKS ... vector of length NL of times of interior landmarks for
#                 target curve
#  WBASIS  ... optional basis object used for estimating warp
#                 functions.  If not supplied the basis for FDOBJ is used.
#  LFD     ... integer or functional data object defining derivative
#                 or LDO value to be registered.
#  SPARVAL ... smoothing parameter used by smooth.spline
#  MONWRD  ... If T, warping functions are estimated by monotone smoothing,
#                 otherwise by regular smoothing.  The latter is faster, but
#                 not guaranteed to produce a strictly monotone warping
#                 function.  If MONWRD is 0 and an error message results
#                 indicating nonmonotonicity, rerun with MONWRD = 1.
#                 Default:  T
#  Returns:
#  FDREG   ... a functional data object for the registered curves
#  WARPFD  ... a functional data object for the warping functions


 #  Last modified 6 Feb 2001

  if (!(inherits(fd,  "fd"))) stop("Argument FD  not a functional data object.")
  if (!(inherits(fd0, "fd"))) stop("Argument FD0 not a functional data object.")

  coef  <- getcoef(fd)
  coefd <- dim(coef)
  ndim  <- length(coefd)
  nrep  <- coefd[2]

  basis    <- getbasis(fd)
  type     <- basis$type
  nbasis   <- basis$nbasis
  rangeval <- basis$rangeval

  ximarksd <- dim(ximarks)
  if (ximarksd[1] != nrep) stop(
     "Number of rows of third argument wrong.")
  nlandm <- dim(ximarks)[2]
  xmeanmarks <- apply(ximarks,2,mean)
  if (length(x0marks) != nlandm) stop(
     "Number of target landmarks not equal to number of curve landmarks.")

  if (any(ximarks <= rangeval[1]) || any(ximarks >= rangeval[2])) stop(
     "Some landmark values are not within the range.")

  n   <- min(c(101,10*nbasis))
  x   <- seq(rangeval[1],rangeval[2],length=n)
  wtn <- rep(1,n)

  y0  <- eval.fd(x, fd0, Lfd)
  y   <- eval.fd(x, fd,  Lfd)
  yregmat <- y
  hfunmat <- matrix(0,n,nrep)

  xval <- c(rangeval[1],x0marks,rangeval[2])
  nval <- length(xval)
  wval <- rep(1,nval)

  cat("Progress:  Each dot is a curve\n")
  for (irep in 1:nrep) {
    cat(".")
    #  set up landmark times for this curve
    yval   <- c(rangeval[1],ximarks[irep,],rangeval[2])
    #  smooth relation between this curve"s values and target"s values
    if (monwrd) {
       #  use monotone smoother
       result <- warpsmth(xval, yval, wval, Wfd0, wLfd, lambda)
       Wfd    <- result[[1]]
       h      <- monfn(x, Wfd, 1)
       warpfd <- data2fd(h, x, wbasis)
    } else {
       warpfd <- smooth.basis(yval, xval, wbasis, wval, 2, sparval)$fd
       #  set up warping function by evaluating at sampling values
       h <- eval.fd(x, warpfd)
       h <- h*(rangeval[2]-rangeval[1])/(h[n]-h[1])
       h <- h - h[1] + rangeval[1]
       #  check for monotonicity
       deltah <- diff(h)
       if (any(deltah) <= 0) warning(
           paste("Non-increasing warping function estimated for curve",irep))
    }
    hfunmat[,irep] <- h
    #  compute h-inverse
    if (monwrd) {
       wcoef  <- getcoef(Wfd)
       Wfdinv <- fd(-wcoef,wbasis)
       result <- warpsmth(h, x, wtn, Wfdinv, wLfd, lambda)
       Wfdinv <- result[[1]]
       hinv   <- monfn(x, Wfdinv, 1)
    } else {
       hinvfd <- smooth.basis(x, h, wbasis, wtn, 2, sparval)$fd
       hinv   <- eval.fd(x, hinvfd)
       hinv[1] <- x[1]
       hinv[n] <- x[n]
       deltahinv <- diff(hinv)
       if (any(deltahinv) <= 0) warning(
           paste("Non-increasing warping function estimated for curve",irep))
    }

    #  compute registered curves

    if (length(dim(coef)) == 2) {
      #  single variable case
      yregfd <- smooth.basis(y[,irep], hinv, basis, wtn, 2, 1e-10)$fd
      yregmat[,irep] <- eval.fd(x, yregfd)
    }
    if (length(dim(coef)) == 3) {
      #  multiple variable case
      for (ivar in 1:nvar) {
        # evaluate curve as a function of h at sampling points
        yregfd <- smooth.basis(y[,irep,ivar], hinv, basis, wtn, 2,
                      1e-10)$fd
        yregmat[,irep,ivar] <- eval.fd(x, yregfd)
       }
    }
  }

  #  create functional data objects for the registered curves

  yregcoef    <- project.basis(yregmat, x, basis)
  fdregnames  <- getnames(fd)
  names(fdregnames)[3] <- paste("Registered",names(fdregnames)[3])
  regfd       <- create.fd(yregcoef, basis, fdregnames)

  #  create functional data objects for the warping functions

  warpcoef    <- project.basis(hfunmat, x, wbasis)
  warpfdnames <- fdregnames
  names(warpfdnames)[3] <- paste("Warped",names(fdregnames)[1])
  warpfd      <- create.fd(warpcoef, wbasis, warpfdnames)

  return( list("regfd" = regfd, "warpfd" = warpfd, x) )
}
