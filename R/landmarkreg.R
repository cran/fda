landmarkreg <- function(fdobj, ximarks, x0marks=xmeanmarks,
                        WfdPar=fdPar(defbasis), monwrd=FALSE)
{
#  Arguments:
#  FDOBJ   ... functional data object for curves to be registered
#  XIMARKS ... N by NL array of times of interior landmarks for
#                 each observed curve
#  XOMARKS ... vector of length NL of times of interior landmarks for
#                 target curve
#  WFDPAR  ... a functional parameter object defining a warping function
#  MONWRD  ... If TRUE, warping functions are estimated by monotone smoothing,
#                 otherwise by regular smoothing.  The latter is faster, but
#                 not guaranteed to produce a strictly monotone warping
#                 function.  If MONWRD is 0 and an error message results
#                 indicating nonmonotonicity, rerun with MONWRD = 1.
#                 Default:  TRUE
#  Returns:
#  FDREG   ... a functional data object for the registered curves
#  WARPFD  ... a functional data object for the warping functions

 #  Last modified 26 October 2005

  #  check FDOBJ

  if (!(inherits(fdobj,  "fd"))) stop(
		"Argument fdobj  not a functional data object.")

  coef  <- fdobj$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)
  nrep  <- coefd[2]

  basisobj <- fdobj$basis
  type     <- basisobj$type
  nbasis   <- basisobj$nbasis
  rangeval <- basisobj$rangeval
  fdParobj <- fdPar(basisobj, 2, 1e-10)

  #  set up a default basis

  defbasis <- create.bspline.basis(rangeval,5)

  #  check WFDPAR

  if (inherits(WfdPar, "basisfd") & inherits(WfdPar, "fd"))
		WfdPar <- fdPar(WfdPar)
		
  if (!inherits(WfdPar, "fdPar")) stop(
		"WFDPAR is not a fdPar object.")
		
  Wfd0   <- WfdPar$fd
  wLfd   <- WfdPar$Lfd
  wbasis <- Wfd0$basis
  wrange <- wbasis$rangeval
  if (any(rangeval != wrange)) stop(
		"Ranges for FD and WFDPAR do not match.")

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

  y   <- eval.fd(x, fdobj)
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
       Wfd    <- smooth.morph(xval, yval, WfdPar, wval)$Wfdobj
       h      <- monfn(x, Wfd, 1)
       warpfd <- data2fd(h, x, wbasis)
    } else {
       warpfd <- smooth.basis(xval, yval, WfdPar, wval)$fd
       #  set up warping function by evaluating at sampling values
       h <- as.vector(eval.fd(x, warpfd))
       h <- h*(rangeval[2]-rangeval[1])/(h[n]-h[1])
       h <- h - h[1] + rangeval[1]
       #  check for monotonicity
       deltah <- diff(h)
       if (any(deltah <= 0)) warning(
           paste("Non-increasing warping function estimated for curve",irep))
    }
    hfunmat[,irep] <- h
    #  compute h-inverse
    if (monwrd) {
       wcoef  <- Wfd$coefs
       Wfdinv <- fd(-wcoef,wbasis)
       WfdParinv <- fdPar(Wfdinv, wLfd, lambda)
       Wfdinv <- smooth.morph(h, x, WfdParinv, wtn)$Wfdobj
       hinv   <- monfn(x, Wfdinv, 1)
    } else {
       hinvfd <- smooth.basis(h, x, WfdPar)$fd
       hinv   <- as.vector(eval.fd(x, hinvfd))
       hinv[1] <- x[1]
       hinv[n] <- x[n]
       deltahinv <- diff(hinv)
       if (any(deltahinv <= 0)) warning(
           paste("Non-increasing warping function estimated for curve",irep))
    }

    #  compute registered curves

    if (length(dim(coef)) == 2) {
      #  single variable case
      yregfd <- smooth.basis(hinv, y[,irep], fdParobj, wtn)$fd
      yregmat[,irep] <- eval.fd(x, yregfd)
    }
    if (length(dim(coef)) == 3) {
      #  multiple variable case
      for (ivar in 1:nvar) {
        # evaluate curve as a function of h at sampling points
        yregfd <- smooth.basis(hinv, y[,irep,ivar], fdParobj, wtn)$fd
        yregmat[,irep,ivar] <- eval.fd(x, yregfd)
       }
    }
  }

  #  create functional data objects for the registered curves

  yregcoef    <- project.basis(yregmat, x, basisobj)
  fdregnames  <- fdobj$fdnames
  names(fdregnames)[3] <- paste("Registered",names(fdregnames)[3])
  regfdobj    <- fd(yregcoef, basisobj, fdregnames)

  #  create functional data objects for the warping functions

  warpcoef    <- project.basis(hfunmat, x, wbasis)
  warpfdnames <- fdobj$fdnames
  names(warpfdnames)[3] <- paste("Warped",names(fdregnames)[1])
  warpfdobj   <- fd(warpcoef, wbasis, warpfdnames)

  return( list("regfd" = regfdobj, "warpfd" = warpfdobj) )
}
