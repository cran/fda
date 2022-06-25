landmarkreg <- function(unregfd, ximarks, x0marks, x0lim=NULL, 
                        WfdPar=NULL, WfdPar0=NULL, ylambda=1e-10) {
  #  This version of landmarkreg does not assume that the target marks
  #  x0marks are within the same interval as that for the curves to be
  #  registered.  Consequently, it needs a required extra argument X0LIM 
  #  defining the target interval and optional fdPar argument for 
  #  representing the inverse warping function.
  
  #  Arguments:
  #  UNREGFD ... functional data object for curves to be registered
  #  XIMARKS ... N by NL array of times of interior landmarks for
  #                 each observed curve
  #  XOMARKS ... vector of length NL of times of interior landmarks for
  #                 target curve
  #  X0LIM   ... vector of length 2 containing the lower and upper boundary
  #              of the target interval containing x0marks.
  #  WFDPAR  ... a functional parameter object defining a warping function
  #  MONWRD  ... If TRUE, warping functions are estimated by monotone smoothing,
  #                 otherwise by regular smoothing.  The latter is faster, but
  #                 not guaranteed to produce a strictly monotone warping
  #                 function.  If MONWRD is 0 and an error message results
  #                 indicating nonmonotonicity, rerun with MONWRD = 1.
  #                 Default:  TRUE
  #  YLAMBDA ... smoothing parameter to be used in computing the registered
  #                 functions.  For high dimensional bases, local wiggles may be
  #                 found in the registered functions or its derivatives that are
  #                 not seen in the unregistered functions.  In this event, this
  #                 parameter should be increased.
  #  Returns:
  #  FDREG   ... a functional data object for the registered curves
  #  WARPFD  ... a functional data object for the warping functions
  #  WFD     ... a functional data object for the W functions defining the
  #              warping functions
  
  # Warning:  As of March 2022, landmark registration cannot be done using
  # function smooth.basis instead of function smooth.morph.  The 
  # warping function must be strictly monotonic, and we have found that using 
  # smooth.basis too often violates this contraint.  Function 
  # smooth.morph ensures monotonicity.
  
  #  Last modified 2 June 2022  by Jim Ramsay
  
  #  check unregfd containing the curves to be registered
  
  if (!(inherits(unregfd,  "fd"))) stop(
    "Argument unregfd  not a functional data object.")
  
  Ybasis   <- unregfd$basis
  nbasis   <- Ybasis$nbasis
  rangeval <- Ybasis$rangeval
  
  if (is.null(x0lim)) x0lim = rangeval
  
  #   ---------------------------------------------------------------------
  #                  check ximarks and x0marks
  #   ---------------------------------------------------------------------
  
  #  check ximarks being matrix with ncurve rows and nmarks columns
  
  if (is.numeric(ximarks)) {
    nximarks <- length(ximarks)
    # if ximarks is a vector, coerce it to a single row matrix
    if (is.vector(ximarks))     ximarks <- matrix(ximarks,1,nximarks)
    # if ximarks is a data.frame, coerce it to a matrix
    if (is.data.frame(ximarks)) ximarks <- as.matrix(ximarks)
  } else {
    stop("Argument ximarks is not numeric.")
  }
  
  #  check x0marks and coerce it to be a one-row matrix
  
  if (is.numeric(x0marks)) {
    nx0marks <- length(x0marks)
    if (is.vector(x0marks)) x0marks <- matrix(x0marks,1,nx0marks)
  } else {
    stop("Argument x0marks is not numeric.")
  }
  
  #  check that ximarks and x0marks have same number of columns
  
  if (ncol(ximarks) != length(x0marks)) 
    stop("The number of columns in ximarks is not equal to length of x0marks.")
  
  # check that ximarks are within range of unregfd
  
  if (any(ximarks <= rangeval[1]) || any(ximarks >= rangeval[2]))
    stop("Argument ximarks has values outside of range of unregfd.")
  
  # check that x0marks are within range of target interval
  
  if (any(x0marks <= x0lim[1]) || any(x0marks >= x0lim[2]))
    stop("Argument x0marks has values outside of range of target interval.")
  
  #  determine the number of curves to be registered
  
  ncurve   <- dim(ximarks)[1]
  
  #   ---------------------------------------------------------------------
  #                        check WFDPAR
  #   ---------------------------------------------------------------------
  
  #  set up default WfdPar for warping function
  
  if (is.null(WfdPar)) {
    Wnbasis   <- length(x0marks) + 2
    Ybasis    <- create.bspline.basis(rangeval, Wnbasis)
    Wfd       <- fd(matrix(0,Wnbasis,1), Wbasis)
    WfdPar    <- fdPar(Wfd, 2, 1e-10)
  } else {
    WfdPar  <- fdParcheck(WfdPar,  1)
    Wfd     <- WfdPar$fd
    Wbasis  <- Wfd$basis
    Wnbasis <- Wbasis$nbasis
  }

  #  set up default WfdPar0 for inverse warping function
  
  if (is.null(WfdPar0)) {
    Wnbasis0  <- length(x0marks) + 2
    Wbasis0   <- create.bspline.basis(x0lim, Wnbasis0)
    Wfd0      <- fd(matrix(0,Wnbasis0,1), Wbasis0)
    WfdPar0   <- fdPar(Wfd0, 2, 1e-10)
  } else {
    WfdPar0  <- fdParcheck(WfdPar0, 1)
    Wfd0     <- WfdPar0$fd
    Wbasis0  <- Wfd0$basis
    Wnbasis0 <- Wbasis0$nbasis
  }
  
  #   ---------------------------------------------------------------------
  #                        set up analysis
  #   ---------------------------------------------------------------------
  
  nfine   <- min(c(101,10*nbasis))
  xfine   <- seq(rangeval[1], rangeval[2], length=nfine)
  xfine0  <- seq(   x0lim[1],    x0lim[2], length=nfine)
  yfine   <- eval.fd(xfine, unregfd)
  yregmat <- yfine
  hfunmat <- matrix(0,nfine,ncurve)
  hinvmat <- matrix(0,nfine,ncurve)
  
  xval    <- matrix(c(x0lim[1],x0marks,x0lim[2]),nx0marks+2,1)
  Wcoef   <- matrix(0,Wnbasis,ncurve)
  nval    <- length(xval)
  
  #  --------------------------------------------------------------------
  #                  Iterate through curves to register
  #  --------------------------------------------------------------------
  
  if (ncurve > 1) cat("Progress:  Each dot is a curve\n")
  
  for (icurve in 1:ncurve) {
    if (ncurve > 1) cat(".")
    #  set up landmark times for this curve
    yval   <- matrix(c(rangeval[1],ximarks[icurve,],rangeval[2]),nx0marks+2,1)
    #  smooth relation between this curve"s values and target"s values
    #  use monotone smoother
    
    Wfd  <- smooth.morph(xval, yval, rangeval, WfdPar)$Wfdobj
    
    hfun <- monfn(xfine, Wfd)
    b    <- (rangeval[2]-rangeval[1])/(hfun[nfine]-hfun[1])
    a    <- rangeval[1] - b*hfun[1]
    hfun <- a + b*hfun
    hfun[c(1,nfine)] <- rangeval
    Wcoefi           <- Wfd$coef
    Wcoef[,icurve]   <- Wcoefi
    hfunmat[,icurve] <- hfun
    
    #  compute h-inverse  in order to register curves
    
    Wcoefi       <- Wfd$coefs
    Wfdinv       <- smooth.morph(hfun, xfine, x0lim, WfdPar0)$Wfdobj
    hinv         <- monfn(xfine, Wfdinv)
    b            <- (x0lim[2]-x0lim[1])/(hinv[nfine]-hinv[1])
    a            <- x0lim[1] - b*hinv[1]
    hinv         <- a + b*hinv
    hinv[c(1,nfine)] <- rangeval
    hinvmat[,icurve] <- hinv
    
    #  compute registered curves
    
    yregfd <- smooth.basis(hinv, yfine[,icurve], Ybasis)$fd
    yregmat[,icurve] <- eval.fd(xfine, yregfd, 0)
  }
  
  if (ncurve > 1) cat("\n")
  
  #  create functional data objects for the registered curves
  
  regfdPar <- fdPar(Ybasis, 2, ylambda)
  regfd    <- smooth.basis(xfine, yregmat, regfdPar)$fd
  regnames <- unregfd$fdnames
  names(regnames)[3] <- paste("Registered",names(regnames)[3])
  regfd$fdnames <- regnames
  
  #  create functional data objects for the warping functions
  
  warpfd                <- smooth.basis(xfine, hfunmat, Ybasis)$fd
  warpfdnames           <- unregfd$fdnames
  names(warpfdnames)[3] <- paste("Warped",names(regnames)[1])
  warpfd$fdnames        <- warpfdnames

  #  create functional data objects for the inverse warping functions
  
  Ybasis0               <- create.bspline.basis(x0lim, nbasis)
  warpinvfd             <- smooth.basis(xfine0, hinvmat, Ybasis0)$fd
  warpfdnames           <- unregfd$fdnames
  names(warpfdnames)[3] <- paste("Warped",names(regnames)[1])
  warpinvfd$fdnames     <- warpfdnames
  
  #  The core function defining the strictly monotone warping
  
  Wfd <- fd(Wcoef, Wbasis)
  
  return( list(regfd=regfd, warpfd=warpfd, warpinvfd=warpinvfd, Wfd=Wfd) )
}
