bsplinepen <- function(basisfd, Lfd=2)
{

#  Computes the Bspline penalty matrix.
#  Arguments:
#  BASISFD   ... a basis.fd object of type "bspline"
#  LFD ... either the order of derivative or a
#          nonhomogeneous linear differential operator to be penalized.
#  Returns the penalty matrix.

#  Last modified 5 December 2001

if (!(inherits(basisfd, "basis.fd"))) stop(
    "First argument is not a basis.fd object.")

type <- getbasistype(basisfd)
if (type != "bspline") stop("BASISFD not of type bspline")

norder <- basisfd$nbasis - length(basisfd$params)

#  Find the highest order derivative in LFD

if (is.numeric(Lfd)) {
    if (length(Lfd) == 1) {
      	nderiv <- Lfd
      	if (nderiv != as.integer(nderiv)) {
        	stop("Order of derivative must be an integer")
      	}
      	if (nderiv < 0) {
        	stop("Order of derivative must be 0 or positive")
      	}
    } else {
      	stop("Order of derivative must be a single number")
    }
    if (nderiv < 0) stop ("Order of derivative cannot be negative")
} else if (inherits(Lfd, "fd")) {
   	derivcoef <- getcoef(Lfd)
   	if (length(dim(derivcoef))==3) derivcoef <- derivcoef[,,1]
   	nderiv <- dim(derivcoef)[2] - 1
   	if (nderiv < 0) {
   		stop("Order of derivative must be 0 or positive")
   	}
    nderiv <- ncol(derivcoef)
} else {
    stop("Second argument must be an integer or a functional data object")
}

if (nderiv >= norder) {
	cat(paste("Derivative of order", nderiv,
                  "cannot be taken for B-spline of order", norder,"\n"))
	cat("Probable cause is a value of the nbasis argument\n")
	cat(" in function create.basis.fd that is too small.\n")
   	stop()
}

if (nderiv == norder - 1 && is.numeric(Lfd)) {
    breakvals  <- c(basisfd$rangeval[1], basisfd$params,
                    basisfd$rangeval[2])
    nbreakvals <- length(breakvals)
    norder     <- basisfd$nbasis - nbreakvals + 2
    if (nderiv >= norder) stop (
         "NDERIV cannot be as large as order of B-spline.")
    halfseq    <- (breakvals[2:nbreakvals] +
                   breakvals[1:(nbreakvals-1)])/2
    halfmat    <- bsplineS(halfseq, breakvals, norder, nderiv)
    brwidth    <- diff(breakvals)
    penaltymatrix <- t(halfmat) %*% diag(brwidth) %*% halfmat
} else {
    penaltymatrix <- inprod(basisfd, basisfd, Lfd, Lfd)
}

return( penaltymatrix )
}
