polynompen <- function(basisfd, Lfd=2)
{

#  Computes the polynomial penalty matrix for polynomials of the form
#      (x-ctr)^l
#  Arguments:
#  BASISFD   ... a basis.fd object of type "poly"
#  LFD ... either the order of derivative or a
#           nonhomogeneous linear differential operator to be penalized.
#  Returns the penalty matrix.

#  Last modified 5 December 2001

if (!(inherits(basisfd, "basis.fd"))) stop(
    "First argument is not a basis.fd object.")

type <- getbasistype(basisfd)
  if (type != "poly") stop("BASISFD not of type poly")

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

#  Compute penalty matrix

if (is.numeric(Lfd)) {
    nbasis   <- basisfd$nbasis
    rangex   <- basisfd$rangeval
    ctr      <- basisfd$params[1]
    basismat <- getbasismatrix(rangex, basisfd, nderiv)
    penmatl  <- outer(basismat[1,],basismat[1,])*(rangex[1] - ctr)
    penmatu  <- outer(basismat[2,],basismat[2,])*(rangex[2] - ctr)
    penaltymatrix <- matrix(0,nbasis,nbasis)
    for (i in (nderiv+1):nbasis) for (j in (nderiv+1):i) {
      	penaltymatrix[i,j] <- (penmatu[i,j] - penmatl[i,j])/(i + j - 2*nderiv - 1)
      	penaltymatrix[j,i] <- penaltymatrix[i,j]
    }
} else {
    penaltymatrix <- inprod(basisfd, basisfd, Lfd, Lfd)
}

  return( penaltymatrix )
}
