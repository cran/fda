exponpen <- function(basisfd, Lfd=2)
{

#  Computes the Exponential penalty matrix.
#  Argument:
#  BASISFD ... a basis.fd object of type "expon"
#  LFD     ... either the order of derivative or a
#                linear differential operator to be penalized.
#  Returns the penalty matrix.

#  Last modified 5 December 2001

if (!(inherits(basisfd, "basis.fd"))) stop(
    "First argument is not a basis.fd object.")

type <- getbasistype(basisfd)
if (type != "expon") stop ("Wrong basis type")

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
    ratevec <- basisfd$params
    nrate   <- length(ratevec)
    penaltymatrix <- matrix(0,nrate,nrate)
    tl <- basisfd$rangeval[1]
    tu <- basisfd$rangeval[2]
    for (irate in 1:nrate) {
      	ratei <- ratevec[irate]
      	for (jrate in 1:irate) {
        	ratej <- ratevec[jrate]
        	ratesum <- ratei + ratej
        	if (ratesum != 0) {
          		penaltymatrix[irate,jrate] <- (ratei*ratej)^nderiv *
              	(exp(ratesum*tu) - exp(ratesum*tl)) / ratesum
        	} else {
          		if (nderiv == 0) penaltymatrix[irate,jrate] <- tu - tl
        	}
        	penaltymatrix[jrate,irate] <- penaltymatrix[irate,jrate]
      	}
	}
} else {
    penaltymatrix <- inprod(basisfd, basisfd, Lfd, Lfd)
}

return( penaltymatrix )
}
