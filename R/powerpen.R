powerpen <- function(basisobj, Lfdobj=int2Lfd(2)) {
#  POWERPEN  Computes the power basis penalty matrix.
#  Arguments:
#  BASISFD  ... a monomial basis object
#  Lfd     ... either the order of derivative or a
#               linear differential operator to be penalized.
#  Returns a list the first element of which is the basis matrix
#   and the second element of which is the diagonal of the penalty matrix.

#  Last modified:  26 October 2005

	if (!inherits(basisobj, "basisfd")) stop(
    	"First argument is not a basis object.")

  	type <- basisobj$type
  	rang <- basisobj$rangeval
  	if (type != "power") stop("BASISOBJ not of type POWER.")


  	exponents <- basisobj$params

  Lfdobj <- int2Lfd(Lfdobj)

  if (is.integerLfd(Lfdobj)) {
      nderiv <- Lfdobj$nderiv

    	if (any(exponents - nderiv < 0) && rang[1] == 0)
        	stop(paste("A negative exponent is needed and",
                    "an argument value is 0."))
    	nbasis     <- basisobj$nbasis
    	penaltymat <- matrix(0,nbasis,nbasis)
    	xrange     <- basisobj$rangeval
    	for (ibasis in 1:nbasis) {
      		ideg <- exponents[ibasis]
      		if (nderiv == 0) ifac <- 1 else ifac <- ideg
			if (nderiv > 1)
				for (k in 2:nderiv) ifac <- ifac*(ideg - k + 1)
      		for (jbasis in 1:ibasis) {
        		jdeg <- exponents[jbasis]
        		if (nderiv == 0) jfac <- 1 else jfac <- jdeg
				if (nderiv > 1)
	  				for (k in 2:nderiv) jfac <- jfac*(jdeg - k + 1)
				if (ideg >= nderiv && jdeg >= nderiv) {
	  				penaltymat[ibasis,jbasis] <- ifac*jfac*
	      				(xrange[2]^(ideg+jdeg-2*nderiv+1) -
	       		 	 xrange[1]^(ideg+jdeg-2*nderiv+1))
	  				penaltymat[jbasis,ibasis] <- penaltymat[ibasis,jbasis]
				}
      		}
    	}
  	} else {
    	penaltymat <- inprod(basisobj, basisobj, Lfd, Lfd)
	}
	penaltymat
}

