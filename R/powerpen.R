powerpen <- function(basisobj, Lfd=2) {
#  POWERPEN  Computes the power bais penalty matrix.
#  Arguments:
#  BASISFD  ... a monomial basis object
#  Lfd     ... either the order of derivative or a
#               linear differential operator to be penalized.
#  Returns a list the first element of which is the basis matrix
#   and the second element of which is the diagonal of the penalty matrix.

#  Last modified:  13 December 2002

	if (!(inherits(basisfd, "basis.fd"))) stop(
    	"First argument is not a basis.fd object.")

  	type <- getbasistype(basisobj)
  	rang <- basisobj$rangeval
  	if (type != "power") stop("BASISOBJ not of type POWER.")


  	exponents <- basisobj$params

  	if (!is.Lfd(Lfd))
    	stop (paste("Argument Lfd is neither a functional data object", 
             " nor an integer."))

	if (is.numeric(Lfd)) {
    	if (length(Lfd) == 1) {
      		nderiv <- Lfd
      		if (nderiv != as.integer(nderiv)) 
        		stop("Order of derivative must be an integer")
      		if (nderiv < 0) 
        		stop("Order of derivative must be 0 or positive")
    	} else {
      		stop("Order of derivative must be a single number")
    	}
    	if (nderiv < 0) 
			stop ("Order of derivative cannot be negative")

    	if (any(exponents - nderiv < 0) && rang[1] == 0)
        	stop("A negative exponent is needed and an argument value is 0.")
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
	return(penaltymat)
}

