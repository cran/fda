getcoef <- function(fd) {
  	#  Extracts the coefficient array from a functional data object FD, or,
  	#    if FD is a basis object, assigns the identity matrix as the coefficient
  	#    array

  	#  Last modified 10 Feb 2003
  
  	if (inherits(fd, "fd")) {
    	coef <- as.array(fd$coefs)
    	if (length(dim(coef)) == 1) coef <- matrix(coef, length(coef), 1)
  	}  else {
    	if (inherits(fd, "basis.fd")) {
      		coef <- diag(rep(1,fd$nbasis))
    	} else {
      		stop("An object of class fd or basis.fd expected, but not found.")
    	}
  	}
  	return(coef)
}
