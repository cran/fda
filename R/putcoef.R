putcoef <- function(coef, fd) {
	#  replace the coefficient matrix in a functional data object
	  
	#  Last modified 10 Feb 2003
  
  	if (inherits(fd, "fd")) {
		fd$coefs <- as.array(coef)
  	}  else {
      	stop("Argument FD is not a functional data object.")
  	}
  	return(fd)
}
