derivs <- function(tnow, y, wfd) {
	#  Sets up a linear differential equation of order m 
	#  as a system of m first order equations
	#  Arguments:
	#  TNOW ... A vector of values of the independent variable t
	#  Y    ... A matrix of m values of Y corresponding to TNOW
	#  WFD  ... A functional data object containing coefficient functions
	#  Returns:
	#  DY   ... A matrix of derivative values corresponding to Y
	
	#  Last modified:  24 March 2003
	
  w  <- eval.fd(tnow,wfd)
  m  <- length(w) - 1;
  wmat <- matrix(0, m, m)
  wmat[1:(m-1),2:m] <- diag(rep(1,m-1))
  wmat[m,] <- -w[2:(m+1)]
  dy <- wmat %*% y
  return(dy)
}
