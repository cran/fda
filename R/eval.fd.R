#  ----------------------------------------------------------------------------

eval.fd <- function(evalarg, fd, Lfd=0) {

#  Evaluates a functional data observation, or the value of a linear
#  differential operator LFD applied to the object, at argument values in an array 
#  EVALARGS.
#
#  Note that this function is preferred to eval.fd() since it is so easy to confuse
#  eval.fd() with the function eval(), which is generic and will process the argument 
#  list for eval.fd, but give inappropriate results.
#
#  If LFD is a functional data object with m + 1 functions c_1, ... c_{m+1}, then it
#    is assumed to define the order m NONHOMOGENEOUS linear differential operator
#  Lx(t) = c_1(t) + c_2(t)x(t) + c_3(t)Dx(t) + ... + c_{m+1}D^{m-1}x(t) + D^m x(t).
#  This is a change from previous usage where LFD was assumed to define a HOMOGONEOUS
#  differential operator, for which the forcing function c_1(t) = 0.  
#
#  Arguments:
#  EVALARG ... Either a vector of values at which all functions are to evaluated,
#              or a matrix of values, with number of columns corresponding to
#              number of functions in argument FD.  If the number of evaluation
#              values varies from curve to curve, pad out unwanted positions in
#              each column with NA.  The number of rows is equal to the maximum
#              of number of evaluation points.
#  FD      ... Functional data object
#  LFD     ... If an integer, defines NDERIV, the order of derivative to be evaluated
#              If a functional data object, defines weight
#              functions for computing the value of a nonhomogeneous linear 
#              differential operator applied to the functions that are evaluated.
#  Note that the first two arguments may be interchanged.

#  Returns:  An array of function values corresponding to the evaluation 
#              arguments in EVALARG

#  Last modified 5 December 2001

#  Exchange the first two arguments if the first is an BASIS.FD object 
#    and the second numeric

if (is.numeric(fd) && inherits(evalarg, "fd")) {
    temp    <- fd
    fd      <- evalarg
    evalarg <- temp
}

#  Check the arguments

evaldim <- dim(evalarg)

if (!(is.numeric(evalarg)))  stop("Argument EVALARG is not numeric.")
	
if (!(length(evaldim) < 3)) stop(
   "Argument EVALARG is not a vector or a matrix.")

if (!(inherits(fd, "fd"))) stop(
     "Argument FD is not a functional data object.")

if (!(is.Lfd(Lfd))) stop(
     "Argument LFD is not a linear differential operator.")

#  Extract information about the basis

basisfd  <- getbasis(fd)
nbasis   <- basisfd$nbasis
rangeval <- basisfd$rangeval
onerow   <- rep(1,nbasis)

temp <- c(evalarg)
temp <- temp[!(is.na(temp))]
if (min(temp) < rangeval[1] || max(temp) > rangeval[2]) 
	warning(paste(
    "Values in argument EVALARG are outside of permitted range,",
    "and will be ignored."))

#  get maximum number of evaluation values

if (is.vector(evalarg)) {
	n <- length(evalarg)
} else {
	n <- evaldim[1]
}

#  Set up coefficient array for FD

coef  <- getcoef(fd)
coefd <- dim(coef)
ndim  <- length(coefd)
if (ndim <= 1) nrep <- 1 else nrep <- coefd[2]
if (ndim <= 2) nvar <- 1 else nvar <- coefd[3]

#  Set up array for function values

if (ndim <= 2) evalarray <- matrix(0,n,nrep)
	else        evalarray <- array(0,c(n,nrep,nvar))
if (ndim == 2) dimnames(evalarray) <- list(NULL,dimnames(coef)[[2]])
if (ndim == 3) dimnames(evalarray) <- list(NULL,dimnames(coef)[[2]],
	                                             dimnames(coef)[[3]])

#  Case where EVALARG is a vector of values to be used for all curves

if (is.vector(evalarg)) {

    evalarg[evalarg < rangeval[1]-1e-10] <- NA
    evalarg[evalarg > rangeval[2]+1e-10] <- NA
    basismat <- getbasismatrix(evalarg, basisfd, Lfd)

    #  evaluate the functions at arguments in EVALARG

   	if (inherits(Lfd, "fd")) {
		Lfdmat <- eval.fd(evalarg, Lfd)
		if (length(dim(Lfdmat))==3) Lfdmat <- Lfdmat[,,1]
		force <- outer(Lfdmat[,1],rep(1,nrep))
	} else {
		force <- outer(rep(0,length(evalarg)),rep(1,nrep))
	}
    if (ndim <= 2) {
	    evalarray <- force + basismat %*% coef
    } else {
       evalarray <- array(0,c(n,nrep,nvar))
       for (ivar in 1:nvar) evalarray[,,ivar] <- force + basismat %*% coef[,,ivar]
    }

} else {
	
	#  case of evaluation values varying from curve to curve
	
	for (i in 1:nrep) {
		evalargi <- evalarg[,i]
       if (all(is.na(evalargi))) stop(
            paste("All values are NA for replication",i))

		index    <- !(is.na(evalargi) | evalargi < rangeval[1] | evalargi > rangeval[2])
		evalargi <- evalargi[index]
       basismat <- getbasismatrix(evalargi, basisfd, Lfd)

   	  	if (inherits(Lfd, "fd")) {
			Lfdmat <- eval.fd(evalargi, Lfd)
			if (length(dim(Lfdmat))==3) Lfdmat <- Lfdmat[,,1]
			force <- outer(Lfdmat[,1],rep(1,nrep))
		} else {
		   	force <- outer(rep(0,length(evalargi)),rep(1,nrep))
		}
       #  evaluate the functions at arguments in EVALARG

       if (ndim == 2) {
           evalarray[  index, i] <- force + basismat %*% coef[,i]
           evalarray[!(index),i] <- NA
       }
       if (ndim == 3) {
           for (ivar in 1:nvar) {
	           evalarray[   index,i,ivar] <- force + basismat %*% coef[,i,ivar]
               evalarray[!(index),i,ivar] <- NA
           }
       }
	}
	
}

return(evalarray)

}

