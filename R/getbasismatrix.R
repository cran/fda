getbasismatrix <- function(evalarg, basisfd, Lfd=0) {
#  Computes the basis matrix evaluated at arguments in EVALARG associated
#    with basis.fd object BASISFD.  The basis matrix contains the values   
#    at argument value vector EVALARG of applying the nonhomogeneous 
#    linear differential operator LFD to the basis functions.  By default
#    LFD is 0, and the basis functions are simply evaluated at argument
#    values in EVALARG.
#
#  If LFD is a functional data object with m + 1 functions c_1, ... c_{m+1}, then it
#    is assumed to define the order m NONHOMOGENEOUS linear differential operator
#  Lx(t) = c_1(t) + c_2(t)x(t) + c_3(t)Dx(t) + ... + c_{m+1}D^{m-1}x(t) + D^m x(t).
#  This is a change from previous usage where LFD was assumed to define a HOMOGONEOUS
#  differential operator, for which the forcing function c_1(t) = 0.  
#
#  If the basis type is either polygonal or constant, LFD is ignored.
#
#  Arguments:
#  EVALARG ... Either a vector of values at which all functions are to evaluated,
#              or a matrix of values, with number of columns corresponding to
#              number of functions in argument FD.  If the number of evaluation
#              values varies from curve to curve, pad out unwanted positions in
#              each column with NA.  The number of rows is equal to the maximum
#              of number of evaluation points.
#  BASISFD ... A basis object
#  LFD     ... If an integer, defines NDERIV, the order of derivative to be evaluated
#              If a functional data object, defines weight
#              functions for computing the value of a nonhomogeneous linear 
#              differential operator applied to the functions that are evaluated.
#
#  Note that the first two arguments may be interchanged.
#
#  Last modified 13 December 2002
 
#  Exchange the first two arguments if the first is an BASIS.FD object 
#    and the second numeric

if (is.numeric(basisfd) && inherits(evalarg, "basis.fd")) {
    temp    <- basisfd
    basisfd <- evalarg
    evalarg <- temp
}

if (!(inherits(basisfd, "basis.fd"))) stop(
    "Second argument is not a basis object.")

type   <- getbasistype(basisfd)
nbasis <- basisfd$nbasis

#  determine the highest order of derivative NDERIV required

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
   	Lfd <- NULL
   	if (nderiv < 0) stop ("Order of derivative cannot be negative")
} else if (inherits(Lfd, "fd")) {
   	derivcoef <- getcoef(Lfd)
   	if (length(dim(derivcoef))==3) derivcoef <- derivcoef[,,1]
   	nderiv <- dim(derivcoef)[2] - 1
   	if (nderiv < 0) {
   		stop("Order of derivative must be 0 or positive")
   	}
} else {
   	stop("Argument LFD must be an integer or a functional data object")
}

onerow <- rep(1,nbasis)

#  -------------------------------  Fourier basis  -------------------

if        (type == "fourier") {
   	period <- basisfd$params[1]
   	basis  <- fourier(evalarg, nbasis, period, nderiv)
   	if (nderiv > 0 && !is.null(Lfd)) {
        Lfdmat <- eval.fd(evalarg, Lfd)
        if (length(dim(Lfdmat)) == 3) Lfdmat <- Lfdmat[,,1]
        for (j in 1:nderiv) {
            if (any(abs(Lfdmat[,j+1])) > 1e-7) basis <- 
                basis + outer(Lfdmat[,j+1],onerow)*
                         fourier(evalarg, nbasis, period, j-1)
        }
    }

#  -----------------------------  B-spline basis  -------------------

} else if (type == "bspline") {
   	rangex <- basisfd$rangeval
   	breaks <- c(rangex[1], basisfd$params, rangex[2])
   	norder <- basisfd$nbasis - length(breaks) + 2
   	basis  <- bsplineS(evalarg, breaks, norder, nderiv)
   	if (nderiv > 0 && !is.null(Lfd)) {
        Lfdmat <- eval.fd(evalarg, Lfd)
        if (length(dim(Lfdmat)) == 3) Lfdmat <- Lfdmat[,,1]
        for (j in 1:nderiv) {
            if (any(abs(Lfdmat[,j+1])) > 1e-7) basis <- 
                basis + outer(Lfdmat[,j+1],onerow)*
                         bsplineS(evalarg, breaks, norder, j-1)
        }
    }

#  -----------------------------  Polynomial basis  -------------------

} else if (type == "poly") {
   	norder <- basisfd$nbasis
   	ctr    <- basisfd$params[1]
   	basis  <- polynom(evalarg, norder, nderiv, ctr)
   	if (nderiv > 0 && !is.null(Lfd)) {
        Lfdmat <- eval.fd(evalarg, Lfd)
        if (length(dim(Lfdmat)) == 3) Lfdmat <- Lfdmat[,,1]
        for (j in 1:nderiv) {
            if (any(abs(Lfdmat[,j+1])) > 1e-7) basis <- 
                basis + outer(Lfdmat[,j+1],onerow)*
                         polynom(evalarg, norder, j-1, ctr)
        }
    }

#  -----------------------------  Exponential basis  -------------------

} else if (type == "expon") {
   	basis  <- expon(evalarg, basisfd$params, nderiv)
   	if (nderiv > 0 && !is.null(Lfd)) {
        Lfdmat <- eval.fd(evalarg, Lfd)
        if (length(dim(Lfdmat)) == 3) Lfdmat <- Lfdmat[,,1]
        for (j in 1:nderiv) {
            if (any(abs(Lfdmat[,j+1])) > 1e-7) basis <- 
                basis + outer(Lfdmat[,j+1],onerow)*
                         expon(evalarg, basisfd$params, j-1)
        }
    }

#  -----------------------------  Polygonal basis  -------------------

} else if (type == "polyg") {
    basis  <- polyg(evalarg, basisfd$params)

#  -----------------------------  Power basis  -------------------

} else if (type == "power") {
    basis  <- powerbasis(evalarg, basisfd$params, nderiv)
   	if (nderiv > 0 && !is.null(Lfd)) {
        Lfdmat <- eval.fd(evalarg, Lfd)
        if (length(dim(Lfdmat)) == 3) Lfdmat <- Lfdmat[,,1]
        for (j in 1:nderiv) {
            if (any(abs(Lfdmat[,j+1])) > 1e-7) basis <- 
                basis + outer(Lfdmat[,j+1],onerow)*
                         powerbasis(evalarg, basisfd$params, j-1)
        }
    }

#  -----------------------------  Constant basis  --------------------

} else if (type == "const") {
   	basis  <- rep(1,length(evalarg))

} else {
   	stop("Basis type not recognizable")
}

if(length(evalarg) == 1) basis <- matrix(basis,1,nbasis)
return(basis)

}
