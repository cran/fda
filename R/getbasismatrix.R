getbasismatrix <- function(evalarg, basisobj, nderiv=0) {
#  Computes the basis matrix evaluated at arguments in EVALARG associated
#    with basis.fd object BASISOBJ.  The basis matrix contains the values
#    at argument value vector EVALARG of applying the nonhomogeneous
#    linear differential operator LFD to the basis functions.  By default
#    LFD is 0, and the basis functions are simply evaluated at argument
#    values in EVALARG.
#
#  If LFD is a functional data object with m + 1 functions c_1, ... c_{m+1}, then it
#    is assumed to define the order m HOMOGENEOUS linear differential operator
#  Lx(t) = c_1(t) + c_2(t)x(t) + c_3(t)Dx(t) + ... + c_{m+1}D^{m-1}x(t) + D^m x(t).
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
#  BASISOBJ ... A basis object
#  NDERIV   ... A nonnegative integer indicating a derivative to be evaluated.

#
#  Note that the first two arguments may be interchanged.
#
#  Last modified 26 October 2005

#  Exchange the first two arguments if the first is an BASIS.FD object
#    and the second numeric

if (is.numeric(basisobj) && inherits(evalarg, "basisfd")) {
    temp     <- basisobj
    basisobj <- evalarg
    evalarg  <- temp
}

#  check EVALARG

if (!(is.numeric(evalarg)))  stop("Argument EVALARG is not numeric.")
	
#  check basisobj
	
if (!(inherits(basisobj, "basisfd"))) stop(
    "Second argument is not a basis object.")

#  Extract information about the basis

type     <- basisobj$type
nbasis   <- basisobj$nbasis
params   <- basisobj$params
rangeval <- basisobj$rangeval
dropind  <- basisobj$dropind

#  -----------------------------  B-spline basis  -------------------

if (type == "bspline") {
   	breaks   <- c(rangeval[1], params, rangeval[2])
   	norder   <- nbasis - length(breaks) + 2
   	basismat <- bsplineS(evalarg, breaks, norder, nderiv)

#  -----------------------------  Constant basis  --------------------

} else if (type == "const") {
   	basismat  <- matrix(1,length(evalarg),1)

#  -----------------------------  Exponential basis  -------------------

} else if (type == "expon") {
   	basismat  <- expon(evalarg, params, nderiv)

#  -------------------------------  Fourier basis  -------------------

} else if (type == "fourier") {
   	period   <- params[1]
   	basismat <- fourier(evalarg, nbasis, period, nderiv)

#  -----------------------------  Monomial basis  -------------------

} else if (type == "monom") {
   	basismat  <- monomial(evalarg, params, nderiv)

#  -----------------------------  Polygonal basis  -------------------

} else if (type == "polyg") {
    basismat  <- polyg(evalarg, params)

#  -----------------------------  Polynomial basis  -------------------

} else if (type == "polynom") {
   	norder   <- nbasis
   	ctr      <- params[1]
   	basismat <- polynom(evalarg, norder, nderiv, ctr)

#  -----------------------------  Power basis  -------------------

} else if (type == "power") {
    basismat  <- powerbasis(evalarg, params, nderiv)

#  -----------------------  Unrecognizable basis  --------------------

} else {
   	stop("Basis type not recognizable")
}
	
#  remove columns for bases to be dropped

if (length(dropind) > 0) basismat <- basismat[,-dropind]

return(basismat)

}
