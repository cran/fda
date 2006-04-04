create.power.basis <- function (rangeval=c(0,1), nbasis=length(exponents), exponents=1,
                                dropind=NULL, quadvals=NULL, values=NULL)
{
#  This function creates an power functional data basis
#  Arguments
#  RANGEVAL ... An array of length 2 containing the lower and upper
#               boundaries for the rangeval of argument values
#  NBASIS   ... The number of basis functions.  If this conflicts with
#               the length of exponent, the latter is used.
#  EXPONENT  ... The rate parameters defining x^exponent[i]
#  DROPIND  ... A vector of integers specificying the basis functions to
#               be dropped, if any.  
#  QUADVALS ... A matrix with two columns and a number of rows equal to
#               the number of argument values used to approximate an 
#               integral using Simpson's rule.  
#               The first column contains these argument values.  
#               A minimum of 5 values are required for
#               each inter-knot interval, and that is often enough. These
#               are equally spaced between two adjacent knots.
#               The second column contains the weights used for Simpson's
#               rule.  These are proportional to 1, 4, 2, 4, ..., 2, 4, 1
#   VALUES  ... A list containing the basis functions and their derivatives
#               evaluated at the quadrature points contained in the first 
#               column of QUADVALS.  
#  Returns
#  BASISOBJ  ... a functional data basis object of type "power"

#  Last modified 20 November 2005

type  <- "power"

if (missing(nbasis))  nbasis = length(exponents)
if (missing(exponents)) exponents = c(0:(nbasis-1))
if (!is.numeric(exponents)) stop("Argument exponent parameter is not correct")
params <- sort(as.vector(exponents))
nbasis <- length(params)
if (nbasis<=0) stop("Argument exponent parameter is not correct")
if (nbasis>1) {
    for(i in 1:(nbasis-1)) {
       if (params[i]==params[i+1])
          stop("element in exponent should not equal to each other")
    }
}

if (length(rangeval) == 1){
    if (rangeval <= 0)   stop("RANGEVAL a single value that is not positive.")
    rangeval = c(0,rangeval)
}

if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")

if (rangeval[1] < 0)  stop("Lower limit in RANGEVAL is negative.")

params = exponents

#  check DROPIND

if (missing(dropind)) dropind <- NULL

if (length(dropind) > 0){
    if(length(dropind) >= nbasis)  stop("Too many index values in DROPIND.")
    dropind = sort(dropind)
    if(length(dropind) > 1) {
        if(min(diff(dropind)) == 0) stop("Multiple index values in DROPIND.")
    }
    for(i in 1:length(dropind)) {
   	    if(dropind[i] < 1 || dropind[i] > nbasis)
            stop("An index value is out of range.")
    }
}

#  set up the basis object

basisobj <- basisfd(type=type, rangeval=rangeval, nbasis=nbasis, params=params,
                    dropind=dropind, quadvals=quadvals, values=values)

basisobj

}
