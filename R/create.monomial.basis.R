create.monomial.basis <- function(rangeval=c(0,1), nbasis=2, exponents=NULL,
                                  dropind=NULL, quadvals=NULL, values=NULL)
{
#  CREATE_MONOMIAL_BASIS  Creates a monomial basis:, x^i_1, x^i_2, ...
#  The exponents in this version must be nonnegative integers
#  Argument:
#  RANGEVAL  ... an array of length 2 containing the lower and upper
#                boundaries for the rangeval of argument values.  If a
#                single value is input, it must be positive and the lower
#                limit of the range is set to 0.
#  NBASIS    ... number of basis functions
#  EXPONENTS ... an array of NBASIS nonnegative integer exponents
#                by default this is 0:(NBASIS-1)
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
#  Return:
#  BASISOBJ  ... a functional data basis object of type "monom"
#
#  last modified 20 November 2005

#  check RANGEVAL

if (length(rangeval) == 1) {
    if( rangeval <= 0) stop("RANGEVAL a single value that is not positive.")
    rangeval <- c(0,rangeval)
}
if (!rangechk(rangeval)) stop("RANGEVAL is not a legitimate range.")

if (missing(exponents) || is.null(exponents)) exponents = c(0:(nbasis-1))

#  check whether exponents are nonnegative integers

for(ibasis in 1:nbasis) {
	if (exponents[ibasis] - round(exponents[ibasis]) != 0)
		stop("An exponent is not an integer.")
	if (exponents[ibasis] < 0)  stop("An exponent is negative.")
}

# check if there are duplicate exponents

if (min(diff(sort(exponents))) <= 0) stop("There are duplicate exponents.")

type   = "monom"
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
