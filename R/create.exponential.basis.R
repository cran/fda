create.exponential.basis <- function (rangeval=c(0,1), nbasis=1, ratevec=1,
                                      dropind=NULL, quadvals=NULL,
                                      values=NULL)
{

#  This function creates an exponential functional data basis
#  Arguments
#  RANGEVAL ... An array of length 2 containing the lower and upper
#               boundaries for the rangeval of argument values
#  NBASIS   ... The number of basis functions.  If this conflicts with
#               the length of RATEVEC, the latter is used.
#  RATEVEC  ... The rate parameters defining exp(ratevec[i]*x)
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
#  BASISOBJ  ... a functional data basis object of type "exponential"

#  Last modified 8 December 2005

type <- "expon"

if (length(rangeval) == 1) {
    if(rangeval <= 0)  stop("RANGEVAL a single value that is not positive.")
    rangeval <- c(0,rangeval)
}

if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")

if (length(ratevec) != nbasis) {
    warning("NBASIS is not equal to length of RATEVEC, and is revised.")
    nbasis <- length(ratevec)
}
params <- as.vector(ratevec)

if(missing(nbasis))  nbasis = 1
if(missing(ratevec)) ratevec = c(0:(nbasis-1))

if (length(ratevec)>1 && (min(diff(sort(ratevec)))<= 0))
	stop("There are duplicate ratevec.")

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
