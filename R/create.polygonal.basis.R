create.polygonal.basis <- function (argvals=NULL, dropind=NULL,
                                    quadvals=NULL, values=NULL)
{

#  This function creates a polygonal functional data basis.
#  Arguments
#  ARGVALS  ... A strictly increasing vector of argument values at which 
#               line segments join to form a polygonal line.
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
#  BASISOBJ ... a functional data basis object

#  Last modified 20 November 2005

type <- "polyg"

if(!is.vector(argvals) || !is.numeric(argvals))
    stop("Argument atgvals is not correct")

nbasis <- length(argvals)
if (nbasis < 2) stop("Number of ARGVALS less than two.")
argvals <- sort(argvals)
if (min(diff(argvals))==0)
    stop("element in argvals should not equal to each other")
rangeval <- range(argvals)
params   <- argvals

#  check DROPIND

if (missing(dropind)) dropind <- NULL

if (length(dropind) > 0) {
    if(length(dropind) >= nbasis)  stop("Too many index values in DROPIND.")
    dropind = sort(dropind)
    if(length(dropind) > 1)
        if(min(diff(dropind)) == 0) stop("Multiple index values in DROPIND.")
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
