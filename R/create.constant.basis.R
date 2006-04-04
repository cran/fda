create.constant.basis <- function(rangeval = c(0,1))
{
#  This function creates a constant basis
#  Argument:
#  RANGEVAL ... an array of length 2 containing the lower and upper
#  Return:
#  BASISOBJ  ... a functional data basis object of type "constant"
#

#  Last modified 8 December 2005

type     <- "const"
nbasis   <- 1
params   <- NULL
dropind  <- NULL
quadvals <- NULL
values   <- NULL

if (length(rangeval) == 1){
    if (rangeval <= 0) stop("RANGEVAL a single value that is not postive.")
    rangeval = c(0,rangeval)
}

if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")

basisobj <- basisfd(type=type, rangeval=rangeval, nbasis=nbasis, params=params,
                    dropind=dropind, quadvals=quadvals, values=values)

basisobj

}
