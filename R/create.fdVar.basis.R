
#  ------------------------------------------------------------------------

create.fdVar.basis <- function(rangeval, I, J=I) {
#  CREATE_FDVARIANCE_BASIS Creates a functional basis object with type 
#  'fdVariance' thatcan be used to define a covariance surface.  This 
#  object can consist of one or more trapezoidal domains defined over a 
#  common range.  The range is defined in argument RANGEVAL, which is also
#  the RANGEVAL field of the basis object.
#  Each trapezoidal region is defined by the number of vertical
#  intervals in the corresponding element in I, and the number of
#  horizontal elements in the corresponding element in J.  I and J
#  must have the same lengths, and their common length is the number
#  of trapezoidal domains defined.
#  The PARAMS field for the basis object is a struct object with fields
#  I and J.  The TYPE field is 'fdVariance'.  The NBASIS field is the
#  sum over i of (I(i)+1)*(J(i)+1)
#
#  Arguments:

#  Last modified 102 June 2015 by Jim Ramsay

#  check I and J

nI = length(I)
nJ = length(J)

# check I and J

if (nI != nJ)           stop('I and J do not have same lengths.')
if (any(I <= 0))        stop('I has zero or negative entries.')
if (any(J <  0))        stop('I has negative entries.')
if (any(floor(I) != I)) stop('I has non-integer values.')
if (any(floor(J) != J)) stop('J has non-integer values.')

#  check RANGEVAL

if (length(rangeval) == 1) {
    if (rangeval <= 0) stop('RANGEVAL is a single value that is not positive.')
    rangeval = c(0,rangeval)
}
if (rangechk(rangeval) != 1) stop('RANGEVAL is not a legitimate range.')

#  get sum of number of basis functions

nbasis = 0
for (i in 1:nI) nbasis = nbasis + (I[i]+1)*(J[i]+1)

#  construct basis object

type = 'fdVariance'

params = list(I=I, J=J)

basisobj = basisfd(type, rangeval, nbasis, params)

return(basisobj)

}
