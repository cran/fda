polynom <- function (x, norder=1, nderiv=0, ctr=midrange)
{
#  This computes values of the polynomials,
#        P_l(x) = (x-ctr)^l, l=0,...,NORDER-1
#  or their derivatives.
#  The degree of the highest order polynomial is one less than NORDER.
#  The default is the constant function.

#  Arguments are as follows:
#  X      ... array of values at which the polynomials are to
#             evaluated
#  NORDER ... highest degree plus one
#  NDERIV ... highest order derivative.  0 means only function values
#             are returned.
#  CTR    ... a constant shift that helps to keep the polynomials from
#             getting too ill-conditioned.  A good choice is the mid-range.
#  Return is a matrix with length(X) rows and NORDER columns containing
#  the values of the polynomials

#  last modified 8 June 1999

  x        <- as.vector(x)
  n        <- length(x)
  ndegree  <- norder - 1
  if (nderiv > ndegree) stop('NDERIV exceeds highest degree of polynomial.')
  rangex   <- range(x)
  midrange <- mean(rangex)
  lfac <- 1
  if (nderiv > 1) for (l in 2:nderiv) lfac <- lfac*l
  polyval <- matrix(0,n,norder)
  polyval[,nderiv+1] <- lfac
  if (norder > nderiv+1)
    for (l in (nderiv+2):norder)
       polyval[,l] <- polyval[,l-1]*(x-ctr)*(l-1)/(l-nderiv-1)
  return (polyval)
}
