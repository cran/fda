create.polynomial.basis <- function (rangeval=c(0,1), nbasis=2, ctr=midrange)
{

  #  This function creates a polynomial functional data basis, for
  #    polynomials of the form  (x - c)^j
  #  Arguments
  #  RANGEVAL ... an array of length 2 containing the lower and upper
  #               boundaries for the rangeval of argument values
  #  NBASIS   ... the number of basis functions
  #  CTR      ... The centering constant C.  By default, this is the mid-range
  #  Returns
  #  BASISFD  ... a functional data basis object of type "polynomial"

  #  Last modified 25 March 2003

  type     <- "polynomial"

  if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")
  midrange <- mean(rangeval)
  params   <- as.vector(ctr)
  if ((length(params)<=0) || !is.numeric(params))
    stop("Argument ctr is not correct")
  nbasis <- ceiling(nbasis)
  if (nbasis<=0) stop ("nbasis must be positive integer")

  #  set up basis object

  basisfd        <- list(type, rangeval, nbasis, params)
  names(basisfd) <- c("type", "rangeval", "nbasis", "params")
  class(basisfd) <- "basis.fd"

  return(basisfd)
}
