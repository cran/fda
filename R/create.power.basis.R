create.power.basis <- function (rangeval=c(0,1), exponent=1)
{

  #  This function creates an power functional data basis
  #  Arguments
  #  RANGEVAL ... An array of length 2 containing the lower and upper
  #               boundaries for the rangeval of argument values
  #  NBASIS   ... The number of basis functions.  If this conflicts with
  #               the length of exponent, the latter is used.
  #  EXPONENT  ... The rate parameters defining x^exponent[i]
  #  Returns
  #  BASIS.FD  ... a functional data basis object of type "power"

  #  Last modified 25 March 2003

  type     <- "power"
  if (!rangechk(rangeval)) stop('Argument RANGEVAL is not correct.')
  if (!is.numeric(exponent)) stop('Argument exponent parameter is not correct')
  params   <- sort(as.vector(exponent))
  nbasis   <- length(params)
  if (nbasis<=0) stop('Argument exponent parameter is not correct')
  if (nbasis>1) {
      for(i in 1:(nbasis-1)){
         if (params[i]==params[i+1])
            stop('element in exponent should not equal to each other')
      }
   }

  #  set up basis object

  basisfd        <- list(type, rangeval, nbasis, params)
  names(basisfd) <- c("type", "rangeval", "nbasis", "params")
  class(basisfd) <- "basis.fd"

  return(basisfd)
}
