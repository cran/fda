create.fourier.basis <- function (rangeval=c(0,1), nbasis=3, period=width)
{

  #  This function creates a fourier functional data basis.
  #  Arguments
  #  RANGEVAL ... an array of length 2 containing the lower and upper
  #               boundaries for the rangeval of argument values
  #  NBASIS   ... the number of basis functions.  If the argument value is
  #               even, it is increased by one so both sines and cosines are
  #               present for each period.  A possible underdetermination of
  #               the basis is taken care of in function PROJECT.BASIS.
  #  PERIOD   ... The period.  That is, the basis functions are periodic on
  #                 the interval [0,PARAMS] or any translation of it.
  #  Returns
  #  BASISFD  ... a functional data basis object of type "fourier"

  #  Last modified 25 March 2003

  type <- "fourier"

  if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")

  width <- rangeval[2] - rangeval[1]
  if ((period <= 0) || !is.numeric(period))
    stop ("Period must be positive number for a Fourier basis")
  params <- period

  #  increase the number of basis functions by one if even

  if ((nbasis <= 0) || !is.numeric(nbasis)) stop ("nBasis must be odd positive number for a Fourior basis")
  nbasis <- ceiling(nbasis)
  if (2*floor(nbasis/2) == nbasis) nbasis <- nbasis + 1

  #  set up the basis object

  basisfd         <- list(type, rangeval, nbasis, params)
  names(basisfd)  <- c("type", "rangeval", "nbasis", "params")
   class(basisfd) <- "basis.fd"

  return(basisfd)
}
