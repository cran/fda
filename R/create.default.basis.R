create.default.basis <- function(argvals, nresol, nderiv = 2, periodic = FALSE)
{
#  This function takes a vector or matrix argvals and creates a default
#   basis to be used for data observed on these arguments.
#
#  ARGVALS  ... A vector or matrix of argument values. Missing values allowed.
#
#  NRESOL   ... A number that specifies the number of the finest features or
#               events that are of interest that can occur within the
#               range of the argument values. By feature or event is
#               meant things like peaks, valleys, zero crossings,
#               plateaus, linear slopes, and so on.  NRESOL specifies
#               the amount of resolution required of the functional
#               data object.
#  NDERIV   ... A natural number, 0, 1, 2, ..., specifying the number
#               of derivatives that the functional data object must
#               have.  The default is 2.
#  PERIODIC ... If T, functions are treated as periodic, and in the
#               case of vector ARGVALS the
#               argument domain is extended below by one value to become
#                [ARGVALS[1] - (ARGVALS[N]-ARGVALS[1])/(N-1), ARGVALS[N].
#               The default is F.
#
#  Returns an object of class BASISFD, a functional data basis object

  #  Last modified 25 March 2003

  #  Check values used to set up basis

  if (is.matrix(argvals)) n <- dim(argvals)[1] else n <- length(argvals)
  rangeval <- range(argvals, na.rm = TRUE)
  if (is.integer(nresol) == FALSE) nresol <- as.integer(nresol)
  if(nresol < 1 || nresol > n) stop("NRESOL is not between 1 and N.")
  if (is.integer(nderiv) == FALSE) nderiv <- as.integer(nderiv)
  if (nderiv < 0 || nderiv > n - 1) stop("NDERIV is not between 0 and N-1.")
  if (is.logical(periodic) == FALSE) stop("PERIODIC is not a logical variable.")
  #  Set up basis object.
  if(periodic) {
    rangeval[1] <- rangeval[1] - diff(rangeval)/(n - 1)
    basisfd     <- create.fourier.basis(rangeval, nresol)
  } else {
    if(nresol == 1) {
      basisfd <- create.constant.basis(rangeval)
    } else {
      if(nderiv == 0 & nresol == n & !is.matrix(argvals)) {
        basisfd <- create.polygonal.basis(argvals)
      } else {
        basisfd <- create.bspline.basis(rangeval, nresol, nderiv + 2)
      }
    }
  }
  return(basisfd)
}
