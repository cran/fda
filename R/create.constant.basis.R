create.constant.basis <- function(rangeval = c(0,1))
{
  #  This function creates a constant basis
  #  Argument:
  #  RANGEVAL ... an array of length 2 containing the lower and upper
  #  Return:
  #  BASISFD  ... a functional data basis object of type "constant"
  #

  #  Last modified 25 March 2003

  type <- "constant"

  if (!rangechk(rangeval)) stop('Argument RANGEVAL is not correct.')

  params  <- 0

  nbasis  <- 1

  #  set up basis object

  basisfd <- list(type, rangeval, nbasis, params)
  names(basisfd) <- c("type", "rangeval", "nbasis", "params")

  class(basisfd) <- "basis.fd"

  return(basisfd)
}
