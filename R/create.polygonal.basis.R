create.polygonal.basis <- function (argvals=NULL)
{

  #  This function creates a polygonal functional data basis.
  #  Arguments
  #  Returns
  #  BASISFD  ... a functional data basis object

  #  Last modified 25 March 2003

  type <- "polygonal"

  if(!is.vector(argvals) || !is.numeric(argvals))
    stop("Argument atgvals is not correct")

  nbasis <- length(argvals)
  if (nbasis < 2) stop("Number of ARGVALS less than two.")
  argvals <- sort(argvals)
  if (diff(argvals)==0)
    stop("element in argvals should not equal to each other")
  rangeval <- range(argvals)
  params   <- argvals

  #  set up basis object

  basisfd        <- list(type, rangeval, nbasis, params)
  names(basisfd) <- c("type", "rangeval", "nbasis", "params")
  class(basisfd) <- "basis.fd"

  return(basisfd)
}
