create.exponential.basis <- function (rangeval=c(0,1), ratevec=1)
{

  #  This function creates an exponential functional data basis
  #  Arguments
  #  RANGEVAL ... An array of length 2 containing the lower and upper
  #               boundaries for the rangeval of argument values
  #  NBASIS   ... The number of basis functions.  If this conflicts with
  #               the length of RATEVEC, the latter is used.
  #  RATEVEC  ... The rate parameters defining exp(ratevec[i]*x)
  #  Returns
  #  BASIS.FD  ... a functional data basis object of type "exponential"

  #  Last modified 25 March 2003

  type     <- "exponential"

  if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")

  if(!is.numeric(ratevec)) stop("Rate vector should be numerical vector")

  params <- sort(as.vector(ratevec))
  nbasis <- length(params)
  if(nbasis <=0) stop("RATEVEC is empty.")
  if(nbasis>1) {
        for(i in 1:(nbasis-1)) {
             if (params[i]==params[i+1])
               stop("rate value should not equal to each other")
        }
   }

  #  set up basis object

  basisfd <- list(type, rangeval, nbasis, params)
  names(basisfd) <- c("type", "rangeval", "nbasis", "params")
  class(basisfd) <- "basis.fd"

  return(basisfd)
}
