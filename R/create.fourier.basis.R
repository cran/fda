create.fourier.basis <- function (rangeval=c(0,1), nbasis=3,
          period=width, dropind=NULL, quadvals=NULL, values=NULL,
          longNames=TRUE)
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
#  DROPIND  ... A vector of integers specificying the basis functions to
#               be dropped, if any.  
#  QUADVALS ... A matrix with two columns and a number of rows equal to
#               the number of argument values used to approximate an 
#               integral using Simpson's rule.  
#               The first column contains these argument values.  
#               A minimum of 5 values are required for
#               each inter-knot interval, and that is often enough. These
#               are equally spaced between two adjacent knots.
#               The second column contains the weights used for Simpson's
#               rule.  These are proportional to 1, 4, 2, 4, ..., 2, 4, 1
#   VALUES  ... A list containing the basis functions and their derivatives
#               evaluated at the quadrature points contained in the first 
#               column of QUADVALS.  
#  Returns
#  BASISOBj  ... a functional data basis object of type "fourier"

#  Last modified 21 February 2007 by Spencer Graves
#  Previously modified 20 November 2005

  type <- "fourier"

  if (length(rangeval)==1) {
    if (rangeval<=0) stop("RANGEVAL a single value that is not positive.")
    rangeval <- c(0,rangeval)
  }

  if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")
  
  width <- rangeval[2] - rangeval[1]
  if ((period <= 0) || !is.numeric(period))
    stop ("Period must be positive number for a Fourier basis")
  params <- period

#  increase the number of basis functions by one if even

  if ((nbasis <= 0) || !is.numeric(nbasis))
    stop ("nbasis must be odd positive number for a Fourior basis")
  nbasis <- ceiling(nbasis)
  if (2*floor(nbasis/2) == nbasis) nbasis <- nbasis + 1

#  check DROPIND

  if (missing(dropind)) dropind <- NULL

  if (length(dropind) > 0){
    if(length(dropind) >= nbasis)  stop("Too many index values in 'dropind'.")
    dropind = sort(dropind)
    if(length(dropind) > 1) {
      if(min(diff(dropind)) == 0) stop("Duplicate index values in 'dropind'.")
    }
    for(i in 1:length(dropind)) {
      if(dropind[i] < 1 || dropind[i] > nbasis)
        stop("An index value is out of range.")
    }
  }
  
#  set up the basis object

  basisobj <- basisfd(type=type, rangeval=rangeval, nbasis=nbasis,
               params=params, dropind=dropind, quadvals=quadvals,
                      values=values)
# Names?
  if(!is.na(longNames)){
    nb2 <- floor(nbasis/2)
    Nms <- "const"
    if(nb2>0){
      sinCos <- as.vector(outer(c("sin", "cos"), 1:nb2,
                                paste, sep="") )
      if(longNames)
        sinCos <- paste(sinCos, signif(period, 4), sep=".")
      Nms <- c(Nms, sinCos)
    }
#    
    basisobj$names <- Nms
  }
#  
  basisobj
}
