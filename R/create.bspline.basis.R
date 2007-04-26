create.bspline.basis <- function (rangeval=c(0,1), nbasis=NULL, norder=4,
                                  breaks=NULL, dropind=NULL, quadvals=NULL,
                                  values=NULL, names="bspl")
{

#  This function creates a bspline functional data basis.
#  Arguments
#  RANGEVAL ... an array of length 2 containing the lower and upper
#               boundaries for the rangeval of argument values
#  NBASIS   ... the number of basis functions
#  NORDER   ... order of b-splines (one higher than their degree).  The
#                 default of 4 gives cubic splines.
#  BREAKS   ... also called knots, these are a strictly increasing sequence
#               of junction points between piecewise polynomial segments.
#               They must satisfy BREAKS[1] = RANGEVAL[1] and
#               BREAKS[NBREAKS] = RANGEVAL[2], where NBREAKS is the total
#               number of BREAKS.  There must be at least 3 BREAKS.
#  There is a potential for inconsistency among arguments NBASIS, NORDER, and
#  BREAKS.  It is resolved as follows:
#     If BREAKS is supplied, NBREAKS = length(BREAKS), and
#     NBASIS = NBREAKS + NORDER - 2, no matter what value for NBASIS is
#     supplied.
#     If BREAKS is not supplied but NBASIS is, NBREAKS = NBASIS - NORDER + 2,
#        and if this turns out to be less than 3, an error message results.
#     If neither BREAKS nor NBASIS is supplied, NBREAKS is set to 21.
#  DROPIND  ... A vector of integers specificying the basis functions to
#               be dropped, if any.  For example, if it is required that
#               a function be zero at the left boundary, this is achieved
#               by dropping the first basis function, the only one that
#               is nonzero at that point.
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
#  BASISFD  ... a functional data basis object

# Last modified 3 March 2007 by Spencer Graves
#  Last modified 20 November 2005

  type <- "bspline"

#  check RANGE

  if (length(rangeval) == 1){
    if(rangeval <= 0) stop("RANGEVAL a single value that is not positive.")
    rangeval = c(0,rangeval)
  }
  if (!rangechk(rangeval)) stop("Argument RANGEVAL is not correct.")

#  check NORDER

  if (!is.numeric(norder) || norder<=0) stop("norder must be positive integer")
  norder <- ceiling(norder)

#  check for various argument configurations

# case of complete argument
  if (!is.null(nbasis) & !is.null(breaks)) {
    nbreaks <- length(breaks)
    if(!is.numeric(nbasis) || nbasis<=0)
      stop("Argument nbasis is not correct")
    else nbasis <- ceiling(nbasis)
  }


# case of NULL NBASIS

  if (is.null(nbasis) && !is.null(breaks)) {
    nbreaks <- length(breaks)
    nbasis  <- nbreaks + norder - 2
    if (nbasis<=0) stop("not enough break points")
  }

# case of NULL BREAKS

  if (!(is.null(nbasis)) && is.null(breaks)) {
    if(!is.numeric(nbasis) || nbasis<=0) stop("Argument nbasis is not correct")
    else nbasis <- ceiling(nbasis)
    nbreaks <- nbasis - norder + 2
    if (nbreaks < 2) stop("Number of knots is less than 2")
    breaks  <- seq(rangeval[1], rangeval[2], len=nbreaks)
  }

# case of NULL NBASIS and NULL BREAKS

  if (is.null(nbasis) && is.null(breaks)) {
    nbreaks <- 21
    nbasis  <- 19 + norder
    breaks  <- seq(rangeval[1], rangeval[2], len=nbreaks)
  }

#  check NBREAKS

  if (nbreaks < 2) stop ("Number of values in BREAKS less than 2.")

#  check BREAKS

  breaks <- sort(breaks)
  if (breaks[1] != rangeval[1]) stop(
              "Smallest value in BREAKS not equal to RANGEVAL[1].")
  if (breaks[nbreaks] != rangeval[2])   stop(
              "Largest  value in BREAKS not equal to RANGEVAL[2].")
  if (min(diff(breaks)) < 0) stop(
              "Values in BREAKS not increasing")

#  Set up the PARAMS vector.
#  PARAMS contains only the interior knots

  if (nbreaks > 2) {
    params <- breaks[2:(nbreaks-1)]
  } else {
    params <- NULL
  }

#  check DROPIND

  if (missing(dropind))  dropind = NULL
  if (length(dropind) > 0) {
    if(length(dropind) >= nbasis)  stop("Too many index values in DROPIND.")
    dropind <- sort(dropind)
    if(length(dropind) > 1) {
      if(min(diff(dropind)) == 0) stop("Multiple index values in DROPIND.")
    }
    for(i in 1:(length(dropind))) {
      if((dropind[i] < 1) || (dropind[i] > nbasis))
        stop("An index value of dropind is out of range.")
    }
  }

  basisobj <- basisfd(type=type, rangeval=rangeval, nbasis=nbasis, params=params,
                    dropind=dropind, quadvals=quadvals, values=values)
  basisobj$names <- {
    if(length(names) == nbasis) names
    else paste(names, norder, ".", 1:nbasis, sep="")
  }
  basisobj

}
