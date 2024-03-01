bsplineS <- function (x, breaks, norder=4, nderiv=0, returnMatrix=FALSE)
{
#  This is a wrapper function for the S-PLUS spline.des function.
#  The number of spline functions is equal to the number of
#     discrete break points, length(BREAKVALUES), plus the order, NORDER,
#           minus 2.
#  Arguments are as follows:
#  X      ... array of values at which the spline functions are to evaluated
#  BREAKS ... a STRICTLY INCREASING sequence of break points or knot
#             values.  It must contain all the values of X within its
#             range.
#  NORDER ... order of spline (1 more than degree), so that 1 gives a
#             step function, 2 gives polygonal functions, 3 gives
#             quadratic functions between knots and default 4 gives 
#             cubic splines
#  NDERIV ... highest order derivative.  default 0 means function are returned
#  Return is a matrix with length(X) rows and number of columns equal to
#                   number of b-splines

#  previously modified 31 May 2023 by Jim Ramsay

  x <- as.vector(x)
  n <- length(x)  # number of evaluation values
  tol <- 1e-14
  nbreaks <- length(breaks)  # number of knots, including end points
  if (nbreaks < 2) stop('Number of knots less than 2.')
  #  check that breaks are strictly increasing
  if (min(diff(breaks)) < 0 ) stop('Knots are not increasing')
  #  check that all evaluation values are within break sequence
  if ( max(x) > max(breaks) + tol ||min(x) < min(breaks) - tol ) 
    stop('Knots do not span the values of X')
  #  modify break limits if small deviations from interval exist
  if ( x[n] > breaks[nbreaks]) breaks[nbreaks] <- x[n]
  if ( x[1] < breaks[1]      ) breaks[1]       <- x[1]
  #  check that norder is within [1,20]
  if (norder > 20) stop('NORDER exceeds 20.')
  if (norder <  1) stop('NORDER less than 1.')
  #  check that order of derivative with [1,19]
  if (nderiv > 19) stop('NDERIV exceeds 19.')
  if (nderiv <  0) stop('NDERIV is negative.')
  # check that order of drivative is less than order of spline
  # if (nderiv >= norder) 
  #     stop ('NDERIV cannot be as large as order of B-spline.')
  #  compute the number of basis functions
  nbasis <- nbreaks + norder - 2
  #  if nderiv not less than norder, return a matrix of zeros
  #  note: this is not strictly legitimate since if an 
  #  evaluation value falls on a break value, the output value 
  #  should be NA instead of 0 since it does not exist
  if (nderiv >= norder) {
    return(matrix(0,n,nbasis))
  }
  #  now set up the whole knot sequence including norder
  #  equal initial and final and final values
  knots  <- c(rep(breaks[1      ],norder-1), breaks,
              rep(breaks[nbreaks],norder-1)  )
  #  numeric object to contain the derivative values
  derivs <- rep(nderiv,n)
  #  check that the number of basis values is not less than
  #  the order of the spline
  if (nbasis >= norder) {
    if (nbasis > 1) {
      basismat <-           spline.des(knots, x, norder, derivs)$design
    } else {
      basismat <- as.matrix(spline.des(knots, x, norder, derivs)$design)
    }
    if((!returnMatrix) && (length(dim(basismat)) == 2)) {
      return(as.matrix(basismat))
    } else {
      return(basismat)
    }
  } else {
    stop("NBASIS is less than NORDER.")
  }
}
