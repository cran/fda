\name{bsplineS}
\alias{bsplineS}
\title{
  B-spline Basis Function Values
}
\description{
Evaluates a set of B-spline basis functions, or a derivative of these
functions, at a set of arguments.
}
\usage{
bsplineS(x, breaks, norder=4, nderiv=0, returnMatrix=FALSE)
}
\arguments{
  \item{x}{
    A vector of argument values at which the B-spline basis functions
    are to be evaluated.
  }
  \item{breaks}{
    A strictly increasing set of break values defining the B-spline
    basis.  The argument values \code{x} should be within the interval
    spanned by the break values.
  }
  \item{norder}{
    The order of the B-spline basis functions.  The order less one is
    the degree of the piece-wise polynomials that make up any B-spline
    function. The default is order 4, meaning piece-wise cubic.
  }
  \item{nderiv}{
    A nonnegative integer specifying the order of derivative to be
    evaluated.  The derivative must not exceed the order.  The default
    derivative is 0, meaning that the basis functions themselves are
    evaluated.
  }
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
}
\value{
  a matrix of function values.  The number of rows equals the number of
  arguments, and the number of columns equals the number of basis
  functions.
}
\examples{
# Minimal example:  A B-spline of order 1 (i.e., a step function)
# with 0 interior knots:
bS <- bsplineS(seq(0, 1, .2), 0:1, 1, 0)

# check
\dontshow{stopifnot(}
all.equal(bS, matrix(1, 6))
\dontshow{)}

#  set up break values at 0.0, 0.2,..., 0.8, 1.0.
breaks <- seq(0,1,0.2)
#  set up a set of 11 argument values
x <- seq(0,1,0.1)
#  the order willl be 4, and the number of basis functions
#  is equal to the number of interior break values (4 here)
#  plus the order, for a total here of 8.
norder <- 4
#  compute the 11 by 8 matrix of basis function values
basismat <- bsplineS(x, breaks, norder)

# use sparse Matrix representation to save memory
Basismat <- bsplineS(x, breaks, norder, returnMatrix=TRUE)

# check
\dontshow{stopifnot(}
class(Basismat)=='dgCMatrix'
\dontshow{)}

\dontshow{stopifnot(}
all.equal(basismat, as.matrix(Basismat))
\dontshow{)}

}
\keyword{smooth}
