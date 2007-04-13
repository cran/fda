\name{create.polygonal.basis}
\alias{create.polygonal.basis}
\title{
  Create a Polygonal Basis
}
\description{
  A basis is set up for constructing polygonal lines, consisting of 
  straight line segments that join together.
}
\usage{
create.polygonal.basis(argvals=NULL, dropind=NULL,
                       quadvals=NULL, values=NULL)
}
\arguments{
  \item{argvals}{
    a strictly increasing vector of argument values at which line
    segments join to form a polygonal line.
  }
  \item{dropind}{
    a vector of integers specifiying the basis functions to
    be dropped, if any.  For example, if it is required that
    a function be zero at the left boundary, this is achieved
    by dropping the first basis function, the only one that
    is nonzero at that point. Default value NULL.
  }
  \item{quadvals}{
    a matrix with two columns and a number of rows equal to the number
    of argument values used to approximate an integral using Simpson's
    rule.  The first column contains these argument values.  
    A minimum of 5 values are required for
    each inter-knot interval, and that is often enough. These
    are equally spaced between two adjacent knots.
    The second column contains the weights used for Simpson's
    rule.  These are proportional to 1, 4, 2, 4, ..., 2, 4, 1.
  }
  \item{values}{
    a list containing the basis functions and their derivatives
    evaluated at the quadrature points contained in the first
    column of \code{ quadvals }.
  }
}
\value{
  a basis object with the type \code{polyg}.
}
\details{
  The actual basis functions consist of triangles, each with its apex
  over an argument value. Note that in effect the polygonal basis is
  identical to a B-spline basis of order 2 and a knot or break value at
  each argument value.  The range of the polygonal basis is set to the
  interval defined by the smallest and largest argument values.
}
\seealso{
\code{\link{basisfd}}, 
\code{\link{create.bspline.basis}}, 
\code{\link{create.constant.basis}}, 
\code{\link{create.exponential.basis}}, 
\code{\link{create.fourier.basis}}, 
\code{\link{create.monomial.basis}}, 
\code{\link{create.polynomial.basis}}, 
\code{\link{create.power.basis}}
}
\examples{
#  Create a polygonal basis over the interval [0,1]
#  with break points at 0, 0.1, ..., 0.95, 1
basisobj <- create.polygonal.basis(seq(0,1,0.1))
#  plot the basis
plot(basisobj)
}
% docclass is function
\keyword{smooth}
