\name{create.polynomial.basis}
\alias{create.polynomial.basis}
\title{
  Create a Polynomial Basis
}
\description{
  Creates a set of basis functions consisting of powers
  of the argument shifted by a constant.
}
\usage{
create.polynomial.basis(rangeval=c(0, 1), nbasis=2,
                        ctr=midrange, dropind=NULL,
                        quadvals=NULL, values=NULL)
}
\arguments{
  \item{rangeval}{
    a vector of length 2 defining the range.
  }
  \item{nbasis}{
    the number of basis functions. The default is 2,
    which defines a basis for straight lines.
  }
  \item{ctr}{
    this value is used to shift the argument prior to taking
    its power.
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
  a basis object with the type \code{polynom}.
}
\details{
  The only difference between a monomial and a polynomial basis
  is the use of a shift value.  This helps to avoid rounding error
  when the argument values are a long way from zero.
}
\seealso{
\code{\link{basisfd}}, 
\code{\link{create.bspline.basis}}, 
\code{\link{create.constant.basis}}, 
\code{\link{create.fourier.basis}}, 
\code{\link{create.exponential.basis}}, 
\code{\link{create.monomial.basis}}, 
\code{\link{create.polygonal.basis}}, 
\code{\link{create.power.basis}}
}
\examples{
#  Create a polynomial basis over the years in the 20th century
#  and center the basis functions on 1950.
basisobj <- create.polynomial.basis(c(1900, 2000), nbasis=3, ctr=1950)
#  plot the basis
plot(basisobj)
}
% docclass is function
\keyword{smooth}
