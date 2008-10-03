\name{create.monomial.basis}
\alias{create.monomial.basis}
\title{
  Create a Monomial Basis
}
\description{
  Creates a set of basis functions consisting of powers
  of the argument.
}
\usage{
create.monomial.basis(rangeval=c(0, 1), nbasis=2,
                      exponents=NULL, dropind=NULL,
                      quadvals=NULL, values=NULL)
}
\arguments{
  \item{rangeval}{
    a vector of length 2 containing the initial and final
    values of the interval over which the functional
    data object can be evaluated.
  }
  \item{nbasis}{
    the number of basis functions. The default is 2,
    which defines a basis for straight lines.
  }
  \item{exponents}{
    the nonnegative integer powers to be used.  By default,
    these are 0, 1, 2, ..., (nbasis-1).
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
  a basis object with the type \code{monom}.
}
\seealso{
  \code{\link{basisfd}},
  \code{\link{create.bspline.basis}},
  \code{\link{create.constant.basis}},
  \code{\link{create.fourier.basis}},
  \code{\link{create.exponential.basis}},
  \code{\link{create.polygonal.basis}},
  \code{\link{create.polynomial.basis}},
  \code{\link{create.power.basis}}
}
\examples{
##
## simplest example: one constant 'basis function'
##
m0 <- create.monomial.basis(nbasis=1, exponents=0)
plot(m0)

##
## Create a monomial basis over the interval [-1,1]
##  consisting of the first three powers of t
##
basisobj <- create.monomial.basis(c(-1,1), 3)
#  plot the basis
plot(basisobj)
}
% docclass is function
\keyword{smooth}
