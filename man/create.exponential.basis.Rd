\name{create.exponential.basis}
\alias{create.exponential.basis}
\title{
  Create an Exponential Basis
}
\description{
  Create an exponential basis object defining a set of exponential
  functions with rate constants in argument ratevec.
}
\usage{
create.exponential.basis(rangeval=c(0,1), nbasis=1,
                         ratevec=1, dropind=NULL,
                         quadvals=NULL, values=NULL)
}
\arguments{
  \item{rangeval}{
    a vector of length 2 containing the initial and final values of the
    interval over which the functional data object can be evaluated.
  }
  \item{nbasis}{
    the number of basis functions.
  }

  \item{ratevec}{
    a vector of length \code{nbasis} of rate constants
    defining basis functions of the form \code{exp(rate*x)}.
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
  a basis object with the type \code{expon}.
}
\details{
  Exponential functions are of the type $exp(bx)$ where $b$
  is the rate constant.  If $b = 0$, the constant function is
  defined.
}
\seealso{
  \code{\link{basisfd}}, 
  \code{\link{create.bspline.basis}}, 
  \code{\link{create.constant.basis}}, 
  \code{\link{create.fourier.basis}}, 
  \code{\link{create.monomial.basis}}, 
  \code{\link{create.polygonal.basis}}, 
  \code{\link{create.polynomial.basis}}, 
  \code{\link{create.power.basis}}
}
\examples{

#  Create an exponential basis over interval [0,5]
#  with basis functions 1, exp(-t) and exp(-5t)
basisobj <- create.exponential.basis(c(0,5),3,c(0,-1,-5))
#  plot the basis
plot(basisobj)

}
% docclass is function
\keyword{smooth}
