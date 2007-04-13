\name{eval.basis}
\alias{eval.basis}
\title{
  Values of Basis Functions or their Derivatives
}
\description{
  A set of basis functions are evaluated at a vector of
  argument values.  If a linear differential object is provided, the 
  values are the result of applying the the operator to each basis
  function. 
}
\usage{
eval.basis(evalarg, basisobj, Lfdobj=0)
}
\arguments{
  \item{evalarg}{
    a vector of argument values.
  }
  \item{basisobj}{
    a basis object defining basis functions whose values
    are to be computed.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential.
    operator object.
  }
}
\value{
  a matrix of basis function values with rows corresponding
  to argument values and columns to basis functions.
}
\details{
  If a linear differential operator object is supplied, the basis must
  be such that the highest order derivative can be computed. If a
  B-spline basis is used, for example, its order must be one larger than
  the highest order of derivative required.
}
\seealso{
  \code{\link{getbasismatrix}}, 
  \code{\link{eval.fd}}, 
  \code{\link{plot.basisfd}}
}
% docclass is function
\keyword{smooth}
