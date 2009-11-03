\name{getbasismatrix}
\alias{getbasismatrix}
\title{
  Values of Basis Functions or their Derivatives
}
\description{
Evaluate a set of basis functions or their derivatives at
a set of argument values.
}
\usage{
getbasismatrix(evalarg, basisobj, nderiv=0)
}
\arguments{
  \item{evalarg}{
    a vector of arguments values.
  }
  \item{basisobj}{
    a basis object.
  }
  \item{nderiv}{
    a nonnegative integer specifying the derivative to be evaluated.
  }
}
\value{
  a matrix of basis function or derivative values.  Rows correspond
  to argument values and columns to basis functions.
}
\seealso{
  \code{\link{eval.fd}}
}
\examples{
# Minimal example:  a B-spline of order 1, i.e., a step function
# with 0 interior knots:
bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)
getbasismatrix(seq(0, 1, .2), bspl1.1)

}
% docclass is function
\keyword{smooth}
