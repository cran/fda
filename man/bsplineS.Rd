\name{bsplineS}
\alias{bsplineS}
\title{B-spline Basis Function Values}
\description{
Evaluates a set of B-spline basis functions, or a derivative of these functions, at a set of arguments. 
}
\usage{
bsplineS(x, breaks, norder=4, nderiv=0)
}
\arguments{
  \item{x}{A vector of argument values at which the B-spline basis functions are to evaluated.}
  \item{breaks}{A strictly increasing set of break values defining the B-spline basis. The argument values x should be within the interval spanned by the break values.}
  \item{norder}{The order of the B-spline basis functions. The order less one is the degree of the piece-wise polynomials that make up any B-spline function. The default is order 4, meaning piece-wise cubic. }
  \item{nderiv}{The derivative to be evaluated. The derivative must not exceed the order. The default derivative is 0, meaning that the basis functions themselves are evaluated.}
}
\details{
}
\value{
A matrix of function values. The number of rows equals the number of arguments, and the number of columns equals the number of basis functions. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
