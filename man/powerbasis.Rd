\name{powerbasis}
\alias{powerbasis}
\title{Power Basis Function Values }
\description{
Evaluates a set of power basis functions, or a derivative of these functions, at a set of arguments. 
}
\usage{
powerbasis(x, exponents, nderiv=0)
}
\arguments{
  \item{x }{A vector of argument values at which the power basis functions are to evaluated. Since exponents may be negative, for example after differentiation, it is required that all argument values be positive. }
  \item{exponents }{A vector of exponents defining the power basis functions. If y is such a rate value, the corresponding basis function is xy. The number of basis functions is equal to the number of exponents. }
  \item{nderiv }{The derivative to be evaluated. The derivative must not exceed the order. The default derivative is 0, meaning that the basis functions themselves are evaluated. }
}
\details{
}
\value{
A matrix of function values. The number of rows equals the number of arguments, and the number of columns equals the number of basis 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
