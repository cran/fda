\name{expon}
\alias{expon}
\title{Exponential Basis Function Values }
\description{
Evaluates a set of exponential basis functions, or a derivative of these functions, at a set of arguments. 
}
\usage{
expon(x, ratevec=1, nderiv=0)
}
\arguments{
  \item{x }{A vector of argument values at which the exponential basis functions are to evaluated. }
  \item{ratevec }{A vector of rate values defining the exponential basis functions. If b is such a rate value, the corresponding basis function is ebx. The number of basis functions is equal to the number of rate constants. The default is the single constant one. }
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
