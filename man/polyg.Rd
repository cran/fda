\name{polyg}
\alias{polyg}
\title{Polygonal Basis Function Values }
\description{
Evaluates a set of polygonal basis functions, or a derivative of these functions, at a set of arguments. 
}
\usage{
polyg(x, argvals, nderiv=0)
}
\arguments{
  \item{x }{A vector of argument values at which the polygonal basis functions are to evaluated. }
  \item{argvals }{A strictly increasing set of argument values containing the range of x within it that defines the polygonal basis. The default is x itself. }
  \item{nderiv }{The derivative to be evaluated. The derivative must not exceed one. The default derivative is 0, meaning that the basis functions themselves are evaluated. }
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
