\name{fourier}
\alias{fourier}
\title{Fourier Basis Function Values }
\description{
Evaluates a set of Fourier basis functions, or a derivative of these functions, at a set of arguments. 
}
\usage{
fourier(x, nbasis=n, period=span, nderiv=0)
}
\arguments{
  \item{x }{A vector of argument values at which the Fourier basis functions are to evaluated. }
  \item{nbasis }{The number of basis functions in the Fourier basis. The first basis function is the constant function, followed by sets of sine/cosine pairs. Normally the number of basis functions will be an odd. The default number is the number of argument values. }
  \item{period }{The width of an interval over which all sine/cosine basis functions repeat themselves. The default is the difference between the largest and smallest argument values. }
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
\examples{
}
\keyword{}
