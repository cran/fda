\name{fourierpen}
\alias{fourierpen}
\title{Fourier Penalty Matrix }
\description{
Computes the matrix defining the roughness penalty for functions expressed in terms of a Fourier basis. 
}
\usage{
fourierpen(basisobj, Lfdobj=2)
}
\arguments{
  \item{basisobj }{A Fourier basis object. }
  \item{Lfdobj }{Either a nonnegative integer or a linear differential operator object. }
}
\details{
A roughness penalty for a function x(t) is defined by integrating the square of either the derivative of x(t) or, more generally, the result of applying a linear differential operator L to it. The most common roughness penalty is the integral of the square of the second derivative D2x(t), and this is the default. To apply this roughness penalty, the matrix of inner products of the basis functions (possibly after applying the linear differential operator to them) defining this function is necessary. This function just calls the roughness penalty evaluation function specific to the basis involved. 
}
\value{
A symmetric matrix of order equal to the number of basis functions defined by the Fourier basis object. Each element is the inner product of two Fourier basis functions after applying the derivative or linear differential operator defined by Lfdobj. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
