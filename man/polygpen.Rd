\name{polygpen}
\alias{polygpen}
\title{Polygonal Penalty Matrix }
\description{
Computes the matrix defining the roughness penalty for functions expressed in terms of a polygonal basis. 
}
\usage{
polygpen(basisfd, Lfd=1)
}
\arguments{
  \item{basisfd }{A polygonal basis object. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. }
}
\details{
A roughness penalty for a function x(t) is defined by integrating the square of either the derivative of x(t) or, more generally, the result of applying a linear differential operator L to it. The only roughness penalty possible aside from penalizing the size of the function itself is the integral of the square of the first derivative D1x(t), and this is the default. To apply this roughness penalty, the matrix of inner products produced by this function is necessary. 
}
\value{
A symmetric matrix of order equal to the number of basis functions defined by the polygonal basis object. Each element is the inner product of two polygonal basis functions after applying the derivative or linear differential operator defined by Lfdobj. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
