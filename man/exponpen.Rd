\name{exponpen}
\alias{exponpen}
\title{Exponential Penalty Matrix }
\description{
Computes the matrix defining the roughness penalty for functions expressed in terms of a exponential basis. 
}
\usage{
exponpen(basisfd, Lfd=2)
}
\arguments{
  \item{basisfd}{An exponential basis object. }
  \item{Lfd}{Either a nonnegative integer or a linear differential operator object. }
}
\details{
A roughness penalty for a function x(t) is defined by integrating the square of either the derivative of x(t) or, more generally, the result of applying a linear differential operator L to it. The most common roughness penalty is the integral of the square of the second derivative D2x(t), and this is the default. To apply this roughness penalty, the matrix of inner products produced by this function is necessary. 
}
\value{
A symmetric matrix of order equal to the number of basis functions defined by the exponential basis object. Each element is the inner product of two exponential basis functions after applying the derivative or linear differential operator defined by Lfdobj. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
