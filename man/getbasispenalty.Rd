\name{getbasispenalty}
\alias{getbasispenalty}
\title{Evaluate a Roughness Penalty Matrix }
\description{
Computes the matrix defining the roughness penalty for functions expressed in terms of a basis. 
}
\usage{
getbasispenalty(basisfd, Lfd=NULL)
}
\arguments{
  \item{basisfd }{A basis object. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. }
}
\details{
A roughness penalty for a function x(t) is defined by integrating the square of either the derivative of x(t) or, more generally, the result of applying a linear differential operator L to it. The most common roughness penalty is the integral of the square of the second derivative D2x(t), and this is the default. To apply this roughness penalty, the matrix of inner products of the basis functions defining this function is necessary. This function just calls the roughness penalty evaluation function specific to the basis involved. 
}
\value{
A symmetric matrix of order equal to the number of basis functions defined by the B-spline basis object. Each element is the inner product of two B-spline basis functions after applying the derivative or linear differential operator defined by Lfd. 
}
\references{}
\note{}
\author{}
\seealso{
bsplinepen, exponpen, fourierpen, polygpen, powerpen 
}
\examples{}
\keyword{}
