\name{var.fd}
\alias{var.fd}
\title{Variance and Covariance Surfaces for a Functional Data Object }
\description{
Compute the variance and covariance functions for functional data. These are two-argument functions, and therefore define a surface. If there is one functional data object supplied, its variance is computed. If there are two, the covariance between the two is computed. Both functional data objects must be univariate. 
}
\usage{
var.fd(fdx, fdy=fdx)
}
\arguments{
  \item{fdx }{A univariate functional data object. }
  \item{fdy }{A univariate functional data object. }
}
\details{
}
\value{
A two-argument functional data object representing the variance surface for a single functional data object or the covariance between two functional data objects. 
}
\references{}
\note{}
\author{}
\seealso{
meanFd, std.fd 
}
\examples{
}
\keyword{}
