\name{lambda2df}
\alias{lambda2df}
\title{Convert Smoothing Parameter to Degrees of Freedom }
\description{
The degree of roughness of an estimated function is controlled by a smoothing parameter lambda that directly multiplies the penalty. However, it can be difficult to interpret or choose this value, and it is often easier to determine the roughness by choosing a value that is equivalent of the degrees of freedom used by the smoothing procedure. This function converts a multipler lambda into a degrees of freedom value. 
}
\usage{
lambda2df(argvals, basisfd, wtvec=rep(1, n), Lfd=NULL, lambda=0)
}
\arguments{
  \item{argvals }{Argument values associated with the values to be smoothed. }
  \item{basisfd }{A basis object. }
  \item{wtvec }{A vector of weights for the data to be smoothed. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. }
  \item{lambda }{The smoothing parameter to be converted. }
}
\details{
The conversion requires a one-dimensional optimization and may be therefore computationally intensive. 
}
\value{
The equivalent degrees of freedom value. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{
}
\keyword{}
