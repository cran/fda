\name{df2lambda}
\alias{df2lambda}
\title{Convert Degrees of Freedom to Smoothing Parameter }
\description{
The degree of roughness of an estimated function is controlled by a smoothing parameter lambda that directly multiplies the penalty. However, it can be difficult to interpret or choose this value, and it is often easier to determine the roughness by choosing a value that is equivalent of the degrees of freedom used by the smoothing procedure. This function converts a degrees of freedom value into a multipler lambda. 
}
\usage{
df2lambda(argvals, basisobj, wtvec=rep(1, n), Lfd=0, df=nbasis)
}
\arguments{
  \item{argvals }{Argument values associated with the values to be smoothed. }
  \item{basisobj}{A basis object.}
  \item{wtvec}{A vector of weights for the data to be smoothed. }
  \item{Lfd}{Either a nonnegative integer or a linear differential operator object. }
  \item{df}{The degrees of freedom to be converted. }
}
\details{
The conversion requires a one-dimensional optimization and may be therefore computationally intensive. 
}
\value{
The equivalent smoothing parameter value. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
