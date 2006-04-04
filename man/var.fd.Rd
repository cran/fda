\name{var.fd}
\alias{var.fd}
\title{
  Variance and Covariance Surfaces for a Functional Data Object
}
\description{
Compute the variance and covariance functions for
functional data. These are two-argument functions, and therefore define
a surface. If there is one functional data object supplied, its variance
is computed.  If there are two, the covariance between the two is
computed.  Both functional data objects must be univariate.
}
\usage{
var.fd(fdobj1, fdobj2=fdobj1)
}
\arguments{
\item{fdobj1}{
a univariate functional data object.
}
\item{fdobj2}{
a univariate functional data object.
}
}
\value{
a two-argument or bivariate functional data object representing the variance
surface for a single functional data object or the covariance
between two functional data objects.
}
\seealso{
\code{\link{mean.fd}}, 
\code{\link{std.fd}}, 
\code{\link{stddev.fd}}
}
% docclass is function
\keyword{smooth}
