\name{wtcheck}
\alias{wtcheck}
\title{
  Check a vector of weights
}
\description{
  Throws an error if \code{wtvec} is not \code{n} positive numbers, and
  return \code{wtvec} (stripped of any \code{dim} atrribute)
}
\usage{
wtcheck(n, wtvec)
}
\arguments{
  \item{n}{ the required length of \code{wtvec} }
  \item{wtvec}{ an object to be checked }
}
\value{
  a vector of \code{n} positive numbers
}
%\seealso{}
\examples{
wtcheck(3, 1:3)

wtcheck(2, matrix(1:2, 2))

}
% docclass is function
\keyword{logic}
