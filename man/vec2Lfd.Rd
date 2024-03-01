\name{vec2Lfd}
\alias{vec2Lfd}
\title{
  Make a Linear Differential Operator Object from a Vector
}
\description{
A linear differential operator object of order $m$ is
constructed from the number in a vector of length $m$.
}
\usage{
vec2Lfd(bwtvec, rangeval=c(0,1))
}
\arguments{
\item{bwtvec}{
a vector of coefficients to define the linear differential
operator object
}
\item{rangeval}{
a vector of length 2 specifying the range over which the
operator is defined
}
}
\value{
  a linear differential operator object
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{int2Lfd}},
  \code{\link{Lfd}}
}
\examples{
#  define the harmonic acceleration operator used in the
#  analysis of the daily temperature data
Lcoef <- c(0,(2*pi/365)^2,0)
harmaccelLfd <- vec2Lfd(Lcoef, c(0,365))

hmat <- vec2Lfd(matrix(Lcoef, 1), c(0, 365))

\dontshow{stopifnot(}
all.equal(harmaccelLfd, hmat)
\dontshow{)}

}
% docclass is function
\keyword{smooth}
