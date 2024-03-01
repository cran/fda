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
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\examples{
wtcheck(3, 1:3)
wtcheck(2, matrix(1:2, 2))
}
% docclass is function
\keyword{logic}
