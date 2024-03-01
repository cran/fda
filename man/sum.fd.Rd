\name{sum.fd}
\alias{sum.fd}
\title{
  Sum of Functional Data
}
\description{
Evaluate the sum of a set of functions in a functional data object.
}
\usage{
\method{sum}{fd}(..., na.rm)
}
\arguments{
\item{...}{
a functional data object to sum.
}
\item{na.rm}{Not used.}
}
\value{
a functional data object with a single replication
that contains the sum of the functions in the object \code{fd}.
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
  \code{\link{mean.fd}},
  \code{\link{std.fd}},
  \code{\link{stddev.fd}},
  \code{\link{center.fd}}
}
% docclass is function
\keyword{smooth}
