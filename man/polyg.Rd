\name{polyg}
\alias{polyg}
\title{
  Polygonal Basis Function Values
}
\description{
Evaluates a set of polygonal basis functions, or a derivative of these
functions, at a set of arguments.
}
\usage{
polyg(x, argvals, nderiv=0)
}
\arguments{
\item{x}{
a vector of argument values at which the polygonal basis functions are to
evaluated.
}
\item{argvals}{
a strictly increasing set of argument values containing the range of x
within it that defines the polygonal basis.  The default is x itself.
}
\item{nderiv}{
the order of derivative to be evaluated.  The derivative must not exceed
one.  The default derivative is 0, meaning that the basis functions
themselves are evaluated.
}
}
\value{
a matrix of function values.  The number of rows equals the number of
arguments, and the number of columns equals the number of basis
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
  \code{\link{create.polygonal.basis}}, 
  \code{\link{polygpen}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
#  set up a set of 21 argument values
x <- seq(0,1,0.05)
#  set up a set of 11 argument values
argvals <- seq(0,1,0.1)
#  with the default period (1) and derivative (0)
basismat <- polyg(x, argvals)
#  plot the basis functions
matplot(x, basismat, type="l")
par(oldpar)
}
% docclass is function
\keyword{smooth}
