\name{powerbasis}
\alias{powerbasis}
\title{
  Power Basis Function Values
}
\description{
Evaluates a set of power basis functions, or a derivative of these
functions, at a set of arguments.
}
\usage{
powerbasis(x, exponents, nderiv=0)
}
\arguments{
\item{x}{
a vector of argument values at which the power basis functions are to
evaluated. Since exponents may be negative, for example after
differentiation, it is required that all argument values be positive.
}
\item{exponents}{
a vector of exponents defining the power basis functions.  If
$y$ is such a rate value, the corresponding basis function is
$x$ to the power $y$.  The number of basis functions is equal to the
number of exponents.
}
\item{nderiv}{
the derivative to be evaluated.  The derivative must not exceed the
order.  The default derivative is 0, meaning that the basis functions
themselves are evaluated.
}
}
\value{
a matrix of function values.  The number of rows equals the number of
arguments, and the number of columns equals the number of basis
functions.
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
  \code{\link{create.power.basis}}, 
  \code{\link{powerpen}}
}
\examples{

#  set up a set of 10 positive argument values.
x <- seq(0.1,1,0.1)
#  compute values for three power basis functions
exponents <- c(0, 1, 2)
#  evaluate the basis matrix
basismat <- powerbasis(x, exponents)

}
% docclass is function
\keyword{smooth}
