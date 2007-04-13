\name{fourier}
\alias{fourier}
\title{
  Fourier Basis Function Values
}
\description{
  Evaluates a set of Fourier basis functions, or a derivative of these
  functions, at a set of arguments.
}
\usage{
fourier(x, nbasis=n, period=span, nderiv=0)
}
\arguments{
\item{x}{
  a vector of argument values at which the Fourier basis functions are
  to be evaluated.
}
\item{nbasis}{
  the number of basis functions in the Fourier basis.  The first basis
  function is the constant function, followed by sets of  sine/cosine
  pairs.  Normally the number of basis functions will be an odd.  The
  default number is the number of argument values.
}
\item{period}{
  the width of an interval over which all sine/cosine basis functions
  repeat themselves. The default is the difference between the largest
  and smallest argument values.
}
\item{nderiv}{
  the derivative to be evaluated.  The derivative must not exceed the
  order.  The default derivative is 0, meaning that the basis functions
  themselves are evaluated.
}
}
\value{
  a matrix of function values.  The number of rows equals the number of
  arguments, and the number of columns equals the number of basis functions.
}
\seealso{
\code{\link{fourierpen}}
}
\examples{

#  set up a set of 11 argument values
x <- seq(0,1,0.1)
names(x) <- paste("x", 0:10, sep="")
#  compute values for five Fourier basis functions
#  with the default period (1) and derivative (0)
(basismat <- fourier(x, 5))

# Create a false Fourier basis, i.e., nbasis = 1
# = a constant function
fourier(x, 1)

}
% docclass is function
\keyword{smooth}
