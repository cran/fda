\name{monomial}
\alias{monomial}
\title{
  Evaluate Monomial Basis
}
\description{
Computes the values of the powers of argument t.
}
\usage{
monomial(evalarg, exponents, nderiv=0)
}
\arguments{
\item{evalarg}{
a vector of argument values.
}
\item{exponents}{
a vector of nonnegative integer values specifying the
powers to be computed.
}
\item{nderiv}{
a nonnegative integer specifying the order of derivative to be
evaluated.
}
}
\value{
a matrix of values of basis functions.  Rows correspond to
argument values and columns to basis functions.
}
\seealso{
\code{\link{power}}, 
\code{\link{expon}}, 
\code{\link{fourier}}, 
\code{\link{polyg}}, 
\code{\link{bsplineS}}
}
\examples{

# set up a monomial basis for the first five powers
nbasis   <- 5
basisobj <- create.monomial.basis(c(-1,1),nbasis)
#  evaluate the basis
tval <- seq(-1,1,0.1)
basismat <- monomial(tval, 1:basisobj$nbasis)

}
% docclass is function
\keyword{smooth}
