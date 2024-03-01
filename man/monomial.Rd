\name{monomial}
\alias{monomial}
\title{
  Evaluate Monomial Basis
}
\description{
Computes the values of the powers of argument t.
}
\usage{
monomial(evalarg, exponents=1, nderiv=0, argtrans=c(0,1))
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
\item{argtrans}{
Linearly transform an argument before constructing a basis.
The first element is the shift in value and the second the scale factor.
}
}
\value{
a matrix of values of basis functions.  Rows correspond to
argument values and columns to basis functions.
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
