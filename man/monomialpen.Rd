\name{monomialpen}
\alias{monomialpen}
\title{
  Evaluate Monomial Roughness Penalty Matrix
}
\description{
The roughness penalty matrix is the set of
inner products of all pairs of a derivative of integer powers of the
argument.
}
\usage{
monomialpen(basisobj, Lfdobj=int2Lfd(2),
            rng=basisobj$rangeval)
}
\arguments{
\item{basisobj}{
a monomial basis object.
}
\item{Lfdobj}{
either a nonnegative integer specifying an order of derivative
or a linear differential operator object.
}
\item{rng}{
the inner product may be computed over a range that is contained
within the range defined in the basis object.  This is a vector
or length two defining the range.
}
}
\value{
a symmetric matrix of order equal to the number of
monomial basis functions.
}
\seealso{
\code{\link{polynompen}}, 
\code{\link{exponpen}}, 
\code{\link{fourierpen}}, 
\code{\link{bsplinepen}}, 
\code{\link{polygpen}}
}
\examples{

# set up a monomial basis for the first five powers
nbasis   <- 5
basisobj <- create.monomial.basis(c(-1,1),nbasis)
#  evaluate the rougness penalty matrix for the
#  second derivative.
penmat <- monomialpen(basisobj, 2)

}
% docclass is function
\keyword{smooth}
