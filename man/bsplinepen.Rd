\name{bsplinepen}
\alias{bsplinepen}
\title{
  B-Spline Penalty Matrix
}
\description{
Computes the matrix defining the roughness penalty for functions
expressed in terms of a B-spline basis.
}
\usage{
bsplinepen(basisobj, Lfdobj=2, rng=basisobj$rangeval)
}
\arguments{
\item{basisobj}{
a B-spline basis object.
}
\item{Lfdobj}{
either a nonnegative integer or a linear differential operator object.
}
\item{rng}{
a vector of length 2 defining range over which the basis penalty is to
be computed.
}
}
\value{
a symmetric matrix of order equal to the number of basis functions
defined by the B-spline basis object.  Each element is the inner product
of two B-spline basis functions after applying the derivative or linear
differential operator defined by \code{Lfdobj}.
}
\details{
A roughness penalty for a function $x(t)$ is defined by
integrating the square of either the derivative of $x(t)$ or,
more generally, the result of applying a linear differential operator
$L$ to it.  The most common roughness penalty is the integral of
the square of the second derivative, and
this is the default. To apply this roughness penalty, the matrix of
inner products of the basis functions (possibly after applying the
linear differential operator to them) defining this function
is necessary. This function just calls the roughness penalty evaluation
function specific to the basis involved.
}
\examples{
#  set up a B-spline basis of order 4 with 13 basis functions
#  and knots at 0.0, 0.1,..., 0.9, 1.0.
basisobj <- create.bspline.basis(c(0,1),13)
#  compute the 13 by 13 matrix of inner products of second derivatives
penmat <- bsplinepen(basisobj)
}
% docclass is function
\keyword{smooth}
