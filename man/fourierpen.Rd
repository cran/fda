\name{fourierpen}
\alias{fourierpen}
\title{
  Fourier Penalty Matrix
}
\description{
Computes the matrix defining the roughness penalty for functions
expressed in terms of a Fourier basis.
}
\usage{
fourierpen(basisobj, Lfdobj=int2Lfd(2))
}
\arguments{
\item{basisobj}{
a Fourier basis object.
}
\item{Lfdobj}{
either a nonnegative integer or a linear differential operator object.
}
}
\value{
a symmetric matrix of order equal to the number of basis functions
defined by the Fourier basis object.  Each element is the inner product
of two Fourier basis functions after applying the derivative or linear
differential operator defined by Lfdobj.
}
\details{
A roughness penalty for a function $x(t)$ is defined by
integrating the square of either the derivative of  $ x(t) $ or,
more generally, the result of applying a linear differential operator
$L$ to it.  The most common roughness penalty is the integral of
the square of the second derivative, and
this is the default. To apply this roughness penalty, the matrix of
inner products of the basis functions (possibly after applying the
linear differential operator to them) defining this function
is necessary. This function just calls the roughness penalty evaluation
function specific to the basis involved.
}
\seealso{
\code{\link{fourier}}, 
\code{\link{eval.penalty}}, 
\code{\link{getbasispenalty}}
}
\examples{

#  set up a Fourier basis with 13 basis functions
#  and and period 1.0.
basisobj <- create.fourier.basis(c(0,1),13)
#  compute the 13 by 13 matrix of inner products
#  of second derivatives
penmat <- fourierpen(basisobj)

}
% docclass is function
\keyword{smooth}
