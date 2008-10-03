\name{polygpen}
\alias{polygpen}
\title{
  Polygonal Penalty Matrix
}
\description{
Computes the matrix defining the roughness penalty for functions
expressed in terms of a polygonal basis.
}
\usage{
polygpen(basisobj, Lfdobj=int2Lfd(1))
}
\arguments{
\item{basisobj}{
a polygonal functional basis object.
}
\item{Lfdobj}{
either an integer that is either 0 or 1, or a
linear differential operator object of degree 0 or 1.
}
}
\value{
a symmetric matrix of order equal to the number of basis functions
defined by the polygonal basis object.  Each element is the inner product
of two polygonal basis functions after applying the derivative or linear
differential operator defined by Lfdobj.
}
\details{
a roughness penalty for a function $ x(t) $ is defined by
integrating the square of either the derivative of  $ x(t) $ or,
more generally, the result of applying a linear differential operator
$ L $ to it.  The only roughness penalty possible aside from
penalizing the size of the function itself is the integral
of the square of the first derivative, and
this is the default. To apply this roughness penalty, the matrix of
inner products produced by this function is necessary.
}
\seealso{
\code{\link{create.polygonal.basis}},
\code{\link{polyg}}
}
\examples{

#  set up a sequence of 11 argument values
argvals <- seq(0,1,0.1)
#  set up the polygonal basis
basisobj <- create.polygonal.basis(argvals)
#  compute the 11 by 11 penalty matrix

penmat <- polygpen(basisobj)

}
% docclass is function
\keyword{smooth}
