\name{inprod}
\alias{inprod}
\title{
  Inner products of Functional Data Objects.
}
\description{
Computes a matrix of inner products for each pairing of a
replicate for the first argument with a replicate for the
second argument.  This is perhaps the most important function
in the functional data library.  Hardly any analysis fails
to use inner products in some way, and many employ multiple
inner products.  While in certain cases
these may be computed exactly, this is a more general function that
approximates the inner product approximately when required.
The inner product is defined by two derivatives or linear
differential operators that are applied to the
first two arguments.  The range used to compute the inner
product may be contained within the range over which the
functions are defined.  A weight functional data object may
also be used to define weights for the inner product.
}
\usage{
inprod(fdobj1, fdobj2,
       Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0), rng = range1, wtfd = 0)
}
\arguments{
\item{fdobj1}{
a functional data object or a basis object.  If the object is
of the basis class, it is converted to a functional data object
by using the identity matrix as the coefficient matrix.
}
\item{fdobj2}{
a functional data object or a basis object.  If the object is
of the basis class, it is converted to a functional data object
by using the identity matrix as the coefficient matrix.
}
\item{Lfdobj1}{
either a nonnegative integer specifying the derivative of
the first argument to be used, or a linear differential operator
object to be applied to the first argument.
}
\item{Lfdobj2}{
either a nonnegative integer specifying the derivative of
the second argument to be used, or a linear differential operator
object to be applied to the second argument.
}
\item{rng}{
a vector of length 2 defining a restricted range contained
within the range over which the arguments are defined.
}
\item{wtfd}{
a univariate functional data object with a single replicate
defining weights to be used in computing the inner product.
}
}
\value{
a matrix of inner products.  The number of rows is the number of
functions or basis functions in argument \code{fd1}, and the number of
columns is the same thing for argument \code{fd2}.
}
\details{

The approximation method is Richardson extrapolation using numerical
integration by the trapezoidal rule.  At each iteration, the number of
values at which the functions are evaluated is doubled, and a polynomial
extrapolation method is used to estimate the converged integral values
as well as an error tolerance.  Convergence is declared when the
relative error falls below \code{EPS} for all products.  The
extrapolation method generally saves at least one and often two
iterations relative to un-extrapolated trapezoidal integration.
Functional data analyses will seldom need to use \code{inprod}
directly, but code developers should be aware of its pivotal role.
Future work may require more sophisticated and specialized numerical
integration methods.
\code{inprod} computes the definite integral, but some functions
such as \code{smooth.monotone} and \code{register.fd} also need to
compute indefinite integrals.  These use the same approximation scheme,
but usually require more accuracy, and hence more iterations.
When one or both arguments are basis objects, they are converted to
functional data objects using identity matrices as the coefficient
matrices.
\code{inprod} is only called when there is no faster or exact
method available.  In cases where there is, it has been found that the
approximation is good to about four to five significant digits, which is
sufficient for most applications.  Perhaps surprisingly, in the case of
B-splines, the exact method is not appreciably faster, but of course is
more accurate.
\code{inprod} calls function \code{eval.fd} perhaps thousands
of times, so high efficiency for this function and the functions that
it calls is important.
}
\section{References}{
Press, et, al, $Numerical Recipes$.
}
\seealso{
\code{\link{eval.penalty}},
}
% docclass is function
\keyword{smooth}
