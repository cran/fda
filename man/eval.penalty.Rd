\name{eval.penalty}
\alias{eval.penalty}
\title{
  Evaluate a Basis Penalty Matrix
}
\description{
A basis roughness penalty matrix is the matrix containing
the possible inner products of pairs of basis functions.
These inner products are typically defined in terms of
the value of a derivative or of a linear differential
operator applied to the basis function.  The basis penalty
matrix plays an important role in the computation of
functions whose roughness is controlled by a roughness
penalty.
}
\usage{
eval.penalty(basisobj, Lfdobj=int2Lfd(0), rng=rangeval)
}
\arguments{
  \item{basisobj}{a basis object.}
  \item{Lfdobj}{
    either a nonnegative integer defining an order of a
    derivative or a linear differential operator.
  }
  \item{rng}{
    a vector of length 2 defining a restricted range.
    Optionally, the inner products can be computed over
    a range of argument values that lies within the
    interval covered by the basis function definition.
  }
}
\value{
  a square symmetric matrix whose order is equal
  to the number of basis functions defined by
  the basis function object \code{ basisobj }.
  If \code{Lfdobj} is $m$ or a linear
  differential operator of order $m$, the rank
  of the matrix should be at least approximately equal to
  its order minus  $m$.
}
\details{

 The inner product can be computed exactly for many
 types of bases if $m$ is an integer.  These include
 B-spline, fourier, exponential, monomial, polynomial and power bases.
 In other cases, and for noninteger operators, the
 inner products are computed by an iterative numerical
 integration method called Richard extrapolation using the
 trapezoidal rule.

 If the penalty matrix must be evaluated repeatedly,
 computation can be greatly speeded up by avoiding the use
 of this function, and instead using quadrature points and
 weights defined by Simpson's rule.  
}
\seealso{
  \code{\link{getbasispenalty}}, 
  \code{\link{eval.basis}}, 
}
% docclass is function
\keyword{smooth}
