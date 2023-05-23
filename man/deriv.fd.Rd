\name{deriv.fd}
\alias{deriv.fd}
\title{
  Compute a Derivative of a Functional Data Object
}
\description{
  A derivative of a functional data object, or the result of applying
  a linear differential operator to a functional data object, is then
  converted to a functional data object. This is intended for situations
  where a derivative is to be manipulated as a functional data object
  rather than simply evaluated.
}
\usage{
\method{deriv}{fd}(expr, Lfdobj=int2Lfd(1), ...)
}
\arguments{
  \item{expr}{
    a functional data object.  It is assumed that the basis for
    representing the object can support the order of derivative
    to be computed.  For B-spline bases, this means that the
    order of the spline must be at least one larger than the order of
    the derivative to be computed.
  }
  \item{Lfdobj}{
    either a positive integer or a linear differential operator object.
  }
  \item{...}{Other arguments to match generic for 'deriv'}
}
\value{
  a functional data object for the derivative
}
\details{
  Typically, a derivative has more high frequency variation or detail
  than the function itself.  The basis defining the function is used,
  and therefore this must have enough basis functions to represent
  the variation in the derivative satisfactorily.
}
\seealso{
  \code{\link{getbasismatrix}},
  \code{\link{eval.basis}}
  \code{\link{deriv}}
}
\examples{

#  Estimate the acceleration functions for growth curves
#  See the analyses of the growth data.
#  Set up the ages of height measurements for Berkeley data
age <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))
#  Range of observations
rng <- c(1,18)
#  Set up a B-spline basis of order 6 with knots at ages
knots  <- age
norder <- 6
nbasis <- length(knots) + norder - 2
hgtbasis <- create.bspline.basis(rng, nbasis, norder, knots)
#  Set up a functional parameter object for estimating
#  growth curves.  The 4th derivative is penalyzed to
#  ensure a smooth 2nd derivative or acceleration.
Lfdobj <- 4
lambda <- 10^(-0.5)   #  This value known in advance.
growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)
#  Smooth the data.  The data for the boys and girls
#  are in matrices hgtm and hgtf, respectively.
hgtmfd <- smooth.basis(age, growth$hgtm, growfdPar)$fd
hgtffd <- smooth.basis(age, growth$hgtf, growfdPar)$fd
#  Compute the acceleration functions
accmfd <- deriv.fd(hgtmfd, 2)
accffd <- deriv.fd(hgtffd, 2)
#  Plot the two sets of curves
oldpar <- par(mfrow=c(2,1))
plot(accmfd)
plot(accffd)
par(oldpar)
}
% docclass is function
\keyword{smooth}
