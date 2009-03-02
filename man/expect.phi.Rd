\name{expect.phi}
\alias{normint.phi}
\alias{normden.phi}
\alias{expect.phi}
\alias{expectden.phi}
\alias{expectden.phiphit}
\title{
  Expectation of basis functions
}
\description{
    Computes expectations of basis functions with respect to a density
  by numerical integration using Romberg integration
}
\usage{
normint.phi(basisobj, cvec, JMAX=15, EPS=1e-7)
normden.phi(basisobj, cvec, JMAX=15, EPS=1e-7)
expect.phi(basisobj, cvec, nderiv=0, rng=rangeval,
                     JMAX=15, EPS=1e-7)
expectden.phi(basisobj, cvec, Cval=1, nderiv=0, rng=rangeval,
                     JMAX=15, EPS=1e-7)
expectden.phiphit(basisobj, cvec, Cval=1, nderiv1=0,
                 nderiv2=0, rng=rangeval, JMAX=15, EPS=1e-7)
}
\arguments{
\item{basisobj}{
  a basis function object
}
\item{cvec}{
  coefficient vector defining density, of length NBASIS
}
\item{Cval}{
  normalizing constant defining density
}
\item{nderiv, nderiv1, nderiv2}{
  order of derivative required for basis function expectation

}
\item{rng}{
  a vector of length 2 giving the interval over which the integration is
  to take place
}
\item{JMAX}{
  maximum number of allowable iterations
}
\item{EPS}{
  convergence criterion for relative stop
}
}
\value{
  A vector SS of length NBASIS of integrals of functions.
}
\details{
  normint.phi computes integrals of
      p(x) = exp phi'(x) %*% cvec

   normdel.phi computes integrals of
   p(x) = exp phi"(x) %*% cvec

   expect.phi computes expectations of basis functions with respect to
   intensity
        p(x) <- exp t(c)*phi(x)


   expectden.phi computes expectations of basis functions with respect
to density

      p(x) <- exp(t(c)*phi(x))/Cval

   expectden.phiphit computes expectations of cross product of basis
functions with respect to density

      p(x) <- exp(t(c)*phi(x))/Cval
}
\keyword{smooth}
\seealso{
\code{\link{plot.basisfd}},
}
% docclass is function
\keyword{smooth}

