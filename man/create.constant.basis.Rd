\name{create.constant.basis}
\alias{create.constant.basis}
\title{
  Create a Constant Basis
}
\description{
  Create a constant basis object, defining a single basis function
  whose value is everywhere 1.0.
}
\usage{
create.constant.basis(rangeval=c(0, 1), names="const", axes=NULL)
}
\arguments{
  \item{rangeval}{
    a vector of length 2 containing the initial and final
    values of argument t defining the interval over which the functional
    data object can be evaluated.  However, this is seldom used
    since the value of the basis function does not depend on the range
    or any argument values.
  }
  \item{names}{
    a character vector of length 1.
  }
  \item{axes}{
    an optional list used by selected \code{plot} functions to create
    custom \code{axes}.  If this \code{axes} argument is not
    \code{NULL}, functions \code{plot.basisfd}, \code{plot.fd},
    \code{plot.fdSmooth} \code{plotfit.fd}, \code{plotfit.fdSmooth}, and
    \code{plot.Lfd} will create axes via \code{do.call(x$axes[[1]],
      x$axes[-1])}.  The primary example of this uses
    \code{list("axesIntervals", ...)}, e.g., with \code{Fourier} bases
    to create \code{CanadianWeather} plots
  }
}
\value{
  a basis object with type component \code{const}.
}
\seealso{
  \code{\link{basisfd}},
  \code{\link{create.bspline.basis}},
  \code{\link{create.exponential.basis}},
  \code{\link{create.fourier.basis}},
  \code{\link{create.monomial.basis}},
  \code{\link{create.polygonal.basis}},
  \code{\link{create.power.basis}}
}
\examples{

basisobj <- create.constant.basis(c(-1,1))

}
% docclass is function
\keyword{smooth}
