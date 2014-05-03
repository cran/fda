\name{TaylorSpline}
\alias{TaylorSpline}
\alias{TaylorSpline.fd}
\alias{TaylorSpline.fdPar}
\alias{TaylorSpline.fdSmooth}
\title{
  Taylor representation of a B-Spline
}
\description{
  Convert B-Spline coefficients into a local Taylor series
  representation expanded about the midpoint between each pair of
  distinct knots.
}
\usage{
TaylorSpline(object, ...)
\method{TaylorSpline}{fd}(object, ...)
\method{TaylorSpline}{fdPar}(object, ...)
\method{TaylorSpline}{fdSmooth}(object, ...)
}
\arguments{
  \item{ object }{
    a spline object possibly of some other class 
  }
  \item{...}{ optional arguments }
}
\details{
  1.  Is \code{object} a spline object with a B-spline basis?  If no,
  throw an error.

  2.  Find \code{knots} and \code{midpoints}.

  3.  Obtain coef(object).

  4.  Determine the number of dimensions of coef(object) and create
  empty \code{coef} and \code{deriv} arrays to match.  Then fill the
  arrays.
}
\value{
  a list with the following components:

  \item{knots}{
    a numeric vector of knots(object, interior=FALSE)
  }
  \item{midpoints}{
    midpoints of intervals defined by unique(knots)
  }
  \item{coef}{
    A matrix of dim = c(nKnots-1, norder) containing the coeffients of
    a polynomial in (x-midpoints[i]) for interval i, where nKnots =
    length(unique(knots)).
  }
  \item{deriv}{
    A matrix of dim = c(nKnots-1, norder) containing the derivatives
    of the spline evaluated at \code{midpoints}.
  }
%    \item{Min, Max}{
%      matrices with one row for each segment giving \code{x} at which
%      local extreme are found, with the corresponding function values.
%      If a local minimum or maximum is not achieved in that interval,
%      then that row of Min or Max is NA.
%
%      If the minimum is achieved at multiple locations in an interval,
%      the smallest such \code{x} value is returned.  Similarly, if the
%      maximum is achieved at multiple locations, the largest such
%      \code{x} is returned.
%    }
}
\author{ Spencer Graves }
\seealso{
  \code{\link{fd}}
  \code{\link{create.bspline.basis}}
}
\examples{
##
## The simplest b-spline basis:  order 1, degree 0, zero interior knots:
##       a single step function
##
bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)
# ... jump to pi to check the code
fd.bspl1.1pi <- fd(pi, basisobj=bspl1.1)
bspl1.1pi <- TaylorSpline(fd.bspl1.1pi)
\dontshow{
bsp1.1ref <- list(knots=0:1, midpoints=0.5,
  coef=array(pi, dim=rep(1, 3), dimnames=list(NULL, 'b0', 'reps 1')),
  deriv=array(pi, dim=rep(1, 3), dimnames=list(NULL, 'D0', 'reps 1')) )
class(bsp1.1ref) <- "Taylor"
stopifnot(all.equal(bspl1.1pi, bsp1.1ref))
}

##
## Cubic spline:  4  basis functions
##
bspl4 <- create.bspline.basis(nbasis=4)
plot(bspl4)
parab4.5 <- fd(c(3, -1, -1, 3)/3, bspl4)
# = 4*(x-.5)
TaylorSpline(parab4.5)

##
## A more realistic example
##
data(titanium)
#  Cubic spline with 5 interior knots (6 segments)
# do not run this example on CRAN to save test time
#if(!CRAN()){
#titan10 <- with(titanium, curfit.free.knot(x, y))
#(titan10T <- TaylorSpline(titan10) )
#}

}
\keyword{ smooth }
\keyword{ manip }
