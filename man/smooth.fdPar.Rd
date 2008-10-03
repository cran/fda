\name{smooth.fdPar}
\alias{smooth.fdPar}
\title{
  Smooth a functional data object using a directly specified roughness
  penalty
}
\description{
  Smooth data already converted to a functional data object, fdobj,
  using directly specified criteria.
}
\usage{
smooth.fdPar(fdobj, Lfdobj=NULL, lambda=0,
             estimate=TRUE, penmat=NULL)
}
\arguments{
  \item{fdobj}{
    a functional data object to be smoothed.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object.

    If \code{NULL}, Lfdobj depends on fdobj[['basis']][['type']]:

    \itemize{
      \item{bspline}{
	Lfdobj <- int2Lfd(max(0, norder-2)), where norder =
	norder(fdobj).
      }
      \item{fourier}{
	Lfdobj = a harmonic acceleration operator:

	\code{Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng)}

	where rng = fdobj[['basis']][['rangeval']].
      }
      \item{anything else}{Lfdobj <- int2Lfd(0)}
    }
  }
  \item{lambda}{
    a nonnegative real number specifying the amount of smoothing
    to be applied to the estimated functional parameter.
  }
  \item{estimate}{
    a logical value:  if \code{TRUE}, the functional parameter is
    estimated, otherwise, it is held fixed.
  }
  \item{penmat}{
    a roughness penalty matrix.  Including this can eliminate the need
    to compute this matrix over and over again in some types of
    calculations.
  }
}
\details{
  1.  fdPar

  2.  smooth.fd
}
\value{
  a functional data object.
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\seealso{
\code{\link{smooth.fd}},
\code{\link{fdPar}},
\code{\link{smooth.basis}},
\code{\link{smooth.pos}},
\code{\link{smooth.morph}}
}
\examples{
##
## 1.  B-spline example
##
#  Shows the effects of two levels of smoothing
#  where the size of the third derivative is penalized.
#  The null space contains quadratic functions.
x <- seq(-1,1,0.02)
y <- x + 3*exp(-6*x^2) + rnorm(rep(1,101))*0.2
#  set up a saturated B-spline basis
basisobj <- create.bspline.basis(c(-1,1),81)
#  convert to a functional data object that interpolates the data.
result <- smooth.basis(x, y, basisobj)
yfd  <- result$fd
#  set up a functional parameter object with smoothing
#  parameter 1e-6 and a penalty on the 2nd derivative.
yfdPar <- fdPar(yfd, 2, 1e-6)
yfd1 <- smooth.fd(yfd, yfdPar)

yfd1. <- smooth.fdPar(yfd, 2, 1e-6)
all.equal(yfd1, yfd1.)
# TRUE

#  set up a functional parameter object with smoothing
#  parameter 1 and a penalty on the 2nd derivative.
yfd2 <- smooth.fdPar(yfd, 2, 1)

#  plot the data and smooth
plot(x,y)           # plot the data
lines(yfd1, lty=1)  #  add moderately penalized smooth
lines(yfd2, lty=3)  #  add heavily  penalized smooth
legend(-1,3,c("0.000001","1"),lty=c(1,3))
#  plot the data and smoothing using function plotfit.fd
plotfit.fd(y, x, yfd1)  # plot data and smooth
##
## 2.  Fourier basis with harmonic acceleration operator
##
daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65)
daySmooth <- with(CanadianWeather, smooth.basis(day.5,
       dailyAv[,,"Temperature.C"],
       daybasis65, fdnames=list("Day", "Station", "Deg C")) )

daySmooth2 <- smooth.fdPar(daySmooth$fd, lambda=1e5)
op <- par(mfrow=c(2,1))
plot(daySmooth)
plot(daySmooth2)
par(op)
}
% docclass is function
\keyword{smooth}
