\name{eval.monfd}
\alias{eval.monfd}
\alias{predict.monfd}
\alias{fitted.monfd}
\alias{residuals.monfd}
\title{
  Values of a Monotone Functional Data Object
}
\description{
  Evaluate a monotone functional data object at specified argument
  values, or evaluate a derivative of the functional object.
}
\usage{
eval.monfd(evalarg, Wfdobj, Lfdobj=int2Lfd(0), returnMatrix=FALSE)
\method{predict}{monfd}(object, newdata=NULL, Lfdobj=0, returnMatrix=FALSE, ...)
\method{fitted}{monfd}(object, ...)
\method{residuals}{monfd}(object, ...)
}
\arguments{
  \item{evalarg, newdata}{
    a vector of argument values at which the functional data object is
    to be evaluated.
  }
  \item{Wfdobj}{
    an object of class \code{fd} that defines the monotone function
    to be evaluated.  Only univariate functions are permitted.
  }
  \item{Lfdobj}{
    a nonnegative integer specifying a derivative to be evaluated.  At
    this time of writing, permissible derivative values are 0, 1, 2, or
    3.  A linear differential operator is not allowed.
  }
  \item{object}{
    an object of class \code{monfd} that defines the monotone function
    to be evaluated.  Only univariate functions are permitted.
  }
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
  \item{\dots}{
    optional arguments required by \code{predict};  not currently used.
  }
}
\details{
  A monotone function data object $h(t)$ is defined by $h(t) = [D^\{-1\}
  exp Wfdobj](t)$.  In this equation, the operator  $D^\{-1\}$ means
  taking the indefinite integral of the function to which it applies.
  Note that this equation implies that the monotone function has a value
  of zero at the lower limit of the arguments.  To actually fit monotone
  data, it will usually be necessary to estimate an intercept and a
  regression coefficient to be applied to $h(t)$, usually with the least
  squares regression function \code{lsfit}.  The function \code{Wfdobj}
  that defines the monotone function is usually estimated by monotone
  smoothing function \code{smooth.monotone.}

  \code{eval.monfd} only computes the standardized monotone form.
  \code{predict.monfd} computes the scaled version using
  \code{with(object, beta[1] + beta[2]*eval.monfd(...))} if Lfdobj = 0
  or beta[2]*eval.monfd(...) if Lfdobj > 0.
}
\value{
  a matrix containing the monotone function values.  The first dimension
  corresponds to the argument values in \code{evalarg} and the second to
  replications.
}
\seealso{
  \code{\link{eval.fd}},
  \code{\link{smooth.monotone}}
  \code{\link{eval.posfd}}
}
\examples{

#  Estimate the acceleration functions for growth curves
#  See the analyses of the growth data.
#  Set up the ages of height measurements for Berkeley data
age <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))
#  Range of observations
rng <- c(1,18)
#  First set up a basis for monotone smooth
#  We use b-spline basis functions of order 6
#  Knots are positioned at the ages of observation.
norder <- 6
nage   <- length(age)
nbasis <- nage + norder - 2
wbasis <- create.bspline.basis(rng, nbasis, norder, age)
#  starting values for coefficient
cvec0 <- matrix(0,nbasis,1)
Wfd0  <- fd(cvec0, wbasis)
#  set up functional parameter object
Lfdobj    <- 3          #  penalize curvature of acceleration
lambda    <- 10^(-0.5)  #  smoothing parameter
growfdPar <- fdPar(Wfd0, Lfdobj, lambda)
#  Smooth the data for the first girl
hgt1 <- growth$hgtf[,1]
#   set conv = 0.1 and iterlim=1 to reduce the compute time
#   required for this test on CRAN;
#   We would not do this normally.
result <- smooth.monotone(age, hgt1, growfdPar, conv=0.1,
                          iterlim=1)
#  Extract the functional data object and regression
#  coefficients
Wfd  <- result$Wfdobj
beta <- result$beta
#  Evaluate the fitted height curve over a fine mesh
agefine <- seq(1,18,len=60)
hgtfine <- beta[1] + beta[2]*eval.monfd(agefine, Wfd)
#  Plot the data and the curve
plot(age, hgt1, type="p")
lines(agefine, hgtfine)
#  Evaluate the acceleration curve
accfine <- beta[2]*eval.monfd(agefine, Wfd, 2)
#  Plot the acceleration curve
plot(agefine, accfine, type="l")
lines(c(1,18),c(0,0),lty=4)

##
## using predict.monfd
##
hgtfit <- with(result, beta[1]+beta[2]*eval.monfd(argvals, Wfdobj))
hgtfit. <- fitted(result)
\dontshow{stopifnot(}
all.equal(hgtfit, hgtfit.)
\dontshow{)}

accfine. <- predict(result, agefine, Lfdobj=2)
\dontshow{stopifnot(}
all.equal(accfine, accfine.)
\dontshow{)}

growthResid <- resid(result)
\dontshow{stopifnot(}
all.equal(growthResid, with(result, y-hgtfit.))
\dontshow{)}

}
% docclass is function
\keyword{smooth}
