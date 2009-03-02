\name{eval.fd}
\alias{eval.fd}
\alias{predict.fd}
\alias{predict.fdPar}
\alias{predict.fdSmooth}
\alias{fitted.fdSmooth}
\alias{residuals.fdSmooth}
\title{
  Values of a Functional Data Object
}
\description{
  Evaluate a functional data object at specified argument values, or
  evaluate a derivative or the result of applying a linear differential
  operator to the functional object.
}
\usage{
eval.fd(evalarg, fdobj, Lfdobj=0)
\method{predict}{fd}(object, newdata=NULL, Lfdobj=0, ...)
\method{predict}{fdPar}(object, newdata=NULL, Lfdobj=0, ...)
\method{predict}{fdSmooth}(object, newdata=NULL, Lfdobj=0, ...)
\method{fitted}{fdSmooth}(object, ...)
\method{residuals}{fdSmooth}(object, ...)
}
\arguments{
  \item{evalarg, newdata}{
    a vector of argument values at which the functional data object is
    to be evaluated.
  }
  \item{fdobj}{
    a functional data object to be evaluated.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object.  If present, the derivative or the value of applying the
    operator is evaluated rather than the functions themselves.
  }
  \item{object}{
    an object of class \code{fd}
  }
  \item{\dots}{
    optional arguments for \code{predict}, not currently used
  }
}
\details{
  \code{eval.fd} evaluates \code{Lfdobj} of \code{fdobj} at
  \code{evalarg}.  \code{predict.fd} is a convenience wrapper for
  \code{eval.fd}.  If \code{newdata} is NULL and
  \code{fdobj[['basis']][['type']]} is \code{bspline}, \code{newdata} =
  \code{unique(knots(fdojb,interior=FALSE))};  otherwise, \code{newdata}
  = \code{fdobj[['basis']][['rangeval']]}.

  \code{predict.fdSmooth}, \code{fitted.fdSmooth} and
  \code{residuals.fdSmooth} are other wrapper for \code{eval.fd}.
}
\value{
  an array of 2 or 3 dimensions containing the function
  values.  The first dimension corresponds to the argument values in
  \code{evalarg},
  the second to replications, and the third if present to functions.
}
\author{
  Soren Hosgaard wrote an initial version of \code{predict.fdSmooth},
  \code{fitted.fdSmooth}, and \code{residuals.fdSmooth}.
}
\seealso{
  \code{\link{getbasismatrix}},
  \code{\link{eval.bifd}},
  \code{\link{eval.penalty}},
  \code{\link{eval.monfd}},
  \code{\link{eval.posfd}}
}
\examples{
##
## eval.fd
##
#    set up the fourier basis
daybasis <- create.fourier.basis(c(0, 365), nbasis=65)
#  Make temperature fd object
#  Temperature data are in 12 by 365 matrix tempav
#  See analyses of weather data.
#  Set up sampling points at mid days
#  Convert the data to a functional data object
tempfd <- data2fd(CanadianWeather$dailyAv[,,"Temperature.C"],
                   day.5, daybasis)
#   set up the harmonic acceleration operator
Lbasis  <- create.constant.basis(c(0, 365))
Lcoef   <- matrix(c(0,(2*pi/365)^2,0),1,3)
bfdobj  <- fd(Lcoef,Lbasis)
bwtlist <- fd2list(bfdobj)
harmaccelLfd <- Lfd(3, bwtlist)
#   evaluate the value of the harmonic acceleration
#   operator at the sampling points
Ltempmat <- eval.fd(day.5, tempfd, harmaccelLfd)
#  Plot the values of this operator
matplot(day.5, Ltempmat, type="l")

##
## predict.fd
##
predict(tempfd) # end points only at 35 locations
str(predict(tempfd, day.5)) # 365 x 35 matrix
str(predict(tempfd, day.5, harmaccelLfd))

# cublic splie with knots at 0, .5, 1
bspl3 <- create.bspline.basis(c(0, .5, 1))
plot(bspl3) # 5 bases
fd.bspl3 <- fd(c(0, 0, 1, 0, 0), bspl3)
pred3 <- predict(fd.bspl3)

pred3. <- matrix(c(0, .5, 0), 3)
dimnames(pred3.) <- list(NULL, 'reps 1')
\dontshow{stopifnot(}
all.equal(pred3, pred3.)
\dontshow{)}

pred.2 <- predict(fd.bspl3, c(.2, .8))

pred.2. <- matrix(.176, 2, 1)
dimnames(pred.2.) <- list(NULL, 'reps 1')
\dontshow{stopifnot(}
all.equal(pred.2, pred.2.)
\dontshow{)}

##
## predict.fdSmooth
##
lipSm9 <- smooth.basisPar(liptime, lip, lambda=1e-9)
plot(lipSm9)

pred.lipSm9 <- predict(lipSm9)
\dontshow{stopifnot(}
all.equal(pred.lipSm9, fitted(lipSm9))
\dontshow{)}

resid.lipSm9 <- resid(lipSm9)
\dontshow{stopifnot(}
all.equal(resid.lipSm9, residuals(lipSm9))
\dontshow{)}

}
% docclass is function
\keyword{smooth}
