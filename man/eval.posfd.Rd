\name{eval.posfd}
\alias{eval.posfd}
\alias{predict.posfd}
\alias{fitted.posfd}
\alias{residuals.posfd}
\title{
  Evaluate a Positive Functional Data Object
}
\description{
  Evaluate a positive functional data object at specified argument
  values, or evaluate a derivative of the functional object.
}
\usage{
eval.posfd(evalarg, Wfdobj, Lfdobj=int2Lfd(0), returnMatrix=FALSE)
\method{predict}{posfd}(object, newdata=NULL, Lfdobj=0,
                       returnMatrix=FALSE, ...)
\method{fitted}{posfd}(object, ...)
\method{residuals}{posfd}(object, ...)
}
\arguments{
  \item{evalarg, newdata}{
    a vector of argument values at which the functional data object is
    to be evaluated.
  }
  \item{Wfdobj}{
    a functional data object that defines the positive function to be
    evaluated.  Only univariate functions are permitted.
  }
  \item{Lfdobj}{
    a nonnegative integer specifying a derivative to be evaluated.  At
    this time of writing, permissible derivative values are 0, 1 or 2.
    A linear differential operator is not allowed.
  }
  \item{object}{
    an object of class \code{posfd} that defines the positive function
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
  A positive function data object $h(t)$ is defined by $h(t) =[exp
  Wfd](t)$.  The function \code{Wfdobj} that defines the positive
  function is usually estimated by positive smoothing function
  \code{smooth.pos}
}
\value{
  a matrix containing the positive function values.  The first dimension
  corresponds to the argument values in \code{evalarg} and the second to
  replications.
}
\seealso{
  \code{\link{eval.fd}},
  \code{\link{eval.monfd}}
}
\examples{
harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
smallbasis   <- create.fourier.basis(c(0, 365), 65)
index        <- (1:35)[CanadianWeather$place == "Vancouver"]
VanPrec      <- CanadianWeather$dailyAv[,index, "Precipitation.mm"]
lambda       <- 1e4
dayfdPar     <- fdPar(smallbasis, harmaccelLfd, lambda)
VanPrecPos   <- smooth.pos(day.5, VanPrec, dayfdPar)
#  compute fitted values using eval.posfd()
VanPrecPosFit1 <- eval.posfd(day.5, VanPrecPos$Wfdobj)
#  compute fitted values using predict()
VanPrecPosFit2 <- predict(VanPrecPos, day.5)
\dontshow{stopifnot(}
all.equal(VanPrecPosFit1, VanPrecPosFit2)
\dontshow{)}
#  compute fitted values using fitted()
VanPrecPosFit3 <- fitted(VanPrecPos)
#  compute residuals
VanPrecRes <- resid(VanPrecPos)
\dontshow{stopifnot(}
all.equal(VanPrecRes, VanPrecPos$y-VanPrecPosFit3)
\dontshow{)}

}
% docclass is function
\keyword{smooth}
