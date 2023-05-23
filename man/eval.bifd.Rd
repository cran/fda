\name{eval.bifd}
\alias{eval.bifd}
\title{
  Values a Two-argument Functional Data Object
}
\description{
  A vector of argument values for the first argument \code{s} of the
  functional data object to be evaluated.
}
\usage{
eval.bifd(sevalarg, tevalarg, bifd, sLfdobj=0, tLfdobj=0)
}
\arguments{
  \item{sevalarg}{
    a vector of argument values for the first argument \code{s} of the
    functional data object to be evaluated.
  }
  \item{tevalarg}{
    a vector of argument values for the second argument \code{t} of the
    functional data object to be evaluated.
  }
  \item{bifd}{
    a two-argument functional data object.
  }
  \item{sLfdobj}{
    either a nonnegative integer or a linear differential operator
    object.  If present, the derivative or the value of applying the
    operator to the object as a function of the first argument \code{s}
    is evaluated rather than the functions themselves.
  }
  \item{tLfdobj}{
    either a nonnegative integer or a linear differential operator
    object.  If present, the derivative or the value of applying the
    operator to the object as a function of the second argument \code{t}
    is evaluated rather than the functions themselves.
  }
}
\value{
  an array of 2, 3, or 4 dimensions containing the function values.  The
  first dimension corresponds to the argument values in sevalarg, the
  second to argument values in tevalarg, the third if present to
  replications, and the fourth if present to functions.
}
\examples{
# every-other-day basis to save test time
daybasis   <- create.fourier.basis(c(0,365), 183)
harmLcoef  <- c(0,(2*pi/365)^2,0)
harmLfd    <- vec2Lfd(harmLcoef, c(0,365))
templambda <- 1.0
tempfdPar  <- fdPar(daybasis, harmLfd, lambda=1)
tempfd     <- smooth.basis(day.5,
          CanadianWeather$dailyAv[,,"Temperature.C"], tempfdPar)$fd
#    define the variance-covariance bivariate fd object
tempvarbifd <- var.fd(tempfd)
#    evaluate the variance-covariance surface and plot
weektime    <- seq(0,365,len=53)
tempvarmat  <- eval.bifd(weektime,weektime,tempvarbifd)
#    make a perspective plot of the variance function
oldpar <- par(no.readonly=TRUE)
persp(tempvarmat)
par(oldpar)
}
\keyword{smooth}
% Converted by Sd2Rd version 1.21.

