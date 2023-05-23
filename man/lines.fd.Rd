\name{lines.fd}
\alias{lines.fd}
\alias{lines.fdSmooth}
\title{
  Add Lines from Functional Data to a Plot
}
\description{
Lines defined by functional observations are added to an existing plot.
}
\usage{
\method{lines}{fd}(x, Lfdobj=int2Lfd(0), nx=201, ...)
\method{lines}{fdSmooth}(x, Lfdobj=int2Lfd(0), nx=201, ...)
}
\arguments{
  \item{x}{
    a univariate functional data object to be evaluated at \code{nx}
    points over \code{xlim} and added as a line to an existing plot.  
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object.  If present, the derivative or the value of applying the
    operator is evaluated rather than the functions themselves.
  }
  \item{nx}{
    Number of points within \code{xlim} at which to evaluate \code{x}
    for plotting.  
  }
  \item{\dots }{
    additional arguments such as axis titles and so forth that can be
    used in plotting programs called by \code{lines.fd} or
    \code{lines.fdSmooth}.  
  }
}
\section{Side Effects}{
Lines added to an existing plot.
}
\seealso{
  \code{\link{plot.fd}}, 
  \code{\link{plotfit.fd}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
##
## plot a fit with 3 levels of smoothing
##
x <- seq(-1,1,0.02)
y <- x + 3*exp(-6*x^2) + sin(1:101)/2
# sin not rnorm to make it easier to compare
# results across platforms 

result4.0  <- smooth.basisPar(argvals=x, y=y, lambda=1)
result4.m4 <- smooth.basisPar(argvals=x, y=y, lambda=1e-4)

plot(x, y)
lines(result4.0$fd)
lines(result4.m4$fd,      col='blue')
lines.fdSmooth(result4.0, col='red') 

plot(x, y, xlim=c(0.5, 1))
lines.fdSmooth(result4.0)
lines.fdSmooth(result4.m4, col='blue')
# no visible difference from the default?  
par(oldpar)
}
\keyword{smooth}
