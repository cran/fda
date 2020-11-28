\name{as.fd}
\alias{as.fd}
\alias{as.fd.fdSmooth}
\alias{as.fd.function}
\alias{as.fd.smooth.spline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert a spline object to class 'fd'
}
\description{
  Translate a spline object of another class into the 
  Functional Data (class \code{fd}) format.
}
\usage{
as.fd(x, ...)
\method{as.fd}{fdSmooth}(x, ...)
\method{as.fd}{function}(x, ...)
\method{as.fd}{smooth.spline}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an object to be converted to class \code{fd}.
  }
  \item{\dots}{
    optional arguments passed to specific methods, currently
    unused.
  }
}
\details{
  The behavior depends on the \code{class} and nature of \code{x}.

  \itemize{
    \item{as.fd.fdSmooth}{extract the \code{fd} component}

    \item{as.fd.function}{
      Create an \code{fd} object from a function of the form
      created by \code{splinefun}.  This will translate method
      = 'fmn' and 'natural' but not 'periodic':  'fmn' splines
      are isomorphic to standard B-splines with coincident
      boundary knots, which is the basis produced by 
      \code{create.bspline.basis}.  'natural' splines occupy 
      a subspace of this space, with the restriction that the
      second derivative at the end points is zero (as noted in
      the Wikipedia \code{spline} article).  'periodic' splines
      do not use coincident boundary knots and are not 
      currently supported in \code{fda};  instead, \code{fda}
      uses finite Fourier bases for periodic phenomena.
    }

    \item{as.fd.smooth.spline}{
      Create an \code{fd} object from a \code{smooth.spline}
      object.
    }
  }
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), 
  \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), 
  \emph{Applied Functional Data Analysis}, Springer, New York.

  \code{spline} entry in \emph{Wikipedia}
  \url{https://en.wikipedia.org/wiki/Spline_(mathematics)}
}
\author{Spencer Graves}
\seealso{
  \code{\link{fd}}
  \code{\link{splinefun}}
}
\examples{
##
## as.fd.fdSmooth
##
girlGrowthSm <- with(growth, 
  smooth.basisPar(argvals=age, y=hgtf, lambda=0.1))
girlGrowth.fd <- as.fd(girlGrowthSm)

##
## as.fd.function(splinefun(...), ...)
##
x2 <- 1:7
y2 <- sin((x2-0.5)*pi)
f <- splinefun(x2, y2)
fd. <- as.fd(f)
x. <- seq(1, 7, .02)
fx. <- f(x.)
fdx. <- eval.fd(x., fd.)

# range(y2, fx., fdx.) generates an error 2012.04.22

rfdx <- range(fdx.)

plot(range(x2), range(y2, fx., rfdx), type='n')
points(x2, y2)
lines(x., sin((x.-0.5)*pi), lty='dashed')
lines(x., f(x.), col='blue')
lines(x., eval.fd(x., fd.), col='red', lwd=3, lty='dashed')
# splinefun and as.fd(splineful(...)) are close
# but quite different from the actual function
# apart from the actual 7 points fitted,
# which are fitted exactly
# ... and there is no information in the data
# to support a better fit!

# Translate also a natural spline
fn <- splinefun(x2, y2, method='natural')
fn. <- as.fd(fn)
lines(x., fn(x.), lty='dotted', col='blue')
lines(x., eval.fd(x., fn.), col='green', lty='dotted', lwd=3)

\dontrun{
# Will NOT translate a periodic spline
fp <- splinefun(x, y, method='periodic')
as.fd(fp)
#Error in as.fd.function(fp) :
#  x (fp)  uses periodic B-splines, and as.fd is programmed
#   to translate only B-splines with coincident boundary knots.
}

##
## as.fd.smooth.spline
##
cars.spl <- with(cars, smooth.spline(speed, dist))
cars.fd <- as.fd(cars.spl)

plot(dist~speed, cars)
lines(cars.spl)
sp. <- with(cars, seq(min(speed), max(speed), len=101))
d. <- eval.fd(sp., cars.fd)
lines(sp., d., lty=2, col='red', lwd=3)
}
\keyword{smooth}
\keyword{manip}
