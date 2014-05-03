\name{landmarkreg}
\alias{landmarkreg}
\title{
  Landmark Registration of Functional Observations
}
\description{
It is common to see that among a set of functions certain prominent
features such peaks and valleys, called $landmarks$, do not occur
at the same times, or other
argument values.  This is called $phase variation$, and it can be
essential to align these features before proceeding with further
functional data analyses.  This function uses the timings of these
features to align or register the curves.  The registration involves
estimating a nonlinear transformation of the argument continuum for each
functional observation.  This transformation is called a warping
function. It must be strictly increasing and smooth.
}
\usage{
landmarkreg(fdobj, ximarks, x0marks=xmeanmarks,
            WfdPar, monwrd=FALSE, ylambda=1e-10, returnMatrix=FALSE)
}
\arguments{
  \item{fdobj}{
    a functional data object containing the curves to be registered.
  }
  \item{ximarks}{
    a matrix containing the timings or argument values associated with
    the landmarks for the observations in \code{fd} to be registered.
    The number of rows N equals the number of observations, and the
    number of columns NL equals the number of landmarks. These landmark
    times must be in the interior of the interval over which the
    functions are defined.
  }
  \item{x0marks}{
    a vector of length NL of times of landmarks for target curve.  If
    not supplied, the mean of the landmark times in \code{ximarks} is
    used.
  }
  \item{WfdPar}{
    a functional parameter object defining the warping functions that
    transform time in order to register the curves.
  }
  \item{monwrd}{
    A logical value:  if \code{TRUE}, the warping function is estimated
    using a monotone smoothing methhod; otherwise, a regular smoothing
    method is used, which is not guaranteed to give strictly monotonic
    warping functions.
  }
  \item{ylambda}{
    Smoothing parameter controlling the smoothness of the registered
    functions.  It can happen with high dimensional bases that local
    wiggles can appear in the registered curves or their derivatives
    that are not seen in the unregistered versions.  In this case,
    this parameter should be increased to the point where they
    disappear.
  }
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
}
\value{
  a named list of length 2 with components:

  \item{fdreg}{
    a functional data object for the registered curves.
  }
  \item{warpfd}{
    a functional data object for the warping functions.
  }
}
\details{
  It is essential that the location of every landmark be clearly defined
  in each of the curves as well as the template function.  If this is
  not the case, consider using the continuous registration function
  \code{register.fd}.  Although requiring that a monotone smoother be
  used to estimate the warping functions is safer, it adds considerably
  to the computatation time since monotone smoothing is itself an
  iterative process.  It is usually better to try an initial
  registration with this feature to see if there are any failures of
  monotonicity.  Moreover, monotonicity failures can usually be cured
  by increasing the smoothing parameter defining \code{WfdPar}.  Not
  much curvature is usually required in the warping functions, so a
  rather low power basis, usually B-splines, is suitable for defining
  the functional paramter argument \code{WfdPar}.  A registration with
  a few prominent landmarks is often a good preliminary to using the
  more sophisticated but more lengthy process in \code{register.fd}.
}
\seealso{
  \code{\link{register.fd}},
  \code{\link{smooth.morph}}
}
\examples{
#See the analysis for the lip data in the examples.

## setting parameters
library(lattice)
data(landmark.reg.expData) ## containing an simple object called sampleData

# Preferred:
# eps <- .Machine$double.eps
# to reduce compute time:
eps <- 1000*.Machine$double.eps
from <- -1.0187
to <- 9.4551
# Preferred:
# nb <- 201
# to reduce compute time:
nb <- 31
nbreaks <- 11
## assign landmarks
landmark <- matrix(c(0.4999, 0.657, 0.8141, 0.5523, 0.5523,
                     3.3279, 3.066, 3.0137, 3.2231, 3.2231),
                   ncol=2)
wbasis <- create.bspline.basis(rangeval=c(from, to),
	                       norder=4, breaks=seq(from, to, len=nbreaks))
Wfd0   <- fd(matrix(0,wbasis$nbasis,1),wbasis)
WfdPar <- fdPar(Wfd0, 1, 1e-4)
## get the density of the data
x <- split(sampleData, factor(sampleData$which))
# to save time, reduce the number of curves from 5 to 3
k <- 3
densY <- sapply(x[1:k], function(z){
            r <- range(z[, 1])
            z <- z[, 1]
            z <- z[z>r[1]+eps & z<r[2]-eps]
            density(z, from=from, to=to, n=nb, na.rm=TRUE)$y
        })

argvals <- seq(from, to, len=nb)
fdobj   <- smooth.basis(argvals, densY, wbasis,
                   fdnames = c("x", "samples", "density"))$fd

regDens   <- landmarkreg(fdobj, landmark[1:k,], WfdPar=WfdPar, monwrd=TRUE)

warpfdobj <- regDens$warpfd
warpedX   <- as.matrix(eval.fd(warpfdobj, argvals))
funs      <- apply(warpedX, 2, approxfun, argvals)
## newDat <- list()
## for (i in 1:length(funs))
##   newDat[[names(funs)[i]]] <- data.frame(data=funs[[i]](x[[i]][,1]),
##                                          which=names(funs)[i])
## visualization
## lattice::densityplot(~ data | which, sampleData,
##                    layout=c(1,5), col="#A6CEE3", n=120,
##                     main="Original data set")
## dev.new()
## lattice::densityplot( ~ data | which, layout=c(1,length(newDat)),
##                      do.call(rbind, newDat), col="#A6CEE3", n=120,
##                      main="After warping")


## }
% docclass is function
\keyword{smooth}
