\name{register.fd}
\alias{register.fd}
\title{
  Register Functional Data Objects Using a Continuous Criterion
}
\description{
  A function is said to be aligned or registered with a target function
  if its  salient features, such as peaks, valleys and crossings of
  fixed thresholds,  occur at about the same argument values as those of
  the target.  Function \code{register.fd} aligns these features by
  transforming or warping  the argument domain of each function in a
  nonlinear but strictly  order-preserving fashion.  Multivariate
  functions may also be registered.  If the domain is time, we say that
  this transformation transforms clock time to system time.  The
  transformation itself is called a warping function.
}
\usage{
register.fd(y0fd=NULL, yfd=NULL, WfdParobj=NULL,
            conv=1e-04, iterlim=20, dbglev=1, periodic=FALSE, crit=2)
}
\arguments{
  \item{y0fd}{
    a functional data object defining one or more target functions for
    registering the functions in argument \code{yfd}. If the functions
    to be registered are univariate, then \code{y0fd} may contain only a
    single function, or it may contain as many functions as are in
    \code{yfd}.  If \code{yfd} contains multivariate functions, then
    \code{y0fd} may either as many functions as there are variables in
    \code{yfd}, or it may contain as many functions as are in \code{yfd}
    and these functions must then be multivariate and be of the same
    dimension as those in \code{yfd}.

    If \code{yfd} is supplied as a named argument and \code{y0fd} is
    not, then \code{y0fd} is computed inside the function to be the mean
    of the functions in \code{yfd}.

    If the function is called with a single unnamed argument, and there
    is no other argument that is named as \code{y0fd} then this unnamed
    argument is taken to be actually \code{yfd} rather than \code{y0fd},
    and then also \code{y0fd} is computed to be the man of the functions
    supplied.
  }
  \item{yfd}{
    a functional data object defining the functions to be
    registered to target \code{y0fd}.  The functions may be either
    univariate or multivariate.

    If \code{yfd} contains a single multivariate function is to be
    registered, it essential that the coefficient array for \code{y0fd}
    have class \code{array}, have three dimensions, and that its second
    dimension be of length 1.
  }
  \item{WfdParobj}{
    a functional parameter object containing either a single function or
    the same number of functions as are contained in \code{yfd}.  The
    coefficients supply the initial values in the estimation of a
    functions $W(t)$ that defines the warping functions $h(t)$ that
    register the set of curves.  \code{WfdParobj} also defines the
    roughness penalty and smoothing parameter used to control the
    roughness of $h(t)$.

    The basis used for this object must be a B-spline type, and the
    order of the B-spline basis must be at least 2 (piecewise linear).

    If WFDPAROBJ is not supplied, it is constructed from a bspline basis
    of order 2 with 2 basis functions; that is, a basis for piecewise
    linear functions.  The smoothing parameter lambda for this default
    is 0.
  }
  \item{conv}{
    a criterion for convergence of the iterations.
  }
  \item{iterlim}{
    a limit on the number of iterations.
  }
  \item{dbglev}{
    either 0, 1, or 2.  This controls the amount information printed out
    on each iteration, with 0 implying no output, 1 intermediate output
    level, and 2 full output.  R normally postpones displaying these
    results until the entire computation is computed, an option that it
    calls "output buffering."  Since the total computation time may be
    considerable, one may opt for turning this feature off by
    un-checking this box in the "Misc" menu item in the R Console.
  }
  \item{periodic}{
    a logical variable:  if \code{TRUE}, the functions are considered to
    be periodic, in which case a constant can be added to all argument
    values after they are warped.
  }
  \item{crit}{
    an integer that is either 1 or 2 that indicates the nature of the
    continuous registration criterion that is used.  If 1, the criterion
    is least squares, and if 2, the criterion is the minimum eigenvalue
    of a cross-product matrix.  In general, criterion 2 is to be
    preferred.
  }
}
\value{
  a named list of length 4 containing the following components:

  \item{regfd}{
    A functional data object containing the registered functions.
  }
  \item{warpfd}{
    A functional data object containing the warping functions $h(t)$.
  }
  \item{Wfd}{
    A functional data object containing the functions $h W(t)$
    that define the warping functions $h(t)$.
  }
  \item{shift}{
    If the functions are periodic, this is a vector of time shifts.
  }
  \item{y0fd}{
    The target function object y0fd.
  }
  \item{yfd}{
    The function object yfd containing the functions to be registered.
  }
}
\details{
  The warping function that smoothly and monotonely transforms the
  argument is defined by \code{Wfd} is the same as that defines the
  monotone smoothing function in for function \code{smooth.monotone.}
  See the help file for that function for further details.
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York, ch. 6 & 7.
}
\seealso{
  \code{\link{smooth.monotone}},
  \code{\link{smooth.morph}}
  \code{\link{plotreg.fd}}
  \code{\link{register.newfd}}
}
\examples{
#See the analyses of the growth data for examples.

##
## 1.  smooth the growth data for the Berkeley boys
##
# Specify smoothing weight
lambda.gr2.3 <- .03
# Specify what to smooth, namely the rate of change of curvature
Lfdobj.growth    <- 2
# Set up a B-spline basis for smoothing the discrete data
nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
rng.growth <- range(growth$age)
wbasis.growth <- create.bspline.basis(rangeval=rng.growth,
                   nbasis=nbasis.growth, norder=norder.growth,
                   breaks=growth$age)
# Smooth the data
# in afda-ch06.R, and register to individual smooths:
cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth  <- fd(cvec0.growth, wbasis.growth)
growfdPar2.3 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr2.3)
hgtmfd.all   <- with(growth, smooth.basis(age, hgtm, growfdPar2.3)$fd)
# Register the growth velocity rather than the
# growth curves directly
smBv <- deriv.fd(hgtmfd.all, 1)

##
## 2.  Register the first 2 Berkeley boys using the default basis
##     for the warping function
##
# register.fd takes time, so we use only 2 curves as an illustration
# to minimize computing time in these examples
nBoys <- 2
#  Define the target function as the mean of the first nBoys records
smBv0 = mean.fd(smBv[1:nBoys])
#  Register these curves.  The default choice for the functional
#  parameter object WfdParObj is used.
smB.reg.0 <- register.fd(smBv0, smBv[1:nBoys])
#  plot each curve.  Click on the R Graphics window to show each plot.
#  The left panel contains:
#    -- the unregistered curve (dashed blue line)
#    -- the target function (dashed red line)
#    -- the registered curve (solid blue line)
#  The right panel contains:
#    -- the warping function h(t)
#    -- the linear function corresponding to no warping
plotreg.fd(smB.reg.0)
#  Notice that all the warping functions all have simple shapes
#  due to the use of the simplest possible basis

if (!CRAN()) {
##
## 3.  Define a more flexible basis for the warping functions
##
Wnbasis   <- 4
Wbasis    <- create.bspline.basis(rng.growth, Wnbasis)
Wfd0      <- fd(matrix(0,Wnbasis,1),Wbasis)
#  set up the functional parameter object using only
#      a light amount smoothing
WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=0.01)
#  register the curves
smB.reg.1 <- register.fd(smBv0, smBv[1:nBoys], WfdParobj)
plotreg.fd(smB.reg.1)
#  Notice that now the warping functions can have more complex shapes

##
## 4.  Change the target to the mean of the registered functions ...
##     this should provide a better target for registration
##
smBv1 <- mean.fd(smB.reg.1$regfd)
#  plot the old and the new targets
par(mfrow=c(1,1),ask=FALSE)
plot(smBv1)
lines(smBv0, lty=2)
#  Notice how the new target (solid line) has sharper features and
#  a stronger pubertal growth spurt relative to the old target
#  (dashed line).  Now register to the new target
smB.reg.2 <- register.fd(smBv1, smBv[1:nBoys], WfdParobj)
plotreg.fd(smB.reg.2)
#  Plot the mean of these curves as well as the first and second targets
par(mfrow=c(1,1),ask=FALSE)
plot(mean.fd(smB.reg.2$regfd))
lines(smBv0, lty=2)
lines(smBv1, lty=3)
#  Notice that there is almost no improvement over the age of the
#  pubertal growth spurt, but some further detail added in the
#  pre-pubertal region.  Now register the previously registered
#  functions to the new target.
smB.reg.3 <- register.fd(smBv1, smB.reg.1$regfd, WfdParobj)
plotreg.fd(smB.reg.3)
#  Notice that the warping functions only deviate from the straight line
#  over the pre-pubertal region, and that there are some small adjustments
#  to the registered curves as well over the pre-pubertal region.
}

##
## 5.  register and plot the angular acceleration of the gait data
##
gaittime  <- as.matrix(0:19)+0.5
gaitrange <- c(0,20)
#  set up a fourier basis object
gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)
# set up a functional parameter object penalizing harmonic acceleration
harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
gaitfdPar    <- fdPar(gaitbasis, harmaccelLfd, 1e-2)
#  smooth the data
gaitfd <- smooth.basis(gaittime, gait, gaitfdPar)$fd
#  compute the angular acceleration functional data object
D2gaitfd    <- deriv.fd(gaitfd,2)
names(D2gaitfd$fdnames)[[3]] <- "Angular acceleration"
D2gaitfd$fdnames[[3]] <- c("Hip", "Knee")
#  compute the mean angular acceleration functional data object
D2gaitmeanfd  <- mean.fd(D2gaitfd)
names(D2gaitmeanfd$fdnames)[[3]] <- "Mean angular acceleration"
D2gaitmeanfd$fdnames[[3]] <- c("Hip", "Knee")
#  register the functions for the first 2 boys
#  argument periodic = TRUE causes register.fd to estimate a horizontal shift
#  for each curve, which is a possibility when the data are periodic
#  set up the basis for the warping functions
nwbasis   <- 4
wbasis    <- create.bspline.basis(gaitrange,nwbasis,3)
Warpfd    <- fd(matrix(0,nwbasis,nBoys),wbasis)
WarpfdPar <- fdPar(Warpfd)
#  register the functions
gaitreglist <- register.fd(D2gaitmeanfd, D2gaitfd[1:nBoys], WarpfdPar,
                           periodic=TRUE)
#  plot the results
plotreg.fd(gaitreglist)
#  display horizonal shift values
print(round(gaitreglist$shift,1))
}
\keyword{smooth}
