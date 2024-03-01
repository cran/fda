\name{fdPar}
\alias{fdPar}
\title{
  Define a Functional Parameter Object
}
\description{
  Functional parameter objects are used as arguments to functions that
  estimate functional parameters, such as smoothing functions like
  \code{smooth.basis}.  A functional parameter object is a functional
  data object with additional slots specifying a roughness penalty, a
  smoothing parameter and whether or not the functional parameter is to
  be estimated or held fixed.  Functional parameter objects are used as
  arguments to functions that estimate functional parameters.
}
\usage{
  fdPar(fdobj=NULL, Lfdobj=NULL, lambda=0, estimate=TRUE, penmat=NULL)
}
\arguments{
  \item{fdobj}{
    a functional data object, functional basis object, a functional
    parameter object or a matrix.  If it a matrix, it is replaced by
    fd(fdobj).  If class(fdobj) == 'basisfd', it is converted to an
    object of class \code{fd} with a coefficient matrix consisting of a
    single column of zeros.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object.

    If \code{NULL}, Lfdobj depends on fdobj[['basis']][['type']]:

    \describe{
      \item{bspline}{Lfdobj <- int2Lfd(max(0, norder-2)), 
      where norder = norder(fdobj).}
      \item{fourier}{Lfdobj = a harmonic acceleration operator:
	    \code{Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng)}
      where rng = fdobj[['basis']][['rangeval']].}
      \item{anything else}{Lfdobj <- int2Lfd(0)}
    }
  }
  \item{lambda}{
    a nonnegative real number specifying the amount of smoothing
    to be applied to the estimated functional parameter.
  }
  \item{estimate}{not currently used.}
  \item{penmat}{
    a roughness penalty matrix.  Including this can eliminate the need
    to compute this matrix over and over again in some types of
    calculations.
  }
}
\details{
  Functional parameters are often needed to specify initial
  values for iteratively refined estimates, as is the case in
  functions \code{register.fd} and \code{smooth.monotone}.

  Often a list of functional parameters must be supplied to a function
  as an argument, and it may be that some of these parameters are
  considered known and must remain fixed during the analysis.  This is
  the case for functions \code{fRegress} and  \code{pda.fd}, for
  example.
}
\value{
  a functional parameter object (i.e., an object of class \code{fdPar}),
  which is a list with the following components:

  \item{fd}{
    a functional data object (i.e., with class \code{fd})
  }
  \item{Lfd}{
    a linear differential operator object (i.e., with class
    \code{Lfd})
  }
  \item{lambda}{
    a nonnegative real number
  }
  \item{estimate}{not currently used}
  \item{penmat}{
    either NULL or a square, symmetric matrix with penmat[i, j] =
    integral over fd[['basis']][['rangeval']] of basis[i]*basis[j]
  }
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York
}
\seealso{
  \code{\link{cca.fd}},
  \code{\link{density.fd}},
  \code{\link{fRegress}},
  \code{\link{intensity.fd}},
  \code{\link{pca.fd}},
  \code{\link{smooth.fdPar}},
  \code{\link{smooth.basis}},
  \code{\link{smooth.monotone}},
  \code{\link{int2Lfd}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
##
## Simple example
##
#  set up range for density
rangeval <- c(-3,3)
#  set up some standard normal data
x <- rnorm(50)
#  make sure values within the range
x[x < -3] <- -2.99
x[x >  3] <-  2.99
#  set up basis for W(x)
basisobj <- create.bspline.basis(rangeval, 11)
#  set up initial value for Wfdobj
Wfd0 <- fd(matrix(0,11,1), basisobj)
WfdParobj <- fdPar(Wfd0)

WfdP3 <- fdPar(seq(-3, 3, length=11))

##
##  smooth the Canadian daily temperature data
##
#    set up the fourier basis
nbasis   <- 365
dayrange <- c(0,365)
daybasis <- create.fourier.basis(dayrange, nbasis)
dayperiod <- 365
harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), dayrange)
#  Make temperature fd object
#  Temperature data are in 12 by 365 matrix tempav
#    See analyses of weather data.
#  Set up sampling points at mid days
daytime  <- (1:365)-0.5
#  Convert the data to a functional data object
daybasis65 <- create.fourier.basis(dayrange, nbasis, dayperiod)
templambda <- 1e1
tempfdPar  <- fdPar(fdobj=daybasis65, Lfdobj=harmaccelLfd,
                    lambda=templambda)

#FIXME
#tempfd <- smooth.basis(CanadianWeather$tempav, daytime, tempfdPar)$fd
#  Set up the harmonic acceleration operator
Lbasis  <- create.constant.basis(dayrange);
Lcoef   <- matrix(c(0,(2*pi/365)^2,0),1,3)
bfdobj  <- fd(Lcoef,Lbasis)
bwtlist <- fd2list(bfdobj)
harmaccelLfd <- Lfd(3, bwtlist)
#  Define the functional parameter object for
#  smoothing the temperature data
lambda   <- 0.01  #  minimum GCV estimate
#tempPar <- fdPar(daybasis365, harmaccelLfd, lambda)
#  smooth the data
#tempfd <- smooth.basis(daytime, CanadialWeather$tempav, tempPar)$fd
#  plot the temperature curves
#plot(tempfd)

##
## with rangeval of class Date and POSIXct
##
par(oldpar)

}

% docclass is function
\keyword{smooth}
