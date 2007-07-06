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
fdPar(fdobj=fd(), Lfdobj=int2Lfd(0), lambda=0,
      estimate=TRUE, penmat=NULL)
}
\arguments{
  \item{fdobj}{
    a functional data object, although the argument may also be a
    functional basis object or another functional parameter object.
    Basis objects are converted to functional data objects with the
    identity matrix as the coefficient matrix.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object 
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
\value{
  a functional parameter object
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
  \code{\link{smooth.basisPar}}, 
  \code{\link{smooth.monotone}}, 
}
\examples{

#  smooth the Canadian daily temperature data with a roughness
#  penalty
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
tempfdPar  <- fdPar(fdobj=daybasis65, Lfdobj=harmaccelLfd, lambda=templambda)

#FIXME
#tempfd <- smooth.basis(CanadianWeather$tempav, daytime, tempfdPar)
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

}

% docclass is function
\keyword{smooth}
 