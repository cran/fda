\name{smooth.fd}
\alias{smooth.fd}
\title{
  Smooth a Functional Data Object Using a Roughness Penalty
}
\description{
This function is intended to apply a roughness penalized smooth to
data already set up as a functional data object.  For example, data
may have been converted to a functional data object using function
\code{data2fd} using a fairly large set of basis functions, and
subsequently it was desired to smooth the functional data object
that resulted.
}
\usage{
smooth.fd(fdobj, fdParobj)
}
\arguments{
\item{fdobj}{
a functional data object to be smoothed.
}
\item{fdParobj}{
a functional parameter object. This object is defined by a
roughness penalty in slot \code{Lfd} and a smoothing
parameter lambda in slot \code{lambda}, and this
information is used to further smooth argument \code{fdobj}.
}
}
\value{
a functional data object.
}
\seealso{
\code{\link{smooth.basis}}, 
\code{\link{data2fd}}
}
\examples{

#  Shows the effects of two levels of smoothing
#  where the size of the third derivative is penalized.
#  The null space contains quadratic functions.
x <- seq(-1,1,0.02)
y <- x + 3*exp(-6*x^2) + rnorm(rep(1,101))*0.2
#  set up a saturated B-spline basis
basisobj <- create.bspline.basis(c(-1,1),81)
#  convert to a functional data object that interpolates the data.
result <- smooth.basis(x, y, basisobj)
yfd  <- result$fd
#  set up a functional parameter object with smoothing
#  parameter 1e-6 and a penalty on the 3rd derivative.
# FIXME: using 3rd derivative here gave error
yfdPar <- fdPar(yfd, 2, 1e-6)
yfd1 <- smooth.fd(yfd, yfdPar)
#  set up a functional parameter object with smoothing
#  parameter 1 and a penalty on the 3rd derivative.
yfdPar <- fdPar(yfd, 2, 1)
yfd2 <- smooth.fd(yfd, yfdPar)
#  plot the data and smooth
plot(x,y)           # plot the data
lines(yfd1, lty=1)  #  add moderately penalized smooth
lines(yfd2, lty=3)  #  add heavily  penalized smooth
legend(-1,3,c("0.000001","1"),lty=c(1,3))
#  plot the data and smoothing using function plotfit.fd
plotfit.fd(y, x, yfd1)  # plot data and smooth

}
% docclass is function
\keyword{smooth}
