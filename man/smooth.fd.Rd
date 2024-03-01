\name{smooth.fd}
\alias{smooth.fd}
\title{
  Smooth a Functional Data Object Using an Indirectly Specified
  Roughness Penalty
}
\description{
  Smooth data already converted to a functional data object, fdobj,
  using criteria consolidated in a functional data parameter object,
  fdParobj.  For example, data may have been converted to a functional
  data object using function \code{smooth.basis} using a fairly large set of
  basis functions.  This 'fdobj' is then smoothed as specified in
  'fdParobj'.
}
\usage{
smooth.fd(fdobj, fdParobj)
}
\arguments{
  \item{fdobj}{
    a functional data object to be smoothed.
  }
  \item{fdParobj}{
    a functional parameter object. This object is defined by a roughness
    penalty in slot \code{Lfd} and a smoothing parameter lambda in slot
    \code{lambda}, and this information is used to further smooth argument \code{fdobj}.
  }
}
\value{
a functional data object.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{smooth.basis}},
}
\examples{
oldpar <- par(no.readonly=TRUE)
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
yfdPar <- fdPar(yfd, 2, 1e-6)
yfd1 <- smooth.fd(yfd, yfdPar)

#. this code throws an error for. non-cran check
# if (!CRAN()) {
# FIXME: using 3rd derivative here gave error?????
# yfdPar3 <- fdPar(yfd, 3, 1e-6)
# yfd1.3 <- smooth.fd(yfd, yfdPar3)
# Error in bsplinepen(basisobj, Lfdobj, rng) :
#	Penalty matrix cannot be evaluated
#  for derivative of order 3 for B-splines of order 4
# }

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
par(oldpar)
}
\keyword{smooth}
