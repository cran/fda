\name{smooth.basis}
\alias{smooth.basis}
\title{
  Smooth Data with a Roughness Penalty
}
\description{
This is the main function for smoothing data using a roughness penalty.
Unlike function \code{data2fd}, which does not employ a rougness
penalty, this function controls the nature and degree of smoothing by
penalyzing a measure of rougness.  Roughness is definable in a wide
variety of ways using either derivatives or a linear differential
operator.
}
\usage{
smooth.basis(argvals, y, fdParobj, wtvec=rep(1,n),
             dffactor=1,
             fdnames=list(NULL, dimnames(y)[2], NULL))
}
\arguments{
\item{argvals}{
a vector of argument values correspond to the observations in array
\code{y}.
}
\item{y}{
an array containing values of curves at discrete sampling points or
argument values. If the array is a matrix, the rows must correspond to
argumentvalues and columns to replications, and it will be assumed
that there is only one variable per observation.
If \code{y} is a three-dimensional array, the first dimension
corresponds to argument values, the second to replications,
and the third to variables within replications.
If \code{y} is a vector, only one replicate and variable are assumed.
}
\item{fdParobj}{
a functional parameter object, a functional data object or a
functional basis object.  If the object is a functional
parameter object, then the linear differential operator object
and the smoothing parameter in this object define the roughness
penalty.  If the object is a functional data object, the
basis within this object is used without a roughness penalty,
and this is also the case if the object is a functional basis
object.  In these latter two cases, \code{smooth.basis}
is essentially the same as \code{data2fd}.
}
\item{wtvec}{
a vector of the same length as \code{argvals} containing
weights for the values to be smoothed.
}
\item{dffactor}{
Chong Gu in his book \emph{Smoothing Spline ANOVA Models} suggests
a modification of the GCV criterion using a factor modifying
the effective degrees of freedom of the smooth.  He suggests
that values like 1.2 are effective at avoiding undersmoothing
of the data.  The default value of 1 is the classic definition.
}
\item{fdnames}{
a list of length 3 with members containing
\itemize{
	\item a single name for the argument domain, such as 'Time'
	\item a vector of names for the replications or cases
	\item a name for the function, or a vector of names if there are multiple functions.
}
}
}
\value{
a named list of length 7 containing:

\item{fdobj}{
a functional data object that smooths the data.
}
\item{df}{
a degrees of freedom measure of the smooth
}
\item{gcv}{
the value of the generalized cross-validation or GCV criterion.
If there are multiple curves, this is a vector of values, one per
curve.  If the smooth is multivariate, the result is a matrix of
gcv values, with columns corresponding to variables.
}
\item{coef}{
the coefficient matrix or array for the basis function
expansion of the smoothing function
}
\item{SSE}{
the error sums of squares.
SSE is a vector or a matrix of the same size as GCV.
}
\item{penmat}{
the penalty matrix.
}
\item{y2cMap}{
the matrix mapping the data to the coefficients.
}
}
\details{
If the smoothing parameter \code{lambda} is zero,
there is no penalty on roughness.  As
lambda increases, usually in logarithmic terms, the penalty on roughness
increases and the fitted curves become more and more smooth.
Ultimately, the curves are forced to have zero roughness in the sense of
being in the null space of the linear differential operator
object \code{Lfdobj}that is a member of the \code{fdParobj}.

For example, a common choice of roughness penalty is the integrated square
of the second derivative.  This penalizes curvature.  Since the second
derivative of a straight line is zero, very large values of
\code{lambda} will force the fit to become linear.
It is also possible to control the amount of roughness by using a
degrees of freedom measure.  The value equivalent to \code{lambda} is
found in the list returned by the function.  On the other hand, it is
possible to specify a degrees of freedom value, and then use function
\code{df2lambda} to determine the equivalent value of \code{lambda}.
One should not put complete faith in any automatic method for
selecting \code{lambda}, including the GCV method. There are many
reasons for this.  For example, if derivatives are required, then the
smoothing level that is automatically selected may give unacceptably
rough derivatives.  These methods are also highly sensitive to the
assumption of independent errors, which is usually dubious with
functional data.  The best advice is to start with the value minimizing
the \code{gcv} measure, and then explore \code{lambda} values a few
log units up and down from this value to see what the smoothing function
and its derivatives look like.  The function \code{plotfit.fd} was
designed for this purpose.
An alternative to using \code{smooth.basis} is to first represent
the data in a basis system with reasonably high resolution using
\code{data2fd}, and then smooth the resulting functional data object
using function \code{smooth.fd}.
}
\seealso{
\code{\link{data2fd}}, 
\code{\link{df2lambda}}, 
\code{\link{lambda2df}}, 
\code{\link{lambda2gcv}}, 
\code{\link{plot.fd}}, 
\code{\link{project.basis}}, 
\code{\link{smooth.fd}}, 
\code{\link{smooth.monotone}}, 
\code{\link{smooth.pos}}
}
\examples{

#  Shows the effects of three levels of smoothing
#  where the size of the third derivative is penalized.
#  The null space contains quadratic functions.
x <- seq(-1,1,0.02)
y <- x + 3*exp(-6*x^2) + rnorm(rep(1,101))*0.2
#  set up a saturated B-spline basis
basisobj <- create.bspline.basis(c(-1,1),101)
lambda   <- 1
fdParobj <- fdPar(basisobj, 2, lambda)
result1  <- smooth.basis(x, y, fdParobj)
yfd1     <- result1$fd
c(result1$df,result1$gcv) #  display df and gcv measures
lambda   <- 1e-4
fdParobj <- fdPar(basisobj, 2, lambda)
result2  <- smooth.basis(x, y, fdParobj)
yfd2     <- result2$fd
c(result2$df,result2$gcv)
lambda   <- 0
fdParobj <- fdPar(basisobj, 2, lambda)
result3  <- smooth.basis(x, y, fdParobj)
yfd3     <- result3$fd
c(result3$df,result3$gcv)
plot(x,y)           # plot the data
lines(yfd1, lty=2)  #  add heavily penalized smooth
lines(yfd2, lty=1)  #  add reasonably penalized smooth
lines(yfd3, lty=3)  #  add smooth without any penalty
legend(-1,3,c("1","0.0001","0"),lty=c(2,1,3))
plotfit.fd(y, x, yfd2)  # plot data and smooth

}
% docclass is function
\keyword{smooth}
