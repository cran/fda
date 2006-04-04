\name{smooth.monotone}
\alias{smooth.monotone}
\title{
  Monotone Smoothing of Data
}
\description{
When the discrete data that are observed reflect a smooth strictly
increasing or strictly decreasing function, it is often desirable to
smooth the data with a strictly monotone function, even though the data
themselves may not be monotone due to observational error.  An example
is when data are collected on the size of a growing organism over time.
This function computes such a smoothing function, but, unlike other
smoothing functions, for only for one curve at a time.
The smoothing function minimizes a weighted error sum of squares
criterion.  This minimization requires iteration, and therefore is more
computationally intensive than normal smoothing.
}
\usage{
smooth.monotone(x, y, wt=rep(1,nobs), WfdParobj,
                zmat=matrix(1,nobs,1),
                conv=.0001, iterlim=20,
                active=c(FALSE,rep(TRUE,ncvec-1)),
                dbglev=1)
}
\arguments{
\item{x}{
a vector of argument values.
}
\item{y}{
a vector of data values.  This function can only smooth
one set of data at a time.
}
\item{WfdParobj}{
a functional parameter object that provides an initial
value for the coefficients defining function $W(t)$,
and a roughness penalty on this function.
}
\item{wt}{
a vector of weights to be used in the smoothing.
}
\item{zmat}{
a design matrix or a matrix of covariate values that also
define the smooth of the data.
}
\item{conv}{
a convergence criterion.
}
\item{iterlim}{
the maximum number of iterations allowed in the minimization
of error sum of squares.
}
\item{active}{
a logical vector specifying which coefficients defining
$W(t)$ are estimated.  Normally, the first coefficient
is fixed.
}
\item{dbglev}{
either 0, 1, or 2.  This controls the amount information printed out on
each iteration, with 0 implying no output, 1 intermediate output level,
and 2 full output.  If either level 1 or 2 is specified, it can be
helpful to turn off the output buffering feature of S-PLUS.
}
}
\value{
a named list of length 5 containing:

\item{Wfdobj}{
a functional data object defining function $W(x)$ that that
optimizes the fit to the data of the monotone function that it defines.
}
\item{beta}{
the optimal regression coefficient values.
}
\item{Flist}{
a named list containing three results for the final converged solution:
(1)
\bold{f}: the optimal function value being minimized,
(2)
\bold{grad}: the gradient vector at the optimal solution,   and
(3)
\bold{norm}: the norm of the gradient vector at the optimal solution.
}
\item{iternum}{
the number of iterations.
}
\item{iternum}{
the number of iterations.
}
\item{iterhist}{
\code{iternum+1} by 5 matrix containing the iteration
history.
}
}
\details{
The
smoothing function  $f(x)$ is determined by
three objects that need to be estimated from the data:

\itemize{
	\item $W(x)$, a functional data object that is first exponentiated and then the result integrated.  This is the heart of the monotone smooth.  The closer $W(x)$ is to zero, the closer the monotone smooth becomes a straight line.  The closer $W(x)$ becomes a constant, the more the monotone smoother becomes an exponential function.  It is assumed that $W(0) = 0.$

	\item $b0$, an intercept term that determines the value of the smoothing function at $x = 0$. 

	\item $b1$, a regression coefficient that determines the slope of the smoothing function at $x = 0$.
}
In addition, it is possible to have the intercept $b0$
depend in turn on the values of one or more covariates through the
design matrix \code{Zmat} as follows:
$b0 = Z c$. In this case, the single
intercept coefficient is replaced by the regression coefficients
in vector $c$ multipying the design matrix.

}
\seealso{
\code{\link{smooth.basis}}, 
\code{\link{smooth.pos}}, 
\code{\link{smooth.morph}}
}
\examples{

#  Estimate the acceleration functions for growth curves
#  See the analyses of the growth data.
#  Set up the ages of height measurements for Berkeley data
age <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))
#  Range of observations
rng <- c(1,18)
#  First set up a basis for monotone smooth
#  We use b-spline basis functions of order 6
#  Knots are positioned at the ages of observation.
norder <- 6
nage   <- 31
nbasis <- nage + norder - 2
wbasis <- create.bspline.basis(rng, nbasis, norder, age)
#  starting values for coefficient
cvec0 <- matrix(0,nbasis,1)
Wfd0  <- fd(cvec0, wbasis)
#  set up functional parameter object
Lfdobj    <- 3          #  penalize curvature of acceleration
lambda    <- 10^(-0.5)  #  smoothing parameter
growfdPar <- fdPar(Wfd0, Lfdobj, lambda)
#  Set up wgt vector
wgt   <- rep(1,nage)
#  Smooth the data for the first girl
hgt1 = growth$hgtf[,1]
result <- smooth.monotone(age, hgt1, wgt, growfdPar)
#  Extract the functional data object and regression
#  coefficients
Wfd  <- result$Wfdobj
beta <- result$beta
#  Evaluate the fitted height curve over a fine mesh
agefine <- seq(1,18,len=101)
hgtfine <- beta[1] + beta[2]*eval.monfd(agefine, Wfd)
#  Plot the data and the curve
plot(age, hgt1, type="p")
lines(agefine, hgtfine)
#  Evaluate the acceleration curve
accfine <- beta[2]*eval.monfd(agefine, Wfd, 2)
#  Plot the acceleration curve
plot(agefine, accfine, type="l")
lines(c(1,18),c(0,0),lty=4)

}
% docclass is function
\keyword{smooth}
