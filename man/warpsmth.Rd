\name{warpsmth}
\alias{warpsmth}
\title{A Smooth Monotone Warping Function }
\description{
This function does nearly the same thing that function smooth.monotone does, but in this case the constants b0 and b1 are determined so that the initial value of the warping function is the same as the initial argument value, and the final warping function value is the same as the final argument values. This function is mainly used in function landmarkreg to ensure that the warping function is monotone. 
}
\usage{
warpsmth(x, y, wt=rep(1, nobs), Wfd, Lfd=1, lambda=0, conv=0.0001, iterlim=20, active=rep(TRUE, nbasis), dbglev=0)
}
\arguments{
  \item{x }{A vector of argument values associated with the data. }
  \item{y }{A vector of observations to be smoothed. }
  \item{wt }{A vector of positive weights for the observations. The default is all weights one. }
  \item{Wfd }{A functional data object that determines the final monotone smooth. It must be univariate and define a single functional observation. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object that defines a rougness penalty to be applied to functional data object Wfdobj. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
  \item{lambda }{A nonnegative smoothing parameter that controls the degree of roughness in Wfdobj, and consequently in the monotone smoothing function. }
  \item{conv }{A criterion for convergence of the iterations. }
  \item{iterlim }{A limit on the number of iterations. }
  \item{active }{A logical vector of length equal to the number of coefficients defining Wfdobj. If an entry is T, the corresponding coefficient is estimated, and if F, it is held at the value defining the argument Wfdobj. Normally the first coefficient is set to 0 and not estimated, since it is assumed that W(0) = 0. }
  \item{dbglev }{Either 0, 1, or 2. This controls the amount information printed out on each iteration, with 0 implying no output, 1 intermediate output level, and 2 full output. }
}
\details{
See the help file for smooth.monotone. 
}
\value{
See the help file for smooth.monotone.
}
\references{}
\note{}
\author{}
\seealso{}
\examples{
}
\keyword{}
