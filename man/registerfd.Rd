\name{registerfd}
\alias{registerfd}
\title{Register Functional Data Objects }
\description{
Curves are aligned with a target or template curve using a continuous criterion. By aligned is meant that the shape of each curve is matched as closely as possible to that of the target by means of a smooth increasing transformation of the argument, or a warping function. 
}
\usage{
registerfd(y0fd, yfd, Wfd0, Lfdobj=2, lambda=1, conv=0.01, iterlim=10,
dbglev=1, periodic=FALSE, crit=2)
}
\arguments{
  \item{y0fd }{A functional data object defining the target. It must be univariate and it must define a single functional observation. }
  \item{yfd }{A functional data object defining the functions to be registered to yfd0. }
  \item{Wfd0 }{A functional data object defining the initial value of a function W(t) that determines the warping function. }
  \item{Lfdobj }{Either a nonnegative integer or a linear differential operator object that defines a rougness penalty to be applied to functional data object Wfd0. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
  \item{lambda }{A nonnegative smoothing parameter that controls the degree of roughness in Wfdobj, and consequently in the monotone smoothing function. }
  \item{conv }{A criterion for convergence of the iterations. }
  \item{iterlim }{A limit on the number of iterations. }
  \item{dbglev }{Either 0, 1, or 2. This controls the amount information printed out on each iteration, with 0 implying no output, 1 intermediate output level, and 2 full output. }
  \item{periodic }{A logical variable. If T, the functions are considered to be periodic, in which case a constant can be added to all argument values after they are warped. }
  \item{crit }{An integer that is either 1 or 2 that indicates the nature of the continuous registration criterion that is used. If 1, the criterion is least squares, and if 2, the criterion is the minimum eigenvalue of a cross-product matrix. In general, criterion 2 is to be preferred. }
}
\details{
The warping function that smoothly and monotonely transforms the argument is defined by Wfd0 in the same that the monotone smoothing function is for function smooth.monotone. See the help file for that function for further details. 
}
\value{
A list containing: regstr <- list("regfd"=regfd, "Wfd"=Wfd, "shift"=shift) 
regfd 
The registered versions of the original function defined in yfd. 
Wfd 
The functional data object determining the warping functions. 
shift 
A constant that is added to the warped arguments if the data are considered periodic. Otherwise, zero. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{
}
\keyword{}
