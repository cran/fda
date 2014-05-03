\name{intensity.fd}
\alias{intensity.fd}
\title{
  Intensity Function for Point Process
}
\description{
The intensity $mu$ of a series of event times that obey a
homogeneous Poisson process is the mean number of events per unit time.
When this event rate varies over time, the process is said to be
nonhomogeneous, and $mu(t)$, and is estimated by this function
\code{intensity.fd}.
}
\usage{
intensity.fd(x, WfdParobj, conv=0.0001, iterlim=20,
             dbglev=1, returnMatrix=FALSE)
}
\arguments{
\item{x}{
a vector containing a strictly increasing series of event times.
These event times assume that the the events begin to be observed
at time 0, and therefore are times since the beginning of
observation.
}
\item{WfdParobj}{
a functional parameter object estimating the log-intensity function
$W(t) = log[mu(t)]$ .
Because the intensity function $mu(t)$ is necessarily positive,
it is represented by \code{mu(x) = exp[W(x)]}.
}
\item{conv}{
a convergence criterion, required because the estimation
process is iterative.
}
\item{iterlim}{
maximum number of iterations that are allowed.
}
\item{dbglev}{
either 0, 1, or 2.  This controls the amount information printed out on
each iteration, with 0 implying no output, 1 intermediate output level,
and 2 full output.  If levels 1 and 2 are used, turn off the output
buffering option.
}
\item{returnMatrix}{
logical:  If TRUE,  a two-dimensional is returned using a
special class from the Matrix package.
}
}
\value{
a named list of length 4 containing:

\item{Wfdobj}{
a functional data object defining function $W(x)$ that that
optimizes the fit to the data of the monotone function that it defines.
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
\item{iterhist}{
a \code{iternum+1} by 5 matrix containing the iteration
history.
}
}
\details{


The intensity function $I(t)$ is almost the same thing as a
probability density function $p(t)$ estimated by function
\code{densify.fd}.  The only difference is the absence of
the normalizing constant $C$ that a density function requires
in order to have a unit integral.
The goal of the function is provide a smooth intensity function
estimate that approaches some target intensity by an amount that is
controlled by the linear differential operator \code{Lfdobj} and
the penalty parameter in argument \code{WfdPar}.
For example, if the first derivative of
$W(t)$ is penalized heavily, this will force the function to
approach a constant, which in turn will force the estimated Poisson
process itself to be nearly homogeneous.
To plot the intensity function or to evaluate it,
evaluate \code{Wfdobj}, exponentiate the resulting vector.
}
\seealso{
\code{\link{density.fd}}
}
\examples{

#  Generate 101 Poisson-distributed event times with
#  intensity or rate two events per unit time
N  <- 101
mu <- 2
#  generate 101 uniform deviates
uvec <- runif(rep(0,N))
#  convert to 101 exponential waiting times
wvec <- -log(1-uvec)/mu
#  accumulate to get event times
tvec <- cumsum(wvec)
tmax <- max(tvec)
#  set up an order 4 B-spline basis over [0,tmax] with
#  21 equally spaced knots
tbasis <- create.bspline.basis(c(0,tmax), 23)
#  set up a functional parameter object for W(t),
#  the log intensity function.  The first derivative
#  is penalized in order to smooth toward a constant
lambda <- 10
Wfd0 <- fd(matrix(0,23,1),tbasis)
WfdParobj <- fdPar(Wfd0, 1, lambda)
#  estimate the intensity function
Wfdobj <- intensity.fd(tvec, WfdParobj)$Wfdobj
#  get intensity function values at 0 and event times
events <- c(0,tvec)
intenvec <- exp(eval.fd(events,Wfdobj))
#  plot intensity function
plot(events, intenvec, type="b")
lines(c(0,tmax),c(mu,mu),lty=4)

}
% docclass is function
\keyword{smooth}
