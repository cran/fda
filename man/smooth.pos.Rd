\name{smooth.pos}
\alias{smooth.pos}
\title{
  Smooth Data with a Positive Function
}
\description{
A set of data is smoothed with a functional data object that only
takes positive values.  For example, this function can be used
to estimate a smooth variance function from a set of squared residuals.
A function $W(t)$ is estimated such that that the smoothing
function is $exp[W(t)]$.
}
\usage{
smooth.pos(argvals, y, WfdParobj, wt=rep(1,nobs),
           conv=.0001, iterlim=20, dbglev=1)
}
\arguments{
\item{argvals}{
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
\item{conv}{
a convergence criterion.
}
\item{iterlim}{
the maximum number of iterations allowed in the minimization
of error sum of squares.
}
\item{dbglev}{
either 0, 1, or 2.  This controls the amount information printed out on
each iteration, with 0 implying no output, 1 intermediate output level,
and 2 full output.  If either level 1 or 2 is specified, it can be
helpful to turn off the output buffering feature of S-PLUS.
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
\item{iternum}{
the number of iterations.
}
\item{iterhist}{
a \code{iternum+1} by 5 matrix containing the iteration
history.
}
}
\seealso{
\code{\link{smooth.monotone}}, 
\code{\link{smooth.morph}}
}
\examples{
#See the analyses of the daily weather data for examples.
}
% docclass is function
\keyword{smooth}
