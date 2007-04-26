\name{register.fd}
\alias{register.fd}
\title{
  Register Functional Data Objects Using a Continuous Criterion
}
\description{
Curves are aligned with a target or template curve using a continuous
criterion.  By aligned is meant that the shape of each curve is matched
as closely as possible to that of the target by means of a smooth
increasing transformation of the argument, or a warping function.
}
\usage{
register.fd(y0fd, yfd, WfdParobj, conv=1e-04,
            iterlim=20, dbglev=1, periodic=FALSE, crit=2)
}
\arguments{
\item{y0fd}{
a functional data object defining the target for registration.
It must be univariate and it must define a single functional
observation.
}
\item{yfd}{
a functional data object defining the functions to be registered to
target \code{yfd0}.
}
\item{WfdParobj}{
a functional parameter object for a single function.  This function
is used as the initial value in the estimation of a function
$W(t)$ that defines the warping function $h(t)$ that
registers a particular curve. The object also contains information
on a roughness penalty and smoothing parameter to control the
roughness of $h(t)$.
}
\item{conv}{
a criterion for convergence of the iterations.
}
\item{iterlim}{
a limit on the number of iterations.
}
\item{dbglev}{
either 0, 1, or 2.  This controls the amount information printed out on
each iteration, with 0 implying no output, 1 intermediate output level,
and 2 full output. If output is desired, it is helpful to turn off
the output buffering option of S-PLUS.
}
\item{periodic}{
a logical variable:
if \code{TRUE}, the functions are considered to be
periodic, in which case a constant can be added to all argument
values after they are warped.
}
\item{crit}{
an integer that is either 1 or 2 that indicates the nature of the
continuous registration criterion that is used.  If 1, the criterion is
least squares, and if 2, the criterion is the minimum eigenvalue of a
cross-product matrix.  In general, criterion 2 is to be preferred.
}
}
\value{
a named list of length 3 containing the following components:

\item{regfd}{
A functional data object containing the registered functions.
}
\item{Wfd}{
A functional data object containing the functions $h W(t)$
that define the warping functions $h(t)$.
}
\item{shift}{
If the functions are periodic, this is a vector of time shifts.
}
}
\details{
The warping function that smoothly and monotonely transforms the
argument is defined by \code{Wfd} is the same as that
defines the monotone
smoothing function in for function \code{smooth.monotone.}  See the
help file for that function for further details.
}
\seealso{
\code{\link{smooth.monotone}}, 
\code{\link{smooth.morph}}
}
\examples{
#See the analyses of the growth data for examples.
}
% docclass is function
\keyword{smooth}
