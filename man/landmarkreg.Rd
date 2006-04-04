\name{landmarkreg}
\alias{landmarkreg}
\title{
  Landmark Registration of Functional Observations
}
\description{
It is common to see that among a set of functions certain prominent
features such peaks and valleys, called $landmarks$, do not occur
at the same times, or other
argument values.  This is called $phase variation$, and it can be
essential to align these features before proceeding with further
functional data analyses.  This function uses the timings of these
features to align or register the curves.  The registration involves
estimating a nonlinear transformation of the argument continuum for each
functional observation.  This transformation is called a warping
function. It must be strictly increasing and smooth.
}
\usage{
landmarkreg(fdobj, ximarks, x0marks=xmeanmarks,
            WfdPar=fdPar(defbasis), monwrd=FALSE)
}
\arguments{
\item{fdobj}{
a functional data object containing the curves to be registered.
}
\item{ximarks}{
a matrix containing the timings or argument values associated with the
landmarks for the observations in \code{fd} to be registered.  The
number of rows N
equals the number of observations, and the number of columns NL equals the
number of landmarks. These landmark times must be in the interior of the
interval over which the functions are defined.
}
\item{x0marks}{
a vector of length NL of times of landmarks for target curve.  If
not supplied, the mean of the landmark times in \code{ximarks}
is used.
}
\item{WfdPar}{
a functional parameter object defining the warping functions
that transform time in order to register the curves.
}
\item{monwrd}{
A logical value:
if \code{TRUE}, the warping function is estimated using a monotone
smoothing methhod; otherwise, a regular smoothing method is used,
which is not guaranteed to give strictly monotonic warping
functions.
}
}
\value{
a names list of length 2 with components:

\item{fdreg}{
a functional data object for the registered curves.
}
\item{warpfd}{
a functional data object for the warping functions.
}
}
\details{
It is essential that the location of every landmark be clearly
defined
in each of the curves as well as the template function.  If this is not
the case, consider using the continuous registration function
\code{register.fd}.
Although requiring that a monotone smoother be used to estimate the
warping functions is safer, it adds considerably to the computatation
time since monotone smoothing is itself an iterative process.  It is
usually better to try an initial registration with this feature to see
if there are any failures of monotonicity.  Moreover, monotonicity
failures can usually be cured by increasing the smoothing parameter
defining \code{WfdPar.}
Not much curvature is usually required in the warping functions, so
a rather low power basis, usually B-splines, is suitable for
defining the functional paramter argument \code{WfdPar.}
A registration with a few prominent landmarks is often a good
preliminary to using the more sophisticated but more lengthy process in
\code{register.fd.}
}
\seealso{
register.fd, smooth.morph
}
\examples{
#See the analysis for the lip data in the examples.
}
% docclass is function
\keyword{smooth}
