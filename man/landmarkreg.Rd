\name{landmarkreg}
\alias{landmarkreg}
\title{Registration of Functional Observations Using Landmarks }
\description{
It is common to see that among a set of functions certain prominent features such peaks and valleys, called landmarks, do not occur at the same times, or other argument values. This is called phase variation, and it can be essential to align these features before proceeding with further functional data analyses. This function uses the timings of these features to align or register the curves. The registration involves estimating a nonlinear transformation of the argument continuum for each functional observation. This transformation is called a warping function. It must be strictly increasing and smooth. 
}
\usage{
landmarkreg(fd, fd0, ximarks, x0marks=xmeanmarks, wbasis=basis, Lfd=0, sparval=1e-010, monwrd=FALSE)
}
\arguments{
  \item{fd }{A functional data object containing the curves to be registered. }
  \item{fd0 }{A functional data object with a single observation defining a template or gold standard function to which the other functions are registered. }
  \item{ximarks }{A matrix containing the timings or argument values associated with the landmarks for the observations in fd to be registered. The number of rows equals the number of observations, and the number of columns equals the number of landmarks. }
  \item{x0marks }{A vector containing the timings or argument values associated with the landmarks for the template or gold standard function in fd0. }
  \item{wbasis }{A basis object used to define the warping functions. }
  \item{Lfd }{Either a nonnegative integer of a linear differential operator to be applied to the functions before registering them. }
  \item{sparval }{A smoothing parameter controlling the smoothness of the estimated warping functions. }
  \item{monwrd }{If TRUE, the warping function is estimated using a monotone smoothing methhod; otherwise, a regular smoothing method is used, which is not guaranteed to give strictly monotonic warping functions. }
}
\details{
It is essential that the location of every landmark be clearly defined in each of the curves as well as the template function. If this is not the case, consider using the continuous registration function register.fd. 

Although requiring that a monotone smoother be used to estimate the warping functions is safer, it adds considerably to the computatation time since monotone smoothing is itself an iterative process. It is usually better to try an initial registration with this feature to see if there are any failures of monotonicity. Moreover, monotonicity failures can usually be cured by increasing the smoothing parameter sparval. 

Not much curvature is usually required in the warping functions, so a rather low power basis, usually B-splines, is suitable for argument wbasis. 

A registration with a few prominent landmarks is often a good preliminary to using the more sophisticated but more lengthy process in register.fd. 

}
\value{
A list with these entries is returned. 
"fdreg" a functional data object for the registered curves 
"warpfd" A functional data object for the warping functions 
}
\references{
Ramsay, J. O. and Silverman, B. W. (1997) Functional Data Analysis. New York: Springer. Registration is described in Chapter 5. 
}
\note{}
\author{}
\seealso{
register.fd 
}
\examples{
}
\keyword{}
