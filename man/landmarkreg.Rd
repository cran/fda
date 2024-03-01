\name{landmarkreg}
\alias{landmarkreg}
\title{
  Landmark Registration of Functional Observations with Differing Ranges
}
\description{
It is common to see that among a set of functions certain prominent
features such peaks and valleys, called $landmarks$, do not occur
at the same times, or other argument values.  
This is called $phase variation$, and it can be essential to align 
these features before proceeding with further functional data analyses.  

This function uses the timings of these
features to align or register the curves.  The registration involves
estimating a nonlinear transformation of the argument continuum for each
functional observation.  This transformation is called a warping
function. It must be strictly increasing and smooth.

Warning:  As of March 2022, landmark registration cannot be done using
function \code{smooth.basis} instead of function \code{smooth.morph}.  The 
warping function must be strictly monotonic, and we have found that using 
\code{smooth.basis} too often violates this contraint.  Function 
\code{smooth.morph} ensures monotonicity.
}
\usage{
landmarkreg(unregfd, ximarks, x0marks, x0lim, 
             WfdPar=NULL, WfdPar0=NULL, ylambda=1e-10)
}
\arguments{
  \item{unregfd}{
    a functional data object containing the curves to be registered.
  }
  \item{ximarks}{
    A matrix containing the timings or argument values associated with
    the landmarks for the observations in \code{fd} to be registered.
    The number of rows N equals the number of observations, and the
    number of columns NL equals the number of landmarks. These landmark
    times must be in the interior of the interval over which the
    functions are defined.
    Object ximarks may also be a vector or a data.frame.
  }
  \item{x0marks}{
    a vector of length NL of times of landmarks for target curve.  If
    not supplied, the mean of the landmark times in \code{ximarks} is
    used.
    Object x0marks may also be a vector.
  }
  \item{x0lim}{
    A vector of length 2 containing the lower and upper boundary of 
    the interval containing \code{x0marks}.
  }
  \item{WfdPar}{
    a functional parameter object defining the warping functions that
    transform time in order to register the curves.
  }
  \item{WfdPar0}{
    a functional parameter object defining the inverse warping functions 
    that transform time in order to register the curves.
  }
  \item{ylambda}{
    Smoothing parameter controlling the smoothness of the registered
    functions.  It can happen with high dimensional bases that local
    wiggles can appear in the registered curves or their derivatives
    that are not seen in the unregistered versions.  In this case,
    this parameter should be increased to the point where they
    disappear.
  }
}
\value{
  a named list of length 4 with components:

  \item{fdreg}{
    a functional data object for the registered curves.
  }
  \item{warpfd}{
    a functional data object for the warping functions.
  }
  \item{warpinvfd}{
    a functional data object for the inverse warping functions.
  }
  \item{Wfd}{
    a functional data object for the core function defining
    the strictly monotone warping function.
  }
}
\details{
  A value of an arbitrary strictly monotone function at a point \code{x}
  can be defined as the indefinite integral from a fixed lower boundary 
  to \code{x} of the exponential of an unconstrained function value 
  \code{W(x)}.
  
  We use B-spline basis functions to define function \code{W}, and optimize
  the coefficients of its B-spline expansion with respect to the
  data we are fitting.  The optimized core function \code{W} is output 
  along with the registered functions, the warping function qnd the inverse
  warping function.
  
  It is essential that the location of every landmark be clearly defined
  in each of the curves as well as the template function.  If this is
  not the case, consider using the continuous registration function
  \code{register.fd}.  
  
  Not much curvature is usually required in the warping functions, so a
  rather lower power basis, usually B-splines, is suitable for defining
  the functional parameter argument \code{WfdPar}.  A registration with
  a few prominent landmarks is often a good preliminary to using the
  more sophisticated but more lengthy process in \code{register.fd}.
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
  \code{\link{register.fd}},
  \code{\link{landmarkreg}},
  \code{\link{smooth.morph}}
}
\keyword{smooth}
