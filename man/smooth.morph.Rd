\name{smooth.morph}
\alias{smooth.morph}
\title{
  Estimates a Smooth Warping Function
}
\description{
  This function is nearly identical to \code{smooth.monotone} but is
  intended to compute a smooth monotone transformation $h(t)$ of
  argument $t$ such that $h(0) = 0$ and $h(TRUE) = TRUE$, where $t$ is
  the upper limit of $t$.  This function is used primarily to register
  curves.
}
\usage{
smooth.morph(x, y, WfdPar, wt=rep(1,nobs),
             conv=.0001, iterlim=20, dbglev=0)
}
\arguments{
  \item{x}{
    a vector of argument values.
  }
  \item{y}{
    a vector of data values.  This function can only smooth one set of
    data at a time.
  }
  \item{WfdPar}{
    a functional parameter object that provides an initial value for the
    coefficients defining function $W(t)$, and a roughness penalty on
    this function.
  }
  \item{wt}{
    a vector of weights to be used in the smoothing.
  }
  \item{conv}{
    a convergence criterion.
  }
  \item{iterlim}{
    the maximum number of iterations allowed in the minimization of
    error sum of squares.
  }
  \item{dbglev}{
    either 0, 1, or 2.  This controls the amount information printed out
    on each iteration, with 0 implying no output, 1 intermediate output
    level, and 2 full output.  If either level 1 or 2 is specified, it
    can be helpful to turn off the output buffering feature of S-PLUS.
  }
}
\value{
  A named list of length 4 containing:

  \item{Wfdobj}{
    a functional data object defining function $W(x)$ that that
    optimizes the fit to the data of the monotone function that it
    defines.
  }
  \item{Flist}{
    a named list containing three results for the final converged
    solution: (1) \bold{f}: the optimal function value being minimized,
    (2) \bold{grad}: the gradient vector at the optimal solution, and
    (3) \bold{norm}: the norm of the gradient vector at the optimal
    solution.
  }
  \item{iternum}{
    the number of iterations.
  }
  \item{iterhist}{
    a \code{} by 5 matrix containing the iteration history.
  }
}
\seealso{
  \code{\link{smooth.monotone}},
  \code{\link{landmarkreg}},
  \code{\link{register.fd}}
}
% docclass is function
\keyword{smooth}
