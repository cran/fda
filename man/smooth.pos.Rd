\name{smooth.pos}
\alias{smooth.pos}
\title{
  Smooth Data with a Positive Function
}
\description{
  A set of data is smoothed with a functional data object that only
  takes positive values.  For example, this function can be used to
  estimate a smooth variance function from a set of squared residuals.
  A function $W(t)$ is estimated such that that the smoothing function
  is $exp[W(t)]$.
}
\usage{
smooth.pos(argvals, y, WfdParobj, wtvec=rep(1,n),
           conv=.0001, iterlim=50, dbglev=1)
}
\arguments{
  \item{argvals}{
    Argument value array of length N, where N is the number of observed
    curve values for each curve.  It is assumed that that these argument
    values are common to all observed curves.  If this is not the case,
    you will need to run this function inside one or more loops,
    smoothing each curve separately.
  }
  \item{y}{
    Function value array (the values to be fit).  If the functional data
    are univariate, this array will be an N by NCURVE matrix, where N is
    the number of observed curve values for each curve and NCURVE is the
    number of curves observed.  If the functional data are muliivariate,
    this array will be an N by NCURVE by NVAR matrix, where NVAR the
    number of functions observed per case.  For example, for the gait
    data, NVAR = 2, since we observe knee and hip angles.
  }
  \item{WfdParobj}{
    A functional parameter or fdPar object.  This object contains the
    specifications for the functional data object to be estimated by
    smoothing the data.  See comment lines in function fdPar for
    details.  The functional data object WFD in WFDPAROBJ is used to
    initialize the optimization process.  Its coefficient array contains
    the starting values for the iterative minimization of mean squared
    error.
  }
  \item{wtvec}{
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
  an object of class \code{posfd}, being a list with 4 components:

  \item{Wfdobj}{
    a functional data object defining function $W(x)$ that that
    optimizes the fit to the data of the positive function that it
    defines.
  }
  \item{Flist}{
    a named list containing three results for the final converged
    solution: (1) \bold{f}: the optimal function value being minimized,
    (2) \bold{grad}: the gradient vector at the optimal solution, and
    (3) \bold{norm}: the norm of the gradient vector at the optimal
    solution.
  }
  \item{argvals, y}{
    the corresponding input arguments
  }
}
\seealso{
  \code{\link{smooth.monotone}},
  \code{\link{smooth.morph}}
}
\examples{
smallbasis  <- create.fourier.basis(c(0, 365), 65)
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

index <- (1:35)[CanadianWeather$place == "Vancouver"]
VanPrec  <- CanadianWeather$dailyAv[,index, "Precipitation.mm"]

lambda    <- 1e4
dayfdPar <- fdPar(smallbasis, harmaccelLfd365, lambda)
smooth.pos(day.5, VanPrec, dayfdPar)

}
% docclass is function
\keyword{smooth}
