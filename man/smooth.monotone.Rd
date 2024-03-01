\name{smooth.monotone}
\alias{smooth.monotone}
\title{
  Monotone Smoothing of Data
}
\description{
  When the discrete data that are observed reflect a smooth strictly
  increasing or strictly decreasing function, it is often desirable to
  smooth the data with a strictly monotone function, even though the
  data themselves may not be monotone due to observational error.  An
  example is when data are collected on the size of a growing organism
  over time.  This function computes such a smoothing function, but,
  unlike other smoothing functions, for only for one curve at a time.
  The smoothing function minimizes a weighted error sum of squares
  criterion.  This minimization requires iteration, and therefore is
  more computationally intensive than normal smoothing.

  The monotone smooth is beta[1]+beta[2]*integral(exp(Wfdobj)), where
  Wfdobj is a functional data object.  Since exp(Wfdobj)>0, its integral
  is monotonically increasing.
}
\usage{
smooth.monotone(argvals, y, WfdParobj, wtvec=rep(1,n),
                zmat=NULL, conv=.0001, iterlim=50,
                active=rep(TRUE,nbasis), dbglev=1)
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
    a vector of data values.  This function can only smooth one set of
    data at a time.

    Function value array (the values to be fit).  If the functional data
    are univariate, this array will be an N by NCURVE matrix, where N is
    the number of observed curve values for each curve and NCURVE is the
    number of curves observed.  If the functional data are multivariate,
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
  \item{zmat}{
    a design matrix or a matrix of covariate values that also define the
    smooth of the data.
  }
  \item{conv}{
    a convergence criterion.
  }
  \item{iterlim}{
    the maximum number of iterations allowed in the minimization of
    error sum of squares.
  }
  \item{active}{
    a logical vector specifying which coefficients defining W(t) are
    estimated.  Normally, the first coefficient is fixed.
  }
  \item{dbglev}{
    either 0, 1, or 2.  This controls the amount information printed out
    on each iteration, with 0 implying no output, 1 intermediate output
    level, and 2 full output.  If either level 1 or 2 is specified, it
    can be helpful to turn off the output buffering feature of S-PLUS.
  }
}

\details{
  The smoothing function  f(argvals) is determined by three objects
  that need to be estimated from the data:

  \itemize{
    \item W(argvals), a functional data object that is first
        exponentiated and then the result integrated.  This is the heart
	of the monotone smooth.  The closer W(argvals) is to zero, the
	closer the monotone smooth becomes a straight line.  The closer
	W(argvals) becomes a constant, the more the monotone smoother
	becomes an exponential function.  It is assumed that W(0) = 0.

	\item b0, an intercept term that determines the value of the
	smoothing function at argvals = 0.

	\item b1, a regression coefficient that determines the slope
	of the smoothing function at argvals = 0.
  }

  In addition, it is possible to have the intercept b0
  depend in turn on the values of one or more covariates through the
  design matrix \code{Zmat} as follows:
  b0 = Z c. In this case, the single
  intercept coefficient is replaced by the regression coefficients
  in vector c multiplying the design matrix.
}
\value{
  an object of class \code{monfd}, which is a list with the following 5
  components:

  \item{Wfdobj}{
    a functional data object defining function W(argvals) that
    optimizes the fit to the data of the monotone function that it
    defines.
  }
  \item{beta}{
    The regression coefficients b_0 and b_1 for each smoothed
    curve.

    If the curves are univariate and
    ... ZMAT is NULL,       BETA is 2   by NCURVE.
    ... ZMAT has P columns, BETA is P+1 by NCURVE.

    If the curves are multivariate and
    ... ZMAT is NULL,       BETA is 2   by NCURVE by NVAR.
    ... ZMAT has P columns, BETA is P+1 by NCURVE by NVAR.
  }
  \item{yhatfd}{
    A functional data object for the monotone curves that smooth the
    data.  This object is constructed using the basis for WFDOBJ, and
    this basis may  well be too simple to accommodate the curvature
    in the monotone function that Wfdobjnes.  It may be
    necessary to discard this object and use a richer basis
    externally to smooth the values defined by
    beta[1] + beta[2]*eval.monfd(evalarg, Wfdobj).
  }
  \item{Flist}{
    a named list containing three results for the final converged
    solution:  (1) \bold{f}: the optimal function value being minimized,
    (2) \bold{grad}: the gradient vector at the optimal solution, and
    (3) \bold{norm}: the norm of the gradient vector at the optimal
    solution.
  }
  \item{y2cMap}{
    For each estimated curve (and variable if functions are
    multivariate, this is an N by NBASIS matrix containing a linear
    mapping from data to coefficients that can be used for computing
    point-wise confidence intervals.  If NCURVE = NVAR = 1, a matrix is
    returned.  Otherwise an NCURVE by NVAR list is returned, with each
    slot containing this mapping.
  }
  \item{argvals}{
    input \code{argvals}, possibly modified / clarified by
    \code{argcheck}.
  }
  \item{y}{
    input argument \code{y}, possibly modified / clarified by
    \code{ycheck}.
  }
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
  \code{\link{smooth.basis}},
  \code{\link{smooth.pos}},
  \code{\link{smooth.morph}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
#  Estimate the acceleration functions for growth curves
#  See the analyses of the growth data.
#  Set up the ages of height measurements for Berkeley data

age <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))
#  Range of observations
rng <- c(1,18)
#  First set up a basis for monotone smooth
#  We use b-spline basis functions of order 6
#  Knots are positioned at the ages of observation.
norder <- 6
nage   <- length(age)
nbasis <- nage + norder - 2
wbasis <- create.bspline.basis(rng, nbasis, norder, age)
#  starting values for coefficient
cvec0 <- matrix(0,nbasis,1)
Wfd0  <- fd(cvec0, wbasis)
#  set up functional parameter object
Lfdobj    <- 3          #  penalize curvature of acceleration
lambda    <- 10^(-0.5)  #  smoothing parameter
growfdPar <- fdPar(Wfd0, Lfdobj, lambda)
#  Set up wgt vector
wgt   <- rep(1,nage)
#  Smooth the data for the first girl
hgt1 = growth$hgtf[,1]

# conv=0.1 to reduce the compute time,
# required to reduce the test time on CRAN

# delete the test on CRAN because it takes too long

if (!CRAN()) {
result <- smooth.monotone(age, hgt1, growfdPar, wgt,
                          conv=0.1)
#  Extract the functional data object and regression
#  coefficients
Wfd  <- result$Wfdobj
beta <- result$beta
#  Evaluate the fitted height curve over a fine mesh
agefine <- seq(1,18,len=73)
hgtfine <- beta[1] + beta[2]*eval.monfd(agefine, Wfd)
#  Plot the data and the curve
plot(age, hgt1, type="p")
lines(agefine, hgtfine)
#  Evaluate the acceleration curve
accfine <- beta[2]*eval.monfd(agefine, Wfd, 2)
#  Plot the acceleration curve
plot(agefine, accfine, type="l")
lines(c(1,18),c(0,0),lty=4)
}
par(oldpar)
}
\keyword{smooth}
