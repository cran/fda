\name{smooth.basis}
\alias{smooth.basis}
\title{
  Smooth Data with an Indirectly Specified Roughness Penalty
}
\description{
  This is the main function for smoothing data using a roughness
  penalty.  This function controls the nature and degree of smoothing by
  penalyzing a measure of rougness.  Roughness is definable in a wide
  variety of ways using either derivatives or a linear differential
  operator.
}
\usage{
smooth.basis(argvals, y, fdParobj, wtvec=rep(1, length(argvals)),
             fdnames=NULL)
}
\arguments{
  \item{argvals}{
    a vector of argument values correspond to the observations in array
    \code{y}.
  }
  \item{y}{
    an array containing values of curves at discrete sampling points or
    argument values. If the array is a matrix, the rows must correspond
    to argument values and columns to replications, and it will be
    assumed that there is only one variable per observation.  If
    \code{y} is a three-dimensional array, the first dimension
    corresponds to argument values, the second to replications, and the
    third to variables within replications.  If \code{y} is a vector,
    only one replicate and variable are assumed.
  }
  \item{fdParobj}{
    a functional parameter object, a functional data object or a
    functional basis object.  If the object is a functional parameter
    object, then the linear differential operator object and the
    smoothing parameter in this object define the roughness penalty.  If
    the object is a functional data object, the basis within this object
    is used without a roughness penalty, and this is also the case if
    the object is a functional basis object.  In these latter two cases,
    \code{smooth.basis} is essentially the same as \code{data2fd}.
  }
  \item{wtvec}{
    a vector of the same length as \code{argvals} containing weights for
    the values to be smoothed.
  }
  \item{fdnames}{
    a list of length 3 containing character vectors of names for the
    following:

    \itemize{
      \item{args}{
	name for each observation or point in time at which data are
	collected for each 'rep', unit or subject.
      }
      \item{reps}{
	name for each 'rep', unit or subject.
      }
      \item{fun}{
	name for each 'fun' or (response) variable measured repeatedly
	(per 'args') for each 'rep'.
      }
    }
  }
}
\details{
  If the smoothing parameter \code{lambda} is zero, there is no penalty
  on roughness.  As lambda increases, usually in logarithmic terms, the
  penalty on roughness increases and the fitted curves become more and
  more smooth.  Ultimately, the curves are forced to have zero roughness
  in the sense of being in the null space of the linear differential
  operator object \code{Lfdobj}that is a member of the \code{fdParobj}.

  For example, a common choice of roughness penalty is the integrated
  square of the second derivative.  This penalizes curvature.  Since the
  second derivative of a straight line is zero, very large values of
  \code{lambda} will force the fit to become linear.  It is also
  possible to control the amount of roughness by using a degrees of
  freedom measure.  The value equivalent to \code{lambda} is found in
  the list returned by the function.  On the other hand, it is possible
  to specify a degrees of freedom value, and then use function
  \code{df2lambda} to determine the equivalent value of \code{lambda}.
  One should not put complete faith in any automatic method for
  selecting \code{lambda}, including the GCV method. There are many
  reasons for this.  For example, if derivatives are required, then the
  smoothing level that is automatically selected may give unacceptably
  rough derivatives.  These methods are also highly sensitive to the
  assumption of independent errors, which is usually dubious with
  functional data.  The best advice is to start with the value
  minimizing the \code{gcv} measure, and then explore \code{lambda}
  values a few log units up and down from this value to see what the
  smoothing function and its derivatives look like.  The function
  \code{plotfit.fd} was designed for this purpose.

  An alternative to using \code{smooth.basis} is to first represent
  the data in a basis system with reasonably high resolution using
  \code{data2fd}, and then smooth the resulting functional data object
  using function \code{smooth.fd}.
}
\value{
  an object of class \code{fdSmooth}, which is a named list of length 8
  with the following components:

  \item{fd}{
    a functional data object containing a smooth of the data.
  }
  \item{df}{
    a degrees of freedom measure of the smooth
  }
  \item{gcv}{
    the value of the generalized cross-validation or GCV criterion.  If
    there are multiple curves, this is a vector of values, one per
    curve.  If the smooth is multivariate, the result is a matrix of gcv
    values, with columns corresponding to variables.

    \deqn{gcv = n*SSE/((n-df)^2)}
  }
%  \item{coef}{
%    the coefficient matrix or array for the basis function expansion of
%    the smoothing function
%  }
  \item{SSE}{
    the error sums of squares.  SSE is a vector or a matrix of the same
    size as GCV.
  }
  \item{penmat}{
    the penalty matrix.
  }
  \item{y2cMap}{
    the matrix mapping the data to the coefficients.
  }
  \item{argvals, y}{input arguments}
}
\seealso{
  \code{\link{data2fd}}, \code{\link{df2lambda}},
  \code{\link{lambda2df}}, \code{\link{lambda2gcv}},
  \code{\link{plot.fd}}, \code{\link{project.basis}},
  \code{\link{smooth.fd}}, \code{\link{smooth.monotone}},
  \code{\link{smooth.pos}}, \code{\link{smooth.basisPar}}
}
\examples{
##
## Example 1:  Inappropriate smoothing
##
# A toy example that creates problems with
# data2fd:  (0,0) -> (0.5, -0.25) -> (1,1)
b2.3 <- create.bspline.basis(norder=2, breaks=c(0, .5, 1))
# 3 bases, order 2 = degree 1 =
# continuous, bounded, locally linear
fdPar2 <- fdPar(b2.3, Lfdobj=2, lambda=1)

\dontrun{
# Penalize excessive slope Lfdobj=1;
# second derivative Lfdobj=2 is discontinuous,
# so the following generates an error:
  fd2.3s0 <- smooth.basis(0:1, 0:1, fdPar2)
Derivative of order 2 cannot be taken for B-spline of order 2
Probable cause is a value of the nbasis argument
 in function create.basis.fd that is too small.
Error in bsplinepen(basisobj, Lfdobj, rng) :
}

##
## Example 2.  Better
##
b3.4 <- create.bspline.basis(norder=3, breaks=c(0, .5, 1))
# 4 bases, order 3 = degree 2 =
# continuous, bounded, locally quadratic
fdPar3 <- fdPar(b3.4, lambda=1)
# Penalize excessive slope Lfdobj=1;
# second derivative Lfdobj=2 is discontinuous.
fd3.4s0 <- smooth.basis(0:1, 0:1, fdPar3)

plot(fd3.4s0$fd)
# same plot via plot.fdSmooth
plot(fd3.4s0)

##
## Example 3.  lambda = 1, 0.0001, 0
##
#  Shows the effects of three levels of smoothing
#  where the size of the third derivative is penalized.
#  The null space contains quadratic functions.
x <- seq(-1,1,0.02)
y <- x + 3*exp(-6*x^2) + rnorm(rep(1,101))*0.2
#  set up a saturated B-spline basis
basisobj <- create.bspline.basis(c(-1,1), 101)

fdPar1 <- fdPar(basisobj, 2, lambda=1)
result1  <- smooth.basis(x, y, fdPar1)
with(result1, c(df, gcv, SSE))

##
## Example 4.  lambda = 0.0001
##
fdPar.0001 <- fdPar(basisobj, 2, lambda=0.0001)
result2  <- smooth.basis(x, y, fdPar.0001)
with(result2, c(df, gcv, SSE))
# less smoothing, more degrees of freedom,
# smaller gcv, smaller SSE

##
##  Example 5.  lambda = 0
##
fdPar0 <- fdPar(basisobj, 2, lambda=0)
result3  <- smooth.basis(x, y, fdPar0)
with(result3, c(df, gcv, SSE))
# Saturate fit:  number of observations = nbasis
# with no smoothing, so degrees of freedom = nbasis,
# gcv = Inf indicating overfitting;
# SSE = 0 (to within roundoff error)

plot(x,y)           # plot the data
lines(result1[['fd']], lty=2)  #  add heavily penalized smooth
lines(result2[['fd']], lty=1)  #  add reasonably penalized smooth
lines(result3[['fd']], lty=3)  #  add smooth without any penalty
legend(-1,3,c("1","0.0001","0"),lty=c(2,1,3))

plotfit.fd(y, x, result2[['fd']])  # plot data and smooth

##
## Example 6.  Supersaturated
##
basis104 <- create.bspline.basis(c(-1,1), 104)

fdPar104.0 <- fdPar(basis104, 2, lambda=0)
result104.0  <- smooth.basis(x, y, fdPar104.0)
with(result104.0, c(df, gcv, SSE))

plotfit.fd(y, x, result104.0[['fd']], nfine=501)
# perfect (over)fit
# Need lambda > 0.

##
## Example 7.  gait
##
gaittime  <- (1:20)/21
gaitrange <- c(0,1)
gaitbasis <- create.fourier.basis(gaitrange,21)
lambda    <- 10^(-11.5)
harmaccelLfd <- vec2Lfd(c(0, 0, (2*pi)^2, 0))

gaitfdPar <- fdPar(gaitbasis, harmaccelLfd, lambda)
gaitfd <- smooth.basis(gaittime, gait, gaitfdPar)$fd
\dontrun{
# by default creates multiple plots, asking for a click between plots
plotfit.fd(gait, gaittime, gaitfd)
}
}
% docclass is function
\keyword{smooth}
