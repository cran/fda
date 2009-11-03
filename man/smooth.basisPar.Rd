\name{smooth.basisPar}
\alias{smooth.basisPar}
\title{
  Smooth Data Using a Directly Specified Roughness Penalty
}
\description{
  Smooth (argvals, y) data with roughness penalty defined by the
  remaining arguments.  This function acts as a wrapper for those who want 
  to bypass the step of setting up a functional parameter object before invoking 
  function \code{smooth.basis}.  This function simply does this setup for the 
  user.   See the help file for functions \code{smooth.basis} and \code{fdPar} 
  for further details, and more complete descriptions of the arguments.
}
\usage{
smooth.basisPar(argvals, y, fdobj=NULL, Lfdobj=NULL,
      lambda=0, estimate=TRUE, penmat=NULL,
      wtvec=NULL, fdnames=NULL, covariates=NULL, 
                         method="chol", dfscale=1)
}
\arguments{
  \item{argvals}{
    a set of argument values corresponding to the observations in array
    \code{y}.  In most applications these values will be common to all curves
    and all variables, and therefore be defined as a vector or as a matrix
    with a single column.  But it is possible that these argument values
    will vary from one curve to another, and in this case \code{argvals} will
    be input as a matrix with rows corresponding to observation points and
    columns corresponding to curves.  Argument values can even vary from one
    variable to another, in which case they are input as an array with
    dimensions corresponding to observation points, curves and variables,
    respectively.  Note, however, that the number of observation points per
    curve and per variable may NOT vary.  If it does, then curves and variables
    must be smoothed individually rather than by a single call to this function.
    The default value for \code{argvals} are the integers 1 to \code{n}, where
    \code{n} is the size of the first dimension in argument \code{y}.
  }
  \item{y}{
    an set of values of curves at discrete sampling points or
    argument values. If the set is supplied as a matrix object, the rows must
    correspond to argument values and columns to replications, and it will be
    assumed that there is only one variable per observation.  If
    \code{y} is a three-dimensional array, the first dimension
    corresponds to argument values, the second to replications, and the
    third to variables within replications.  If \code{y} is a vector,
    only one replicate and variable are assumed.  If the data
    come from a single replication but multiple vectors, such as data
    on coordinates for a single space curve, then be sure to coerce
    the data into an array object by using the \code{as.array} function
    with one as the central dimension length.
  }
  \item{fdobj}{
    One of the following:

    \itemize{
      \item{fd}{
	a functional data object (class \code{fd})
      }
      \item{basisfd}{
	a functional basis object (class \code{basisfd}), which is
	converted to a functional data object with the identity matrix
	as the coefficient matrix.
      }
      \item{fdPar}{
	a functional parameter object (class \code{fdPar})
      }
      \item{integer}{
	a positive integer giving the order of a B-spline basis, which is
	further converted to a functional data object with the identity
	matrix as the coefficient matrix.
      }
      \item{matrix or array}{replaced by fd(fdobj)}
      \item{NULL}{
	Defaults to fdobj = create.bspline.basis(argvals).
      }
    }
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object.

    If \code{NULL}, Lfdobj depends on fdobj[['basis']][['type']]:

    \itemize{
      \item{bspline}{
	Lfdobj <- int2Lfd(max(0, norder-2)), where norder =
	norder(fdobj).
      }
      \item{fourier}{
	Lfdobj = a harmonic acceleration operator:

	\code{Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng)}

	where rng = fdobj[['basis']][['rangeval']].
      }
      \item{anything else}{Lfdobj <- int2Lfd(0)}
    }
  }
  \item{lambda}{
    a nonnegative real number specifying the amount of smoothing
    to be applied to the estimated functional parameter.
  }
  \item{estimate}{
    a logical value:  if \code{TRUE}, the functional parameter is
    estimated, otherwise, it is held fixed.
  }
  \item{penmat}{
    a roughness penalty matrix.  Including this can eliminate the need
    to compute this matrix over and over again in some types of
    calculations.
  }
  \item{wtvec}{
    typically a vector of length \code{n} that is the length of \code{argvals}
    containing weights for the values to be smoothed, However, it may also
    be a symmetric matrix of order \code{n}.  If \code{wtvec} is a vector,
    all values must be positive, and if it is a symmetric matrix, this must
    be positive definite.  Defaults to all weights equal to 1.
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
  \item{covariates}{
    the observed values in \code{y} are assumed to be primarily determined
    the the height of the curve being estimates, but from time to time
    certain values can also be influenced by other known variables.  For
    example, multi-year sets of climate variables may be also determined by
    the presence of absence of an El Nino event, or a volcanic eruption.
    One or more of these covariates can be supplied as an \code{n} by
    \code{p} matrix, where \code{p} is the number of such covariates.  When
    such covariates are available, the smoothing is called "semi-parametric."
    Matrices or arrays of regression coefficients are then estimated that
    define the impacts of each of these covariates for each cueve and each
    variable.
  }
  \item{method}{
    by default the function uses the usual textbook equations for computing
    the coefficients of the basis function expansions.  But, as in regression
    analysis, a price is paid in terms of rounding error for such
    computations since they involved cross-products of  basis function
    values.  Optionally, if \code{method} is set equal to the string "qr",
    the computation uses an algorithm based on the qr-decomposition which
    is more accurate, but will require substantially more computing time
    when \code{n} is large, meaning more than 500 or so.  The default
    is "chol", referring the Choleski decomposition of a symmetric positive
    definite matrix.
  }
  \item{dfscale}{
    the generalized cross-validation or "gcv" criterion that is often used
    to determine the size of the smoothing parameter involves the
    subtraction of an measue of degrees of freedom from \code{n}.  Chong
    Gu has argued that multiplying this degrees of freedom measure by
    a constant slightly greater than 1, such as 1.2, can produce better
    decisions about the level of smoothing to be used.  The default value
    is, however, 1.0.
  }
}
\details{
  1.  if(is.null(fdobj))fdobj <- create.bspline.basis(argvals).  Else
  if(is.integer(fdobj)) fdobj <- create.bspline.basis(argvals, norder =
  fdobj)

  2.  fdPar

  3.  smooth.basis
}
\value{
  The output of a call to \code{smooth.basis}, which is an object of
  class \code{fdSmooth}, being a list of length 8 with the following
  components:

  \item{fd}{
    a functional data object that smooths the data.
  }
  \item{df}{
    a degrees of freedom measure of the smooth
  }
  \item{gcv}{
    the value of the generalized cross-validation or GCV criterion.  If
    there are multiple curves, this is a vector of values, one per
    curve.  If the smooth is multivariate, the result is a matrix of gcv
    values, with columns corresponding to variables.
  }
%  \item{coef}{
%    the coefficient matrix or array for the basis function expansion of
%    the smoothing function
%  }
  \item{SSE}{
    the error sums of squares.  SSE is a vector or a matrix of the same
    size as 'gcv'.
  }
  \item{penmat}{
    the penalty matrix.
  }
  \item{y2cMap}{
    the matrix mapping the data to the coefficients.
  }
  \item{argvals, y}{input arguments}
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{Data2fd}},
  \code{\link{df2lambda}},
  \code{\link{fdPar}},
  \code{\link{lambda2df}},
  \code{\link{lambda2gcv}},
  \code{\link{plot.fd}},
  \code{\link{project.basis}},
  \code{\link{smooth.basis}},
  \code{\link{smooth.fd}},
  \code{\link{smooth.monotone}},
  \code{\link{smooth.pos}}
}
\examples{
##
## simplest call
##
girlGrowthSm <- with(growth, smooth.basisPar(argvals=age, y=hgtf, lambda=0.1))
plot(girlGrowthSm$fd, xlab="age", ylab="height (cm)",
         main="Girls in Berkeley Growth Study" )
plot(deriv(girlGrowthSm$fd), xlab="age", ylab="growth rate (cm / year)",
         main="Girls in Berkeley Growth Study" )
plot(deriv(girlGrowthSm$fd, 2), xlab="age",
        ylab="growth acceleration (cm / year^2)",
        main="Girls in Berkeley Growth Study" )
#  Undersmoothed with lambda = 0

##
## Another simple call
##
lipSm <- smooth.basisPar(liptime, lip, lambda=1e-9)$fd
plot(lipSm)

##
## A third example
##

x <- seq(-1,1,0.02)
y <- x + 3*exp(-6*x^2) + sin(1:101)/2
# sin not rnorm to make it easier to compare
# results across platforms

#  set up a saturated B-spline basis
basisobj101 <- create.bspline.basis(x)
fdParobj101 <- fdPar(basisobj101, 2, lambda=1)
result101   <- smooth.basis(x, y, fdParobj101)

resultP <- smooth.basisPar(argvals=x, y=y, fdobj=basisobj101, lambda=1)
\dontshow{stopifnot(}
all.equal(result101, resultP)
\dontshow{)}
# TRUE

result4 <- smooth.basisPar(argvals=x, y=y, fdobj=4, lambda=1)
\dontshow{stopifnot(}
all.equal(resultP, result4)
\dontshow{)}
# TRUE

result4. <- smooth.basisPar(argvals=x, y=y, lambda=1)
\dontshow{stopifnot(}
all.equal(resultP, result4.)
\dontshow{)}
# TRUE

with(result4, c(df, gcv)) #  display df and gcv measures

result4.4 <- smooth.basisPar(argvals=x, y=y, lambda=1e-4)
with(result4.4, c(df, gcv)) #  display df and gcv measures
# less smoothing, more degrees of freedom, better fit

plot(result4.4)
lines(result4, col='green')
lines(result4$fd, col='green') # same as lines(result4, ...)

##
## fdnames?
##
girlGrow12 <- with(growth, smooth.basisPar(argvals=age, y=hgtf[, 1:2], 
              fdnames=c('age', 'girl', 'height'), lambda=0.1) )
girlGrow12. <- with(growth, smooth.basisPar(argvals=age, y=hgtf[, 1:2],
    fdnames=list(age=age, girl=c('Carol', 'Sally'), value='height'),
    lambda = 0.1) )

##
## Fourier basis with harmonic acceleration operator
##
daybasis65 <- create.fourier.basis(rangeval=c(0, 365), nbasis=65)
daytemp.fdSmooth <- with(CanadianWeather, smooth.basisPar(day.5,
       dailyAv[,,"Temperature.C"],
       daybasis65, fdnames=list("Day", "Station", "Deg C")) )


}
% docclass is function
\keyword{smooth}
