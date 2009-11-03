\name{Data2fd}
\alias{Data2fd}
\title{
  Create a functional data object from data
}
\description{
  This function converts an array \code{y} of function values plus an
  array \code{argvals} of argument values into a functional data object.
  This function tries to do as much for the user as possible.  NOTE:
  Interpolation with data2fd(...) can be shockingly bad, as illustrated
  in one of the examples.
}
\usage{
Data2fd(argvals=NULL, y=NULL, basisobj=NULL, nderiv=NULL,
        lambda=0, fdnames=NULL)
}
\arguments{
  \item{argvals}{
    a set of argument values.  If this is a vector, the same set of
    argument values is used for all columns of \code{y}.  If
    \code{argvals} is a matrix, the columns correspond to the columns of
    \code{y}, and contain the argument values for that replicate or
    case.

    Dimensions for \code{argvals} must match the first dimensions of
    \code{y}, though \code{y} can have more dimensions.  For example, if
    dim(y) = c(9, 5, 2), \code{argvals} can be a vector of length 9 or a
    matrix of dimenions c(9, 5) or an array of dimensions c(9, 5, 2).
  }
  \item{y}{
    an array containing sampled values of curves.

    If \code{y} is a vector, only one replicate and variable are
    assumed.  If \code{y} is a matrix, rows must correspond to argument
    values and columns to replications or cases, and it will be assumed
    that there is only one variable per observation.  If \code{y} is a
    three-dimensional array, the first dimension (rows) corresponds to
    argument values, the second (columns) to replications, and the third
    (layers) to variables within replications.  Missing values are
    permitted, and the number of values may vary from one replication to
    another.  If this is the case, the number of rows must equal the
    maximum number of argument values, and columns of \code{y} having
    fewer values must be padded out with NA's.
  }
  \item{basisobj}{
    One of the following:

    \itemize{
      \item{basisfd}{
	a functional basis object (class \code{basisfd}).
      }
      \item{fd}{
	a functional data object (class \code{fd}), from which its
	\code{basis} component is extracted.
      }
      \item{fdPar}{
	a functional parameter object (class \code{fdPar}), from which
	its \code{basis} component is extracted.
      }
      \item{integer}{
	an integer giving the order of a B-spline basis,
	create.bspline.basis(argvals, norder=basisobj)
      }
      \item{numeric vector}{
	specifying the knots for a B-spline basis,
	create.bspline.basis(basisobj)
      }
      \item{NULL}{
	Defaults to create.bspline.basis(argvals).
      }
    }
  }
  \item{nderiv}{
    Smoothing typically specified as an integer order for the derivative
    whose square is integrated and weighted by \code{lambda} to smooth.
    By default, if basisobj[['type']] == 'bspline', the smoothing
    operator is int2Lfd(max(0, norder-2)).

    A general linear differential operator can also be supplied.
  }
  \item{lambda}{
    weight on the smoothing operator specified by \code{nderiv}.
  }
  \item{fdnames}{
    Either a charater vector of length 3 or a named list of length 3.
    In either case, the three elements correspond to the following:

    \itemize{
      \item{argname}{
	name of the argument, e.g. "time" or "age".
      }
      \item{repname}{
	a description of the cases, e.g. "reps" or "weather stations"
      }
      \item{value}{
	the name of the observed function value, e.g. "temperature"
      }
    }

    If fdnames is a list, the components provide labels for the levels
    of the corresponding dimension of \code{y}.
  }
}
\details{
  This function tends to be used in rather simple applications where
  there is no need to control the roughness of the resulting curve with
  any great finesse.  The roughness is essentially controlled by how
  many basis functions are used.  In more sophisticated applications, it
  would be better to use the function \code{\link{smooth.basisPar}}.
}
\value{
  an object of the \code{fd} class containing:

  \item{coefs}{
    the coefficient array
  }
  \item{basis}{
    a basis object
  }
  \item{fdnames}{
    a list containing names for the arguments, function values
    and variables
  }
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{smooth.basisPar}},
  \code{\link{smooth.basis}},
  \code{\link{project.basis}},
  \code{\link{smooth.fd}},
  \code{\link{smooth.monotone}},
  \code{\link{smooth.pos}}
  \code{\link{day.5}}
}
\examples{
##
## Simplest possible example:  step function
##
b1.1 <- create.bspline.basis(nbasis=1, norder=1)
# 1 basis, order 1 = degree 0 = step function

y12 <- 1:2
fd1.1 <- Data2fd(y12, basisobj=b1.1)
plot(fd1.1)
# fd1.1 = mean(y12) = 1.5

fd1.1.5 <- Data2fd(y12, basisobj=b1.1, lambda=0.5)
eval.fd(seq(0, 1, .2), fd1.1.5)
# fd1.1.5 = sum(y12)/(n+lambda*integral(over arg=0 to 1 of 1))
#         = 3 / (2+0.5) = 1.2

##
## 3 step functions
##
b1.2 <- create.bspline.basis(nbasis=2, norder=1)
# 2 bases, order 1 = degree 0 = step functions
fd1.2 <- Data2fd(1:2, basisobj=b1.2)

op <- par(mfrow=c(2,1))
plot(b1.2, main='bases')
plot(fd1.2, main='fit')
par(op)
# A step function:  1 to 0.5, then 2

##
## Simple oversmoothing
##
b1.3 <- create.bspline.basis(nbasis=3, norder=1)
fd1.3.5 <- Data2fd(y12, basisobj=b1.3, lambda=0.5)
plot(0:1, c(0, 2), type='n')
points(0:1, y12)
lines(fd1.3.5)
# Fit = penalized least squares with penalty =
#          = lambda * integral(0:1 of basis^2),
#            which shrinks the points towards 0.
# X1.3 = matrix(c(1,0, 0,0, 0,1), 2)
# XtX = crossprod(X1.3) = diag(c(1, 0, 1))
# penmat = diag(3)/3
#        = 3x3 matrix of integral(over arg=0:1 of basis[i]*basis[j])
# Xt.y = crossprod(X1.3, y12) = c(1, 0, 2)
# XtX + lambda*penmat = diag(c(7, 1, 7)/6
# so coef(fd1.3.5) = solve(XtX + lambda*penmat, Xt.y)
#                  = c(6/7, 0, 12/7)

##
## linear spline fit
##
b2.3 <- create.bspline.basis(norder=2, breaks=c(0, .5, 1))
# 3 bases, order 2 = degree 1 =
# continuous, bounded, locally linear

fd2.3 <- Data2fd(0:1, basisobj=b2.3)
round(fd2.3$coefs, 4)
# (0, 0, 1),
# though (0, a, 1) is also a solution for any 'a'
op <- par(mfrow=c(2,1))
plot(b2.3, main='bases')
plot(fd2.3, main='fit')
par(op)

# smoothing?
fd2.3. <- Data2fd(0:1, basisobj=b2.3, lambda=1)
\dontshow{stopifnot(}
all.equal(as.vector(round(fd2.3.$coefs, 4)),
          c(0.0159, -0.2222, 0.8730) )
\dontshow{)}
# The default smoothing with spline of order 2, degree 1
# has nderiv = max(0, norder-2) = 0.
# Direct computations confirm that the optimal B-spline
# weights in this case are the numbers given above.

op <- par(mfrow=c(2,1))
plot(b2.3, main='bases')
plot(fd2.3., main='fit')
par(op)

##
## quadratic spline fit
##
b3.4 <- create.bspline.basis(norder=3, breaks=c(0, .5, 1))
# 4 bases, order 3 = degree 2 =
# continuous, bounded, locally quadratic

fd3.4 <- Data2fd(0:1, basisobj=b3.4)
round(fd3.4$coefs, 4)
# (0, 0, 0, 1),
# but (0, a, b, 1) is also a solution for any 'a' and 'b'
op <- par(mfrow=c(2,1))
plot(b3.4)
plot(fd3.4)
par(op)

#  try smoothing?
fd3.4. <- Data2fd(0:1, basisobj=b3.4, lambda=1)
round(fd3.4.$coef, 4)

op <- par(mfrow=c(2,1))
plot(b3.4)
plot(fd3.4.)
par(op)

##
##  A simple Fourier example
##
gaitbasis3 <- create.fourier.basis(nbasis=3)
# note:  'names' for 3 bases
gaitfd3 <- Data2fd(gait, basisobj=gaitbasis3)
# Note: dimanes for 'coefs' + basis[['names']]
# + 'fdnames'

#    set up the fourier basis
daybasis <- create.fourier.basis(c(0, 365), nbasis=65)
#  Make temperature fd object
#  Temperature data are in 12 by 365 matrix tempav
#    See analyses of weather data.

#  Convert the data to a functional data object
tempfd <- Data2fd(CanadianWeather$dailyAv[,,"Temperature.C"],
                  day.5, daybasis)
#  plot the temperature curves
plot(tempfd)

##
## Terrifying interpolation
##
hgtbasis <- with(growth, create.bspline.basis(range(age),
                                              breaks=age, norder=6))
girl.data2fd <- with(growth, Data2fd(hgtf, age, hgtbasis))
age2 <- with(growth, sort(c(age, (age[-1]+age[-length(age)])/2)))
girlPred <- eval.fd(age2, girl.data2fd)
range(growth$hgtf)
range(growth$hgtf-girlPred[seq(1, by=2, length=31),])
# 5.5e-6 0.028 <
# The predictions are consistently too small
# but by less than 0.05 percent

matplot(age2, girlPred, type="l")
with(growth, matpoints(age, hgtf))
# girl.data2fd fits the data fine but goes berzerk
# between points

# Smooth
girl.data2fd1 <- with(growth, Data2fd(age, hgtf, hgtbasis, lambda=1))
girlPred1 <- eval.fd(age2, girl.data2fd1)

matplot(age2, girlPred1, type="l")
with(growth, matpoints(age, hgtf))

# problems splikes disappear

}

% docclass is function
\keyword{smooth}
