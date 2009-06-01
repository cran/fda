\name{data2fd.old}
\alias{data2fd}
\title{
  Depricated:  use 'Data2fd'
}
\description{
  This function converts an array \code{y} of function values
  plus an array \code{argvals} of argument values into a
  functional data object.  This a function that tries to
  do as much for the user as possible.  A basis function
  expansion is used to represent the curve, but no roughness
  penalty is used.  The data are fit using the least squares
  fitting criterion.  NOTE:  Interpolation with data2fd(...) can be
  shockingly bad, as illustrated in one of the examples.
}
\usage{
data2fd(y, argvals=seq(0, 1, len = n), basisobj,
        fdnames=defaultnames,
        argnames=c("time", "reps", "values"))
}
\arguments{
  \item{y}{
    an array containing sampled values of curves.

    If \code{y} is a vector, only one replicate and variable are
    assumed.

    If \code{y} is a matrix, rows must correspond to argument values and
    columns to replications or cases, and it will be assumed that there
    is only one variable per observation.

    If \code{y} is a three-dimensional array, the first dimension (rows)
    corresponds to argument values, the second (columns) to
    replications, and the third layers) to variables within
    replications.  Missing values are permitted, and the number of
    values may vary from one replication to another.  If this is the
    case, the number of rows must equal the maximum number of argument
    values, and columns of \code{y} having fewer values must be padded
    out with NA's.
  }
  \item{argvals}{
    a set of argument values.

    If this is a vector, the same set of argument values is used for all
    columns of \code{y}.  If \code{argvals} is a matrix, the columns
    correspond to the columns of \code{y}, and contain the argument
    values for that replicate or case.
  }
  \item{basisobj}{
    either:  A \code{basisfd} object created by function
    create.basis.fd(),
    or the value NULL, in which case a \code{basisfd} object is set up
    by the function, using the values of the next three arguments.
  }
  \item{fdnames}{
    A list of length 3, each member being a string vector containing
    labels for the levels of the corresponding dimension of the discrete
    data.  The first dimension is for argument values, and is given the
    default name "time", the second is for replications, and is given
    the default name "reps", and the third is for functions, and is
    given the default name "values".  These default names are
    assigned in function \{tt data2fd\}, which also assigns default
    string vectors by using the dimnames attribute of the discrete data
    array.
  }
  \item{argnames}{
    a character vector of length 3 containing:

    \itemize{
      \item the name of the argument, e.g. "time" or "age"
      \item a description of the cases, e.g. "weather stations"
      \item the name of the observed function value, e.g. "temperature"
    }

    These strings are used as names for the members of list \code{fdnames}.
  }
}
\value{
  an object of the \code{fd} class containing:

  \item{coefs}{
    the coefficient array
  }
  \item{basis}{
    a basis object and
  }
  \item{fdnames}{
    a list containing names for the arguments, function values
    and variables
  }
}
\details{
  This function tends to be used in rather simple applications where
  there is no need to control the roughness of the resulting curve
  with any great finesse.  The roughness is essentially controlled
  by how many basis functions are used.  In more sophisticated
  applications, it would be better to use the function \code{\link{smooth.basis}}
}
\seealso{
  \code{\link{Data2fd}}
  \code{\link{smooth.basis}},
  \code{\link{smooth.basisPar}},
  \code{\link{project.basis}},
  \code{\link{smooth.fd}},
  \code{\link{smooth.monotone}},
  \code{\link{smooth.pos}}
  \code{\link{day.5}}
}
\examples{
# Simplest possible example
b1.2 <- create.bspline.basis(norder=1, breaks=c(0, .5, 1))
# 2 bases, order 1 = degree 0 = step functions

str(fd1.2 <- data2fd(0:1, basisobj=b1.2))
plot(fd1.2)
# A step function:  0 to time=0.5, then 1 after


b2.3 <- create.bspline.basis(norder=2, breaks=c(0, .5, 1))
# 3 bases, order 2 = degree 1 =
# continuous, bounded, locally linear

str(fd2.3 <- data2fd(0:1, basisobj=b2.3))
round(fd2.3$coefs, 4)
# 0, -.25, 1
plot(fd2.3)
# Officially acceptable but crazy:
# Initial negative slope from (0,0) to (0.5, -0.25),
# then positive slope to (1,1).


b3.4 <- create.bspline.basis(norder=3, breaks=c(0, .5, 1))
# 4 bases, order 3 = degree 2 =
# continuous, bounded, locally quadratic

str(fd3.4 <- data2fd(0:1, basisobj=b3.4))
round(fd3.4$coefs, 4)
# 0, .25, -.5, 1
plot(fd3.4)
# Officially acceptable but crazy:
# Initial positive then swings negative
# between 0.4 and ~0.75 before becoming positive again
# with a steep slope running to (1,1).



#  Simple example
gaitbasis3 <- create.fourier.basis(nbasis=3)
str(gaitbasis3) # note:  'names' for 3 bases
gaitfd3 <- data2fd(gait, basisobj=gaitbasis3)
str(gaitfd3)
# Note: dimanes for 'coefs' + basis[['names']]
# + 'fdnames'

#    set up the fourier basis
daybasis <- create.fourier.basis(c(0, 365), nbasis=65)
#  Make temperature fd object
#  Temperature data are in 12 by 365 matrix tempav
#    See analyses of weather data.

#  Convert the data to a functional data object
tempfd <- data2fd(CanadianWeather$dailyAv[,,"Temperature.C"],
                  day.5, daybasis)
#  plot the temperature curves
plot(tempfd)

# Terrifying interpolation
hgtbasis <- with(growth, create.bspline.basis(range(age),
                                              breaks=age, norder=6))
girl.data2fd <- with(growth, data2fd(hgtf, age, hgtbasis))
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

}
% docclass is function
\keyword{smooth}
