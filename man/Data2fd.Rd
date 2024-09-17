\name{Data2fd}
\alias{Data2fd}
\title{
  Create smooth functions that fit scatterplot data.
}
\description{
  The function converts scatter plot data into one or a set of curves that do
  a satisfactory job of representing their corresponding abscissa and ordinate 
  values by smooth curves.  It returns a list object of class \code{fdSmooth} 
  that contains curves defined by functional data objects of class \code{fd}
  as well as a variety of other results related to the data smoothing process.
  
  The function tries to do as much for the user as possible before setting up a 
  call to function \code{smooth.basisPar}. This is achieved by an initial call
  to function \code{argvalsSwap} that examines arguments \code{argvals} that 
  contains abscissa values and \code{y} that contains ordinate values.  
  
  If the arguments are provided in an order that is not the one above, 
  \code{Data2fd} attempts to swap argument positions to provide the correct 
  order. A warning message is returned if this swap takes place.  Any such 
  automatic decision, though, has the possibility of being wrong, and the 
  results should be carefully checked.  

  Preferably, the order of the arguments should be respected: \code{argvals} 
  comes first and \code{y} comes second.
  
  Be warned that the result that the \code{fdSmooth} object that is returned may 
  not define a satisfactory smooth of the data, and consequently that it may be 
  necessary to use the more sophisticated data smoothing function 
  \code{smooth.basis} instead.  Its help file provides a great deal more 
  information than is provided here.}
  
\usage{
  Data2fd(argvals=NULL, y=NULL, basisobj=NULL, nderiv=NULL,
          lambda=3e-8/diff(as.numeric(range(argvals))),
          fdnames=NULL, covariates=NULL, method="chol")}
          
\arguments{
  \item{argvals}{a set of argument values.  If this is a vector, the same set of 
    argument values is used for all columns of \code{y}.  If \code{argvals} is a 
    matrix, the columns correspond to the columns of \code{y}, and contain the 
    argument values for that replicate or case.
    Dimensions for \code{argvals} must match the first dimensions of \code{y}, 
    though \code{y} can have more dimensions.  For example, if dim(y) = 
    c(9, 5, 2), \code{argvals} can be a vector of length 9 or a matrix of 
    dimensions c(9, 5) or an array of dimensions c(9, 5, 2).}
  \item{y}{an array containing sampled values of curves.  If \code{y} is a 
    vector, only one replicate and variable are assumed.  If \code{y} is a 
    matrix, rows must correspond to argument values and columns to replications 
    or cases, and it will be assumed that there is only one variable per 
    observation.  If \code{y} is a three-dimensional array, the first dimension 
    (rows) corresponds to argument values, the second (columns) to replications, 
    and the third (layers) to variables within replications.  Missing values are 
    permitted, and the number of values may vary from one replication to 
    another.  If this is the case, the number of rows must equal the maximum 
    number of argument values, and columns of \code{y} having fewer values 
    must be padded out with NA's.}
  \item{basisobj}{An object of one of the following classes:
    \describe{
      \item{basisfd}{a functional basis object (class \code{basisfd}).}
      \item{fd}{a functional data object (class \code{fd}), from which its 
      \code{basis} component is extracted.}
      \item{fdPar}{a functional parameter object (class \code{fdPar}), from 
        which	its \code{basis} component is extracted.}
      \item{integer}{an integer giving the order of	a B-spline basis, 
      	\code{create.bspline.basis(argvals, norder=basisobj)}}
      \item{numeric vector}{specifying the knots for a B-spline basis, 
      	\code{create.bspline.basis(basisobj)}}
      \item{NULL}{Defaults to create.bspline.basis(argvals).}
    }
  }
  \item{nderiv}{Smoothing typically specified as an integer order for the 
    derivative whose square is integrated and weighted by \code{lambda} to 
    smooth.  By default, \code{if basisobj[['type']] == 'bspline'}, the 
    smoothing operator is \code{int2Lfd(max(0, norder-2)).}

    A general linear differential operator can also be supplied.}
  \item{lambda}{weight on the smoothing operator specified by \code{nderiv}.}
  \item{fdnames}{Either a character vector of length 3 or a named list of length 
    3.  In either case, the three elements correspond to the following:
    \describe{
      \item{argname}{name of the argument, e.g. "time" or "age".}
      \item{repname}{a description of the cases, e.g. "reps" or "weather 
        stations".}
      \item{value}{the name of the observed function value, e.g. "temperature".}
    }
    If fdnames is a list, the components provide labels for the levels of the 
    corresponding dimension of \code{y}.}
  \item{covariates}{the observed values in \code{y} are assumed to be primarily 
    determined by the height of the curve being estimated.  However, from time 
    to time certain values can also be influenced by other known variables.  For 
    example, multi-year sets of climate variables may be also determined by
    the presence of absence of an El Nino event, or a volcanic eruption.One or 
    more of these covariates can be supplied as an \code{n} by\code{p} matrix, 
    where \code{p} is the number of such covariates.  When such covariates are 
    available, the smoothing is called "semi-parametric."  Matrices or arrays of 
    regression coefficients are then estimated that define the impacts of each 
    of these covariates for each curve and each variable.}
  \item{method}{by default the function uses the usual textbook equations for 
    computing the coefficients of the basis function expansions.  But, as in 
    regression analysis, a price is paid in terms of rounding error for such 
    computations since they involved cross-products of  basis function values.  
    Optionally, if \code{method} is set equal to the string "qr", the 
    computation uses an algorithm based on the qr-decomposition which is more 
    accurate, but will require substantially more computing timewhen \code{n} is 
    large, meaning more than 500 or so.  The defaultis "chol", referring the 
    Choleski decomposition of a symmetric positive definite matrix.}
}

\details{
  This function tends to be used in rather simple applications where there is no 
  need to control the roughness of the resulting curve with any great finesse.  
  The roughness is essentially controlled by how many basis functions are used.  
  In more sophisticated applications, it would be better to use the function   
  \code{\link{smooth.basisPar}}.
  
  It may happen that a value in argvals is outside the basis object interval 
  due to rounding error and other causes of small violations.  The code tests
  for this and pulls these near values into the interval if they are within
  1e-7 times the interval width.
}

\value{
  an S3 list object of class \code{fdSmooth} defined in file 
  \code{smooth.basis1.R} containing:
  \describe{
    \item{fd}{the functional data object defining smooth B-spline curves that
      fit the data}
    \item{df}{the degrees of freedom in the fits}
    \item{gcv}{a generalized cross-validation value}
    \item{beta}{the semi-parametric coefficient array}
    \item{SSE}{Error sum of squares}
    \item{penmat}{the symmetric matrix defining roughness penalty}
    \item{argvals}{the argument values}
    \item{y}{the array of ordinate values in the scatter-plot} 
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
  \code{\link{smooth.basisPar}},
  \code{\link{smooth.basis}},
  \code{\link{smooth.basis1}},
  \code{\link{project.basis}},
  \code{\link{smooth.fd}},
  \code{\link{smooth.monotone}},
  \code{\link{smooth.pos}},
  \code{\link{day.5}}
}

\examples{
##
## Simplest possible example:  constant function
##
# 1 basis, order 1 = degree 0 = constant function
b1.1 <- create.bspline.basis(nbasis=1, norder=1)
# data values: 1 and 2, with a mean of 1.5
y12 <- 1:2
# smooth data, giving a constant function with value 1.5
fd1.1 <- Data2fd(y12, basisobj=b1.1)
oldpar <- par(no.readonly=TRUE)
plot(fd1.1)
# now repeat the analysis with some smoothing, which moves the
# toward 0.
fd1.5 <- Data2fd(y12, basisobj=b1.1, lambda=0.5)
#  values of the smooth:
# fd1.1 = sum(y12)/(n+lambda*integral(over arg=0 to 1 of 1))
#         = 3 / (2+0.5) = 1.2
#. JR ... Data2fd returns an fdsmooth object and a member of the
#         is an fd smooth object.  Calls to functions expecting 
#         an fd object require attaching $fd to the fdsmooth object
#         this is required in lines 268, 311 and 337
eval.fd(seq(0, 1, .2), fd1.5)
##
## step function smoothing
##
# 2 step basis functions: order 1 = degree 0 = step functions
b1.2 <- create.bspline.basis(nbasis=2, norder=1)
#  fit the data without smoothing
fd1.2 <- Data2fd(1:2, basisobj=b1.2)
# plot the result:  A step function:  1 to 0.5, then 2
op <- par(mfrow=c(2,1))
plot(b1.2, main='bases')
plot(fd1.2, main='fit')
par(op)
##
## Simple oversmoothing
##
# 3 step basis functions: order 1 = degree 0 = step functions
b1.3 <- create.bspline.basis(nbasis=3, norder=1)
#  smooth the data with smoothing
fd1.3 <- Data2fd(y12, basisobj=b1.3, lambda=0.5)
#  plot the fit along with the points
plot(0:1, c(0, 2), type='n')
points(0:1, y12)
lines(fd1.3)
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
# 3 bases, order 2 = degree 1
b2.3 <- create.bspline.basis(norder=2, breaks=c(0, .5, 1))
# interpolate the values 0, 2, 1
fd2.3 <- Data2fd(c(0,2,1), basisobj=b2.3, lambda=0)
#  display the coefficients
round(fd2.3$coefs, 4)
# plot the results
op <- par(mfrow=c(2,1))
plot(b2.3, main='bases')
plot(fd2.3, main='fit')
par(op)
# apply some smoothing
fd2.3. <- Data2fd(c(0,2,1), basisobj=b2.3, lambda=1)
op <- par(mfrow=c(2,1))
plot(b2.3, main='bases')
plot(fd2.3., main='fit', ylim=c(0,2))
par(op)
all.equal(
 unclass(fd2.3)[-1], 
 unclass(fd2.3.)[-1])
##** CONCLUSION:  
##** The only differences between fd2.3 and fd2.3.
##** are the coefficients, as we would expect.  

##
## quadratic spline fit
##
# 4 bases, order 3 = degree 2 = continuous, bounded, locally quadratic
b3.4 <- create.bspline.basis(norder=3, breaks=c(0, .5, 1))
# fit values c(0,4,2,3) without interpolation
fd3.4 <- Data2fd(c(0,4,2,3), basisobj=b3.4, lambda=0)
round(fd3.4$coefs, 4)
op <- par(mfrow=c(2,1))
plot(b3.4)
plot(fd3.4)
points(c(0,1/3,2/3,1), c(0,4,2,3))
par(op)
#  try smoothing
fd3.4. <- Data2fd(c(0,4,2,3), basisobj=b3.4, lambda=1)
round(fd3.4.$coef, 4)
op <- par(mfrow=c(2,1))
plot(b3.4)
plot(fd3.4., ylim=c(0,4))
points(seq(0,1,len=4), c(0,4,2,3))
par(op)
##
##  Two simple Fourier examples
##
gaitbasis3 <- create.fourier.basis(nbasis=5)
gaitfd3    <- Data2fd(seq(0,1,len=20), gait, basisobj=gaitbasis3)
# plotfit.fd(gait, seq(0,1,len=20), gaitfd3)
#    set up the fourier basis
daybasis <- create.fourier.basis(c(0, 365), nbasis=65)
#  Make temperature fd object
#  Temperature data are in 12 by 365 matrix tempav
#    See analyses of weather data.
tempfd <- Data2fd(CanadianWeather$dailyAv[,,"Temperature.C"],
                  day.5, daybasis)
#  plot the temperature curves
par(mfrow=c(1,1))
plot(tempfd)
##
## argvals of class Date and POSIXct
##
#  These classes of time can generate very large numbers when converted to 
#  numeric vectors.  For basis systems such as polynomials or splines,
#  severe rounding error issues can arise if the time interval for the 
#  data is very large.  To offset this, it is best to normalize the
#  numeric version of the data before analyzing them.
#  Date class time unit is one day, divide by 365.25.
invasion1 <- as.Date('1775-09-04')
invasion2 <- as.Date('1812-07-12')
earlyUS.Canada <- as.numeric(c(invasion1, invasion2))/365.25
BspInvasion <- create.bspline.basis(earlyUS.Canada)
earlyYears  <- seq(invasion1, invasion2, length.out=7)
earlyQuad   <- (as.numeric(earlyYears-invasion1)/365.25)^2
earlyYears  <- as.numeric(earlyYears)/365.25
fitQuad <- Data2fd(earlyYears, earlyQuad, BspInvasion)
# POSIXct: time unit is one second, divide by 365.25*24*60*60
rescale     <- 365.25*24*60*60
AmRev.ct    <- as.POSIXct1970(c('1776-07-04', '1789-04-30'))
BspRev.ct   <- create.bspline.basis(as.numeric(AmRev.ct)/rescale)
AmRevYrs.ct <- seq(AmRev.ct[1], AmRev.ct[2], length.out=14)
AmRevLin.ct <- as.numeric(AmRevYrs.ct-AmRev.ct[1])
AmRevYrs.ct <- as.numeric(AmRevYrs.ct)/rescale
AmRevLin.ct <- as.numeric(AmRevLin.ct)/rescale
fitLin.ct   <- Data2fd(AmRevYrs.ct, AmRevLin.ct, BspRev.ct)
par(oldpar)
}
\keyword{smooth}
