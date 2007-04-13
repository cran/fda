\name{smooth.basisPar}
\alias{smooth.basisPar}
\title{
  Smooth Data Using a Directly Specified Roughness Penalty
}
\description{
  Smooth (argvals, y) data with roughness penalty defined by the
  remaining arguments.   
}
\usage{
smooth.basisPar(argvals, y, fdobj=fd(), Lfdobj=int2Lfd(0), lambda=0,
      estimate=TRUE, penmat=NULL, wtvec=rep(1,n), dffactor=1,
      fdnames=list(NULL, dimnames(y)[[2]], NULL))
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
  \item{fdobj}{
    a functional data object, although the argument may also be a
    functional basis object or another functional parameter object.
    Basis objects are converted to functional data objects with the
    identity matrix as the coefficient matrix.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object 
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
    a vector of the same length as \code{argvals} containing weights for
    the values to be smoothed. 
  }
  \item{dffactor}{
    Chong Gu in his book \emph{Smoothing Spline ANOVA Models} suggests a
    modification of the GCV criterion using a factor modifying the
    effective degrees of freedom of the smooth.  He suggests that values
    like 1.2 are effective at avoiding undersmoothing of the data.  The
    default value of 1 is the classic definition. 
  }
  \item{fdnames}{
    a list of length 3 with members containing the following:
    
    \itemize{
      \item a single name for the argument domain, such as 'Time'
      
      \item a vector of names for the replications or cases
      
      \item a name for the function, or a vector of names if there are
      multiple functions. 
    }
  }
}
\value{
  The output of a call to 'smooth.data', which is a named list of length
  7:  

  \item{fdobj}{
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
  \item{coef}{
    the coefficient matrix or array for the basis function expansion of
    the smoothing function 
  }
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
}
\details{
  1.  dfPar

  2.  smooth.data 
}
\seealso{
  \code{\link{data2fd}}, 
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

#  Shows the effects of three levels of smoothing
#  where the size of the third derivative is penalized.
#  The null space contains quadratic functions.
x <- seq(-1,1,0.02)
y <- x + 3*exp(-6*x^2) + sin(1:101)/2
# sin not rnorm to make it easier to compare
# results across platforms 

#  set up a saturated B-spline basis
basisobj <- create.bspline.basis(c(-1,1),101)
fdParobj <- fdPar(basisobj, 2, lambda=1)
result1.  <- smooth.basis(x, y, fdParobj)

result1 <- smooth.basisPar(argvals=x, y=y,
              fdobj=basisobj, Lfdobj=2, lambda=1)
all.equal(result1, result1.)
# TRUE 

with(result1, c(df, gcv)) #  display df and gcv measures
fdParobj <- fdPar(basisobj, 2, lambda=1e-4)
result2 <- smooth.basisPar(x, y, basisobj, 2, lambda=1e-4)

with(result2, c(df, gcv)) #  display df and gcv measures

result3 <- smooth.basisPar(x, y, basisobj, 2, lambda=0)
with(result3, c(df, gcv)) #  display df and gcv measures

plot(x,y)           # plot the data
lines(result1$fd, lty=2)  #  add heavily penalized smooth
lines(result2$fd, lty=1)  #  add reasonably penalized smooth
lines(result3$fd, lty=3)  #  add smooth without any penalty
legend(-1,3,c("1","0.0001","0"),lty=c(2,1,3))
plotfit.fd(y, x, result2$fd)  # plot data and smooth

}
% docclass is function
\keyword{smooth}
