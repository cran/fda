\name{smooth.Pspline}
\alias{smooth.Pspline}
\title{Smooth Data with Using mth Derivative Penalty }
\description{
This a generalized version of one of the regular S-PLUS functions, smooth.spline, that permits a penalty on the size of the derivative of any order m rather than just m = 2. It uses functions coded in C or Fortran that economize on storage and are much faster than the native S-PLUS code that is used in function smooth.basis. Consequently, it is more suitable for problems involving the need for very large number of basis functions, meaning in practice more than about 100. 
}
\usage{
smooth.Pspline(x, y, w, norder=2, df=norder+2, spar, method=3)
}
\arguments{
  \item{x }{Either a vector or matrix containing the data values to be smoothed. If a matrix, the columns correspond to replications of curves or functions. }
  \item{y }{Either a vector or matrix containing the data values to be smoothed. If a matrix, the columns correspond to replications of curves or functions. }
  \item{w }{Either a vector or matrix containing the data values to be smoothed. If a matrix, the columns correspond to replications of curves or functions. }
  \item{norder }{The order of the derivative whose size is to be penalized. See below for advice on what it should be. }
  \item{df }{The degree of freedom. }
  \item{spar }{A penalty parameter controlling the amount of smoothness or the size of the penalty if method = 3 or 4. }
  \item{method }{An integer with a value from 1 to 4 specifying the way in which the amount of smoothing is determined: 
    the amount of smoothing is controlled by the smoothing parameter lambda. 
    the amount of smoothing is controlled by the degrees of freedom value df. 
    the amount of smoothing optimizes the GCV or generalized cross-validation criterion. 
    the amount of smoothing optimizes the CV or ordinary cross-validation criterion. 
}
}
\details{
The need to choose the order m of the derivative arises if one needs smooth estimates of one or more derivatives. Generally, the advice is to choose m to be 2 larger than the highest order of derivative required. This in effect penalizes the curvature of that derive. So if the second derivative or acceleration is required, then use m = 4. 

Like S-PLUS function smooth.spline, this function also comes with additional functions: 

predictSmoothPspline
calculates values of a derivative of the smoothing function at argument values argvals. The call is 
predictSmoothPspline(splobj, argvals, nderiv = 0)

where splobj is the object of class "smooth.Pspline" returned by smooth.Pspline. 
plotSmoothPspline
plots an object of class "smooth.Pspline". The call is 
plotSmoothPspline(splobj, ...)

where splobj is the object of class "smooth.Pspline" returned by smooth.Pspline. Additional arguments for controlling the plot may also supplied, as in the regular plot function. 
linesSmoothPspline
adds lines to an existing plot showing an object of class "smooth.Pspline". The call is 
linesSmoothPspline(splobj, ...)

where splobj is the object of class "smooth.Pspline" returned by smooth.Pspline. Additional arguments for controlling the plot may also supplied, as in the regular lines function. 
print.smooth.Pspline
prints an object of class "smooth.Pspline". The call is 
print.smooth.Pspline(splobj, ...)

where splobj is the object of class "smooth.Pspline" returned by smooth.Pspline. Additional arguments for controlling the plot may also supplied, as in the regular print function. 
One should not put complete faith in any automatic method for selecting lambda, including the GCV method. There are many reasons for this. For example, if derivatives are required, then the smoothing level that is automatically selected may give unacceptably rough derivatives. These methods are also highly sensitive to the assumption of independent errors, which is usually dubious with functional data. The best advice is to start with the value minimizing the gcv measure, and then explore lambda values a few log units up and down from this value to see what the smoothing function and its derivatives look like. 
}
\value{
An object of class "smooth.Pspline" containing: 
norder
The order of the smoothing spline, that is, the derivative whose size is penalized. 
argvals
The argument value vector. 
ysmth
The values of the smooth function(s) corresponding to the values in y. 
lev
An vector of leverage values. 
gcv
The generalized cross-validation coefficient. 
cv
The ordinary cross-validation coefficient. 
lambda
The final value of the penalty parameter. 
my.call
The calling statement. 
}
\references{}
\note{}
\author{}
\seealso{
smooth.basis, smooth.spline, smooth.fd 
}
\examples{}
\keyword{}
