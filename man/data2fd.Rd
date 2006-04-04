\name{data2fd}
\alias{data2fd}
\title{
  Convert Discrete Data to Functional Data
}
\description{
This function converts an array \code{y} of function values
plus an array \code{argvals} of argument values into a
functional data object.  This a function that tries to
do as much for the user as possible.  A basis function
expansion is used to represent the curve, but no roughness
penalty is used.  The data are fit using the least squares
fitting criterion.
}
\usage{
data2fd(y, argvals=seq(0, 1, len = n), basisobj,
        fdnames=defaultnames,
        argnames=c("time", "reps", "values"))
}
\arguments{
\item{y}{
an array containing sampled values of curves.
If \code{y} is a vector, only one replicate and variable are assumed.
If \code{y} is a matrix, rows must correspond to argument
values and columns to replications or cases, and it will be assumed
that there is only one variable per observation.  If \code{y}
is a three-dimensional array, the first dimension (rows)
corresponds to argument values, the second (columns) to replications,
and the third (layers) to variables within replications.  Missing values
are permitted, and the number of values may vary from one
replication to another.  If this is the case, the number of
rows must equal the maximum number of argument values,
and columns of \code{y} having fewer values must be padded out with
NA's.   
}
\item{argvals}{
a set of argument values.  If this is a vector, the same set of
argument values is used for all columns of \code{y}.  If \code{argvals} is a
matrix, the columns correspond to the columns of \code{y}, and contain the
argument values for that replicate or case.
}
\item{basisobj}{
either:  A \code{basisfd} object created by function create.basis.fd(),
or the value NULL, in which case a \code{basisfd} object is set up by
the function, using the values of the next three arguments.
}
\item{fdnames}{
A list of length 3, each member being a string vector containing labels for the
levels of the corresponding dimension of the discrete data.  The first dimension
is for argument values, and is given the default name "time", the second is
for replications, and is given the default name "reps", and the
third is for
functions, and is given the default name "values".  These default names
are
assigned in function \{tt data2fd\}, which also assigns default string vectors by
using the dimnames attribute of the discrete data array.
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
\code{\link{smooth.basis}}, 
\code{\link{project.basis}}, 
\code{\link{smooth.fd}}, 
\code{\link{smooth.monotone}}, 
\code{\link{smooth.pos}}
}
\examples{

#    set up the fourier basis
nbasis   <- 65
dayrange <- c(0,365)
daybasis <- create.fourier.basis(dayrange, nbasis)
#  Make temperature fd object
#  Temperature data are in 12 by 365 matrix tempav
#    See analyses of weather data.
#  Set up sampling points at mid days
daytime  <- (1:365)-0.5
#  Convert the data to a functional data object
tempfd <- data2fd(daily$tempav, daytime, daybasis)
#  plot the temperature curves
plot(tempfd)

}
% docclass is function
\keyword{smooth}
