\name{data2fd}
\alias{data2fd}
\title{Convert Discrete Data to Functional Data }
\description{
This function converts an array y of function values plus an array argvals of argument values into a functional data object. This a function that tries to do as much for the user as possible. 
}
\usage{
data2fd(y, argvals=seq(0, 1, len = n), basisfd,
        fdnames=defaultnames,
        argnames=c("time", "reps", "values"))
}
\arguments{
  \item{y}{An array containing sampled values of curves. If y is a vector, only one replicate and variable are assumed. If y is a matrix, rows must correspond to argument values and columns to replications or cases, and it will be assumed that there is only one variable per observation. If y is a three-dimensional array, the first dimension (rows) corresponds to argument values, the second (columns) to replications, and the third (layers) to variables within replications. Missing values are permitted, and the number of values may vary from one replication to another. If this is the case, the number of rows must equal the maximum number of argument values, and columns of y having fewer values must be padded out with NA's. }
  \item{argvals }{A set of argument values. If this is a vector, the same set of argument values is used for all columns of y. If argvals is a matrix, the columns correspond to the columns of y, and contain the argument values for that replicate or case. }
  \item{basisfd }{Either: A basisfd object created by function create.basis.fd(), or the value NULL, in which case a basisfd object is set up by the function, using the values of the next three arguments. }
  \item{fdnames }{A list of length 3, each member being a string vector containing labels for the levels of the corresponding dimension of the discrete data. The first dimension is for argument values, and is given the default name "time", the second is for replications, and is given the default name "reps", and the third is for functions, and is given the default name "values". These default names are assigned in function {tt data2fd}, which also assigns default string vectors by using the dimnames attribute of the discrete data array. }
  \item{argnames }{A character vector of length 3 containing: (1) the name of the argument, e.g. "time" or "age" (2) a description of the cases, e.g. "weather stations" (3) the name of the observed function value, e.g. "temperature" These strings are used as names for the members of list fdnames. }
}
\details{
}
\value{
An objecgt of the fd class containing the coefficient array with the name coefs, a basis.fd object with the name basis, and a list with the name fdnames. 
}
\references{
Ramsay, J. O. and Silverman, B. W. 1997. Functional data analysis. Springer-Verlag, New York. 310 pp. 
}
\note{}
\author{}
\seealso{
project.basis, smooth.basis, create.fd 
}
\examples{}
\keyword{}
