\name{create.fd}
\alias{create.fd}
\title{Create a Functional Data Object }
\description{
Create a functional data object containing functional observations. Note that one would normally do this by a call to the data2fd or smooth.basis functions, so that this function usually not needed, except by when writing new functions for functional data objects. 
}
\usage{
create.fd(coef, basisfd, fdnames=defaultnames)
}
\arguments{
  \item{coef}{A 2- or 3-dimensional array, the first dimension corresponding to basis functions, the second to replications, and the third, if present, to functions.}
  \item{basisfd }{An object of the basis.fd class defining the basis to be used. }
  \item{fdnames}{A list of length 3, each element being a string vector containing labels for the levels of the corresponding dimension of the discrete data. The first dimension is for argument values, and is given the default name time, the second is for replications, and is given the default name reps, and the third is for functions, and is given the default name values. These default names are assigned in function data2fd, which also assigns default string vectors by using the dimnames attribute of the discrete data array. }
}
\details{
}
\value{
Returns a list with the fd class attribute containing the coefficient array with the name coefs, the basis.fd object with the name basis, and the dimension names list with the name fdnames. 
}
\references{}
\note{}
\author{}
\seealso{
basisfdobj, fdobj, create.basis.fd, data2fd 
}
\examples{}
\keyword{}
