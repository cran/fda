\name{create.basis}
\alias{create.basis}
\title{Create a Basis Object }
\description{
Create a basis object, used to define a basis for converting discrete data to a functional data object defined by an expansion in terms of the basis functions that this function defines. This is a general purpose function, and special versions of this function are available for each type of basis that is allowed. Normally one would use one of these instead of this function. 
}
\usage{
create.basis(type, rangeval, nbasis, params=paramvec)
}
\arguments{
  \item{type}{One of the following strings: "bspline", "fourier", "expon", "const", "polyg" and "power".}
  \item{rangeval}{A vector of length 2 containing the initial and final values of the interval over which the functional data object can be evaluated. }
  \item{nbasis}{The number of basis functions. }
  \item{params}{A vector of parameters defining the basis. What these mean varies from basis to basis, and the create functions specific to each basis should be consulted. }
}
\details{
See specific versions of this function. 
}
\value{
A basis object. 
}
\references{}
\note{}
\author{}
\seealso{
create.bspline.basis, create.constant.basis, create.exponential.basis, create.fourier.basis, create.polygonal.basis, create.power.basis. 
}
\examples{}
\keyword{}
