\name{create.bifd}
\alias{create.bifd}
\title{Create a Two-Argument Functional Data Object}
\description{
Create a bifd object containing functional observations defined over two arguments. This function is not normally needed; it is called within other functions to create two-argument functions such as covariance functions and bivariate regression functions. 
}
\usage{
create.bifd(coef, sbasisfd, tbasisfd, bifdnames=list(NULL, repnames, NULL))
}
\arguments{
  \item{coef}{A 2-, 3-, or 4-dimensional array, the first two dimensions corresponding to basis functions, the third to replications, and the fourth, if present, to functions. }
  \item{sbasisfd}{An object of the basis.fd class for the first argument. }
  \item{tbasisfd}{An object of the basis.fd class for the second argument. }
  \item{bifdnames}{A list of length 3, each element being a string vector containing labels for the levels of the corresponding dimension of the discrete data. The first dimension is for argument values, and is given the default name time, the second is for replications, and is given the default name reps, and the third is for functions, and is given the default name values. These default names are assigned in function data2fd, which also assigns default string vectors by using the dimnames attribute of the discrete data array. }
}
\details{
}
\value{
Returns a list with the bifd class attribute containing the coefficient array with the name coefs, a basis.fd object with the name sbasis, and a basis.fd object with the name tbasis. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
