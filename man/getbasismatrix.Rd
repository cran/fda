\name{getbasismatrix}
\alias{getbasismatrix}
\title{Get Values of Basis Functions }
\description{
Evaluates a set of basis functions for a vector of argument values. 
}
\usage{
getbasismatrix(evalarg, basisfd, Lfd=0)
}
\arguments{
  \item{evalarg }{A vector of argument values at which the basis object is to be evaluated. }
  \item{basisfd }{A basis object. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is evaluated rather than the basis functions themselves. }
}
\details{
}
\value{
A matrix with number of rows equal to the length of evalarg and number of columns equal to the number of basis functions. The matrix contains values of each basis function at each argument value. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{
}
\keyword{}
