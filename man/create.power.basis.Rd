\name{create.power.basis}
\alias{create.power.basis}
\title{Create a Power Basis Object }
\description{
The basis system is a set of powers of argument x. That is, a basis function would be x**exponent, where exponent is the value of one of the exponents. The power basis would normally only be used for positive values of x. 
}
\usage{
create.power.basis(rangeval=c(0, 1), exponent=1)
}
\arguments{
  \item{rangeval }{A vector of length 2 with the first element being the lower limit of the range of argument values, and the second the upper limit. Of course the lower limit must be less than the upper limit. }
  \item{exponent }{A vector of length nbasis containing the powers. }
}
\details{
}
\value{
A basis object used to define power function bases. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
