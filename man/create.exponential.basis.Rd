\name{create.exponential.basis}
\alias{create.exponential.basis}
\title{Create an Exponential Basis}
\description{
Create an exponential basis object defining a set of exponential functions with rate constants in argument ratevec. 
}
\usage{
create.exponential.basis(rangeval=c(0, 1), ratevec=1)
}
\arguments{
  \item{rangeval}{A vector of length 2 containing the initial and final values of the interval over which the functional data object can be evaluated. }
  \item{ratevec }{A vector of length nbasis of rate constants defining basis functions of the form exp(rate*x). }
}
\details{
Exponential functions are of the type ebx where b is the rate constant. If b = 0, the constant function is defined. 
}
\value{
A basis object with the type "expon". 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
