\name{inprod}
\alias{inprod}
\title{Inner product of two functional data objects or two basis objects}
\description{
This is perhaps the most important function in the functional data library. Hardly any analysis fails to use inner products in some way, and many employ multiple inner products. While in certain cases these may be computed exactly, this is a more general function that approximates the inner product. The inner product is defined by two derivatives or linear differential operators that are applied to the first two arguments. 
}
\usage{
inprod(fd1, fd2, Lfd1=0, Lfd2=0, JMAX=15, EPS=0.0001)
}
\arguments{
  \item{fd1 }{Either a univariate functional data object or a basis object. }
  \item{fd2}{Either a univariate functional data object or a basis object. }
  \item{Lfd1}{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is used rather than the functions themselves. This operator is applied to the first argument fd1. }
  \item{Lfd2}{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is used rather than the functions themselves. This operator is applied to the second argument fd2. }
  \item{JMAX}{This controls the maximum number of function evaluations that will be used in the numerical approximation method (see below for details). The maximum number is 2**(JMAX-1)+1. The default of 15 seems more than adequate for most purposes. }
  \item{EPS}{This controls when convergence is declared for the approximation method. }
}
\details{
The approximation method is Romberg extrapolation using numerical integration by the trapezoidal rule. At each iteration, the number of values at which the functions are evaluated is doubled, and a polynomial extrapolation method is used to estimate the converged integral values as well as an error tolerance. Convergence is declared when the relative error falls below EPS for all products. The Romberg extrapolation method generally saves at least one and often two iterations relative to un-extrapolated trapezoidal integration. 

Functional data analyses will seldom need to use inprod directly, but code developers should be aware of its pivotal role. Future work may require more sophisticated and specialized numerical integration methods. 

inprod computes the definite integral, but some functions such as smooth.monotone and register.fd also need to compute indefinite integrals. These use the same approximation scheme, but usually require more accuracy, and hence more iterations. 

When one or both arguments are basis objects, they are converted to functional data objects using identity matrices as the coefficient matrices. 

inprod is only called when there is no faster or exact method available. In cases where there is, it has been found that the approximation is good to about four to five significant digits, which is sufficient for most applications. Perhaps surprisingly, in the case of B-splines, the exact method is not appreciably faster, but of course is more accurate. 

inprod calls function eval.fd perhaps thousands of times, so high efficiency for this function and the functions that it calls is important. 

}
\value{
A matrix of inner products. The number of rows is the number of functions or basis functions in argument fd1, and the number of columns is the same thing for argument fd2. 
}
\references{
Press, et, al, Numerical Recipes. 
}
\note{}
\author{}
\seealso{
monfn, mongrad, bsplinepen, fourierpen 
}
\examples{
}
\keyword{}
