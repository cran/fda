\name{cca.2fun.fd}
\alias{cca.2fun.fd}
\title{Canonical Correlation for Two Functional Data Objects}
\description{
This is an extension of function cca.fd that computes a functional canonical correlation with regularization using two different functional data samples. These may have different bases, and different penalty functionals. It is assumed that both are univariate. Two differential operators are also possible. 
}
\usage{
cca.2fun.fd(fd1, fd2=fd1, ncan=2, lambda1=0.00025, Lfd1=2, lambda2=lambda1, Lfd2=Lfd1)
}
\arguments{
  \item{fd1}{The first functional data object.}
  \item{fd2}{The second functional data object.}
  \item{ncan}{The number of canonical variables to be estimated.}
  \item{lambda1}{The smoothing parameter for the first set of canonical variables.}
  \item{Lfd1}{The order of derivative to be penalized if it is a scalar, or the linear differential operator whose value is to be penalized if it is a functional data object. This applies to the first functional data object fd1. }
  \item{lambda2}{The smoothing parameter for the second set of canonical variables.}
  \item{Lfd2}{The order of derivative to be penalized if it is a scalar, or the linear differential operator whose value is to be penalized if it is a functional data object. This applies to the second functional data object fd2. }
}
\details{
Function cca.fd is mainly intended for situations where the two functions are reflecting the same variable, such as angle in the gait data. However, one also can want to correlate two variables of an entirely different character, such as temperature and precipitation in the weather data. In this case the two functional data objects will be univariate, but otherwise should be defined separately, and possibly with rather different bases. 
}
\value{
A list with these entries is returned. 
"weight function 1 " A functional data object for the first canonical variable in the "ncan" pairs of canonical variate weight functions. Each replicate corresponds to a pair of weight functions. 
"weight function 2 " A functional data object for the first canonical variable in the "ncan" pairs of canonical variate weight functions. Each replicate corresponds to a pair of weight functions. 
"correlations" The complete set of canonical correlations. The number is equal to the number of basis functions. 
"variates" An array of the values taken by the canonical variates, that is, the scores on the variates. This is a 3-way array with first dimension corresponding to replicates, second to the different variates (dimension NCAN) and third (dimension 2) to the "x" and "y" scores. 
}
\references{Ramsay, J. O. and Silverman, B. W. (1997) Functional Data Analysis. New York: Springer. Canonical correlation analysis is described in Chapter 12. 
}
\note{}
\author{}
\seealso{cca.fd}
\examples{}
\keyword{}
