\name{cca.fd}
\alias{cca.fd}
\title{Canonical correlation analysis }
\description{
Carry out a functional canonical correlation analysis with regularization or roughness penalties on the estimated canonical variables. 
}
\usage{
cca.fd(fd, ncan=2, lambda=0.00025, Lfd=2)
}
\arguments{
  \item{fd}{A bivariate functional data object. That is, it is assumed that there are two functions for each replication. The coefficient matrix for this object will be three-dimensional, with the final dimension being of length 2. }
  \item{ncan}{Number of pairs of canonical variates to be found.}
  \item{lambda}{Smoothing or regularization parameter. If lambda is a 2-vector then the first component will be applied to the "x" and the second to the "y" functions. }
  \item{Lfd}{The order of derivative to be penalized if it is a scalar, or the linear differential operator whose value is to be penalized if it is a functional data object. }
}
\details{
}
\value{
An object of class cca.fd with these entries is returned. 
"weight functions" A functional data object for the "ncan" pairs of canonical variate weight functions. Each replicate corresponds to a pair of weight functions. 
"correlations" The complete set of canonical correlations. The number is equal to the number of basis functions. 
"variates" An array of the values taken by the canonical variates, that is, the scores on the variates. This is a 3-way array with first dimension corresponding to replicates, second to the different variates (dimension NCAN) and third (dimension 2) to the "x" and "y" scores. 
}
\references{
Ramsay, J. O. and Silverman, B. W. (1997) Functional Data Analysis. New York: Springer. Canonical correlation analysis is described in Chapter 12. 
}
\note{}
\author{}
\seealso{plot.cca.fd, cca.2fun.fd}
\examples{}
\keyword{}
