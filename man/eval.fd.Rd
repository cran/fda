\name{eval.fd}
\alias{eval.fd}
\title{Evaluate a Functional Data Object}
\description{
Evaluate a functional data object at specified argument values, or evaluate a derivative or the result of applying a linear differential operator to the functional object. 
}
\usage{
eval.fd(evalarg, fd, Lfd=0)
}
\arguments{
  \item{evalarg }{A vector of argument values at which the functional data object is to be evaluated. }
  \item{fd}{A functional data object to be evaluated. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
}
\details{
}
\value{
Returns an array of 2 or 3 dimensions containing the function values. The first dimension corresponds to the argument values in evalarg, the second to replications, and the third if present to functions. 
}
\references{}
\note{}
\author{}
\seealso{
plotFd, plotfit.fd, fd2data, eval.bifd 
}
\examples{}
\keyword{}
