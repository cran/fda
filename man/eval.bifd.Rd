\name{eval.bifd}
\alias{eval.bifd}
\title{Evaluate a Two-argument Functional Data Object }
\description{
Evaluate a two-argument functional data object at specified argument values s and t. 
}
\usage{
eval.bifd(sevalarg, tevalarg, bifd, sLfd=0, tLfd=0)
}
\arguments{
  \item{sevalarg }{A vector of argument values for the first argument s of the functional data object to be evaluated. }
  \item{tevalarg }{A vector of argument values for the second argument t of the functional data object to be evaluated. }
  \item{bifd }{A two-argument functional data object. }
  \item{sLfd}{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator to the object as a function of the first argument s is evaluated rather than the functions themselves.}
  \item{tLfd }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator to the object as a function of the second argument t is evaluated rather than the functions themselves. }
}
\details{
}
\value{
An array of 2, 3, or 4 dimensions containing the function values. The first dimension corresponds to the argument values in sevalarg, the second to argument values in tevalarg, the third if present to replications, and the fourth if present to functions. 
}
\references{}
\note{}
\author{}
\seealso{
eval.fd 
}
\examples{}
\keyword{}
