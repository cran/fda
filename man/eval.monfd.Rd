\name{eval.monfd}
\alias{eval.monfd}
\title{Evaluate a Monotone Functional Data Object }
\description{
Evaluate a monotone functional data object at specified argument values, or evaluate a derivative of the functional object. 
}
\usage{
eval.monfd(evalarg, Wfd, Lfd=0)
}
\arguments{
  \item{evalarg }{A vector of argument values at which the functional data object is to be evaluated. }
  \item{Wfd}{A functional data object that defines the monotone function to be evaluated. Only univariate functions are permitted. }
  \item{Lfd }{A nonnegative integer specifying a derivative to be evaluated. AT this time of writing, permissible derivative values are 0, 1, 2, or 3. A linear differential operator is not allowed. }
}
\details{
A monotone function data object h(t) is defined by h(t) = [D**{-1} exp Wfd](t). In this equation, the operator D**{-1} means taking the indefinite integral of the function to which it applies. Note that this equation implies that the monotone function has a value of zero at the lower limit of the arguments. To actually fit monotone data, it will usually be necessary to estimate an intercept and a regression coefficient to be applied to h(t), usually with the least squares regression function lsfit. 

The function Wfd that defines the monotone function is usually estimated by monotone smoothing function smooth.monotone. 
}
\value{
Returns an matrix containing the monotone function values. The first dimension corresponds to the argument values in evalarg and the second to replications. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
