\name{eval.posfd}
\alias{eval.posfd}
\title{
  Evaluate a Positive Functional Data Object
}
\description{
Evaluate a positive functional data object at specified argument values,
or evaluate a derivative of the functional object.
}
\usage{
eval.posfd(evalarg, Wfdobj, Lfdobj=int2Lfd(0))
}
\arguments{
\item{evalarg}{
a vector of argument values at which the functional data object is to be
evaluated.
}
\item{Wfdobj}{
a functional data object that defines the positive function to be
evaluated.  Only univariate functions are permitted.
}
\item{Lfdobj}{
a nonnegative integer specifying a derivative to be evaluated.  AT
this time of writing, permissible derivative values are 0, 1 or 2.
A linear differential operator is not allowed.
}
}
\value{
a matrix containing the positive function
values.  The first dimension corresponds to the argument values in
\code{evalarg} and
the second to replications.
}
\details{
A positive function data object $h(t)$ is defined by $h(t) =[exp Wfd](t)$.
The function \code{Wfdobj} that defines the positive function is
usually estimated by positive smoothing function
\code{smooth.positive}
}
\seealso{
\code{\link{eval.fd}}, 
\code{\link{eval.monfd}}
}
% docclass is function
\keyword{smooth}
