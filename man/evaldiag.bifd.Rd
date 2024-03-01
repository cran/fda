\name{evaldiag.bifd}
\alias{evaldiag.bifd}
\title{
  Evaluate the Diagonal of a Bivariate Functional Data Object
}
\description{
Bivariate function data objects are functions of
two arguments, $f(s,t)$.  It can be useful to evaluate
the function for argument values satisfying $s=t$, such
as evaluating the univariate variance function given the
bivariate function that defines the variance-covariance
function or surface.  A linear differential operator can
be applied to function $f(s,t)$ considered as a univariate
function of either object holding the other object fixed.
}
\usage{
evaldiag.bifd(evalarg, bifdobj, sLfd=int2Lfd(0), tLfd=int2Lfd(0))
}
\arguments{
\item{evalarg}{
a vector of values of $s = t$.
}
\item{bifdobj}{
a bivariate functional data object of the \code{bifd} class.
}
\item{sLfd}{
either a nonnegative integer or a linear differential operator
object.
}
\item{tLfd}{
either a nonnegative integer or a linear differential operator
object.
}
}
\value{
a vector or matrix of diagonal function values.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
\code{\link{var.fd}},
\code{\link{eval.bifd}}
}
% docclass is function
\keyword{smooth}
