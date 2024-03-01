\name{expon}
\alias{expon}
\title{
  Exponential Basis Function Values
}
\description{
Evaluates a set of exponential basis functions, or a derivative of these
functions, at a set of arguments.
}
\usage{
expon(x, ratevec=1, nderiv=0)
}
\arguments{
\item{x}{
a vector of values at which the basis functions are to be evaluated.
}
\item{ratevec}{
a vector of rate or time constants defining the exponential
functions.  That is, if $a$ is the value of an
element of this vector, then the corresponding basis function
is $exp(at)$. The number of basis functions is equal
to the length of \code{ratevec}.
}
\item{nderiv}{
a nonnegative integer specifying an order of derivative to
be computed.  The default is 0, or the basis function value.
}
}
\value{
a matrix of basis function values with rows corresponding
to argument values and columns to basis functions.
}
\details{
There are no restrictions on the rate constants.
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
  \code{\link{exponpen}}
}
% docclass is function
\keyword{smooth}
