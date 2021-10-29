\name{polintmat}
\alias{polintmat}
\title{
Polynomial extrapolation for a converging sequence of one or more values}
\description{
For each value in a set \code{xa}, estimates the converged value.  The process
is called Richardson extrapolation.  It terminates when all of two successive 
estimates are within a tolerance of each other.
}
\usage{
polintmat(xa, ya, x)
}
\arguments{
  \item{xa}{Vector of length equal to the maximum number of iterations.}
  \item{ya}{Array with first dimension number of rows equal to length of \code{xa},
  second dimension length equal number of basis functions, and 
  third dimension of length 1.}
  \item{x}{A scalar value which, in our applications, is always zero.}
}
\details{
This function is called with in functions \code{monfn}, \code{mongrad} and 
\code{monhess}.   This function computes the estimate for a single iteration 
of the convergence iterations.
}
\value{
A un-named list of length two containing these objects:
  \item{y:}{Matrix of order number of argument values by number of basis 
  functions.}
  \item{dy:}{Array with same dimensions as argument \code{ya}.}
}
\references{
  Ramsay, James O., Hooker, G. and Graves, S. (2009), \emph{Functional
    Data Analysis with R and Matlab}, Springer, New York.
    
  Ramsay, James O., and Silverman, 
  Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, 
  Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\author{
  J. O. Ramsay
}
\seealso{
  \code{\link{monfn}},
  \code{\link{mongrad}},
  \code{\link{monhess}}
}
