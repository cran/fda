\name{inprod.bspline}
\alias{inprod.bspline}
\title{
  Compute Inner Products B-spline Expansions.
}
\description{
Computes the matrix of inner products when both functions
are represented by B-spline expansions and when both
derivatives are integers.  This function is called by function
\code{inprod}, and is not normally used directly.
}
\usage{
inprod.bspline(fdobj1, fdobj2=fdobj1, nderiv1=0, nderiv2=0)
}
\arguments{
  \item{fdobj1}{
    a functional data object having a B-spline basis function
    expansion.
  }
  \item{fdobj2}{
    a second functional data object with a B-spline basis
    function expansion.  By default, this is the same as
    the first argument.
  }
  \item{nderiv1}{
    a nonnegative integer specifying the derivative for the
    first argument.
  }
  \item{nderiv2}{
    a nonnegative integer specifying the derivative for the
    second argument.
  }
}
\value{
  a matrix of inner products with number of rows equal
  to the number of replications of the first argument and
  number of columns equal to the number of replications
  of the second object.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\keyword{smooth}
