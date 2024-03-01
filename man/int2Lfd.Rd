\name{int2Lfd}
\alias{int2Lfd}
\title{
  Convert Integer to Linear Differential Operator
}
\description{
  This function turns an integer specifying an order of a derivative
  into the equivalent linear differential operator object.  It is also
  useful for checking that an object is of the "Lfd" class.
}
\usage{
int2Lfd(m=0)
}
\arguments{
  \item{m}{
    either a nonnegative integer or a linear differential operator
    object.
  }
}
\details{
  Smoothing is achieved by penalizing the integral of the square of the
  derivative of order \code{m} over \code{rangeval}:

  m = 0 penalizes the squared difference from 0 of the function

  1 = penalize the square of the slope or velocity

  2 = penalize the squared acceleration

  3 = penalize the squared rate of change of acceleration

  4 = penalize the squared curvature of acceleration?
}
\value{
  a linear differential operator object of the "Lfd" class that is
  equivalent to the integer argument.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\examples{
# Lfd to penalize the squared acceleration
# typical for smoothing a cubic spline (order 4)
int2Lfd(2)

# Lfd to penalize the curvature of acceleration
# used with splines of order 6
# when it is desired to study velocity and acceleration
int2Lfd(4)
}
\keyword{smooth}
