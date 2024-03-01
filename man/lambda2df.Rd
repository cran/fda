\name{lambda2df}
\alias{lambda2df}
\title{
  Convert Smoothing Parameter to Degrees of Freedom
}
\description{
The degree of roughness of an estimated function is controlled by a
smoothing parameter $lambda$ that directly multiplies the penalty.
However, it can be difficult to interpret or choose this value, and it
is often easier to determine the roughness by choosing a value that is
equivalent of the degrees of freedom used by the smoothing procedure.
This function converts a multipler $lambda$ into a degrees of freedom value.
}
\usage{
lambda2df(argvals, basisobj, wtvec=rep(1, n), Lfdobj=NULL, lambda=0)
}
\arguments{
  \item{argvals}{
    a vector containing the argument values used in the
    smooth of the data.
  }
  \item{basisobj}{
    the basis object used in the smoothing of the data.
  }
  \item{wtvec}{
    the weight vector, if any, that was used in the smoothing
    of the data.
  }
  \item{Lfdobj}{
    the linear differential operator object used to defining
    the roughness penalty employed in smoothing the data.
  }
  \item{lambda}{
    the smoothing parameter to be converted.
  }
}
\value{
  the equivalent degrees of freedom value.
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
  \code{\link{df2lambda}}
}
% docclass is function
\keyword{smooth}
