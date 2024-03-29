\name{lambda2gcv}
\alias{lambda2gcv}
\title{
  Compute GCV Criterion
}
\description{
  The generalized cross-validation or GCV criterion is often used to
  select an appropriate smoothing parameter value, by finding the
  smoothing parameter that minimizes GCV.  This function locates that
  value.
}
\usage{
  lambda2gcv(log10lambda, argvals, y, fdParobj, wtvec=rep(1,length(argvals)))
}
\arguments{
  \item{log10lambda}{
    the logarithm (base 10) of the smoothing parameter
  }
  \item{argvals}{
    a vector of argument values.
  }
  \item{y}{
    the data to be smoothed.
  }
  \item{fdParobj}{
    a functional parameter object defining the smooth.
  }
  \item{wtvec}{
    a weight vector used in the smoothing.
  }
}
\details{
  Currently, \code{lambda2gcv}
}
\value{
  1.  \eqn{fdParobj[['lambda']] <- 10^log10lambda}

  2.  smoothlist <- smooth.basks(argvals, y, fdParobj, wtvec)

  3.  return(smoothlist[['gcv']])
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
  \code{\link{smooth.basis}}
  \code{\link{fdPar}}
}
% docclass is function
\keyword{smooth}
