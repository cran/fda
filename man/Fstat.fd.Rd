\name{Fstat.fd}
\alias{Fstat.fd}
\title{
  F-statistic for functional linear regression.
}
\description{
  Fstat.fd calculates a pointwise F-statistic for functional linear regression.
}
\usage{ Fstat.fd(y,yhat,argvals=NULL)}
\arguments{
  \item{y}{
    the dependent variable object.  It may be:
    \itemize{
      \item a vector if the dependent variable is scalar.
      \item a functional data object if the dependent variable is functional.
    }
  }
  \item{yhat}{
    The predicted values corresponding to \code{y}. It must be of the same class.
  }
  \item{argvals}{
    If \code{yfdPar} is a functional data object, the points at which to evaluate
    the pointwise F-statistic.
  }
}
\details{
   An F-statistic is calculated as the ratio of residual variance to predicted
  variance.

  If \code{argvals} is not specified and \code{yfdPar} is a \code{fd} object,
  it defaults to 101 equally-spaced points on the range of \code{yfdPar}.
}
\value{
  A list with components

  \item{F}{the calculated pointwise F-statistics.}
  \item{argvals}{
    argument values for evaluating the F-statistic if \code{yfdPar} is
    a functional data object.
  }
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\seealso{
  \code{\link{fRegress}}
  \code{\link{Fstat.fd}}
%  \code{\link{tstat.fd}}
}
%\examples{}
\keyword{smooth}
