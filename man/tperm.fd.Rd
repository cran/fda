\name{tperm.fd}
\alias{tperm.fd}
\title{
  Permutation t-test for two groups of functional data objects.
}
\description{
  tperm.fd creates a null distribution for a test of no difference between two
  groups of functional data objects.
}
\usage{
 tperm.fd(x1fd, x2fd, nperm=200, q=0.95, argvals=NULL, plotres=TRUE)
}
\arguments{
  \item{x1fd}{
    a functional data object giving the first group of functional observations.
  }
  \item{x2fd}{
    a functional data object giving the second group of functional
    observations.
  }
  \item{nperm}{
    number of permutations to use in creating the null distribution.
  }
  \item{argvals}{
    If \code{yfdPar} is a \code{fd} object, the points at which to evaluate
    the point-wise F-statistic.
  }
  \item{q}{
    Critical quantile of the null distribution to compare to the observed
    F-statistic.
  }
  \item{plotres}{
    Argument to plot a visual display of the null distribution displaying the
    \code{q}th quantile and observed F-statistic.
  }
}
\details{
   The usual t-statistic is calculated pointwise and the test based on the
   maximal value. If \code{argvals} is not specified,
  it defaults to 101 equally-spaced points on the range of \code{yfdPar}.
}
\value{
  A list with components
  \itemize{
    \item{pval}{the observed p-value of the permutation test.}
    \item{qval}{the \code{q}th quantile of the null distribution.}
    \item{Tobs}{the observed maximal t-statistic.}
    \item{Tnull}{
      a vector of length \code{nperm} giving the observed values of the
      permutation distribution.
    }
    \item{Tvals}{the pointwise values of the observed t-statistic.}
    \item{Tnullvals}{the pointwise values of of the permutation observations.}
    \item{pvals.pts}{pointwise p-values of the t-statistic.}
    \item{qvals.pts}{pointwise \code{q}th quantiles of the null distribution}
    \item{argvals}{
      argument values for evaluating the F-statistic if \code{yfdPar}is
      a functional data object.
    }
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
%\examples{ }
\keyword{smooth}
