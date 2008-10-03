\name{Fperm.fd}
\alias{Fperm.fd}
\title{
  Permutation F-test for functional linear regression.
}
\description{
  Fperm.fd creates a null distribution for a test of no effect in functional
  linear regression. It makes generic use of \code{fRegress} and permutes the
  \code{yfdPar} input.
}
\usage{
Fperm.fd(yfdPar, xfdlist, betalist,wt=NULL,
            nperm=200,argvals=NULL,q=0.05,plotres=TRUE,...)
}
\arguments{
  \item{yfdPar}{
    the dependent variable object.  It may be an object of
    three possible classes:
    \itemize{
      \item{vector}{ if the dependent variable is scalar.}
      \item{fd}{
    a functional data object if the dependent variable is
    functional.
      }
      \item{fdPar}{
    a functional parameter object if the dependent variable is
    functional, and if it is necessary to smooth the prediction of
    the dependent variable.
      }
    }
  }
  \item{xfdlist}{
    a list of length equal to the number of independent variables. Members
    of this list are the independent variables.  They be objects of either
    of these two classes:

    \itemize{
      \item a vector if the independent dependent variable is scalar.
      \item a functional data object if the dependent variable is functional.
    }

    In either case, the object must have the same number of replications as
    the dependent variable object.  That is, if it is a scalar, it must be
    of the same length as the dependent variable, and if it is functional,
    it must have the same number of replications as the dependent variable.
  }
  \item{betalist}{
    a list of length equal to the number of independent variables. Members
    of this list define the regression functions to be estimated.
    They are functional parameter objects.  Note that even if corresponding
    independent variable is scalar, its regression coefficient will be
    functional if the dependent variable is functional.  Each of these
    functional parameter objects defines a single functional data object,
    that is, with only one replication.
  }
  \item{wt}{
    weights for weighted least squares, defaults to all 1.
  }
  \item{nperm}{
    number of permutations to use in creating the null distribution.
  }
  \item{argvals}{
    If \code{yfdPar} is a \code{fd} object, the points at which to evaluate
    the point-wise F-statistic.
  }
  \item{q}{
    Critical upper-tail quantile of the null distribution to compare to the observed F-statistic.
  }
  \item{plotres}{
    Argument to plot a visual display of the null distribution displaying the
    \code{q}th quantile and observed F-statistic.
  }
  \item{...}{
    Additional plotting arguments that can be used with \code{plot}.
  }
} \details{
   An F-statistic is calculated as the ratio of residual variance to predicted
  variance. The observed F-statistic is returned along with the permutation
  distribution.

  If \code{yfdPar} is a \code{fd} object, the maximal value of the pointwise
  F-statistic is calculated. The pointwise F-statistics are also returned.

  The default of setting \code{q = 0.95} is, by now, fairly standard. The default
  \code{nperm = 200} may be small, depending on the amount of computing time available.

  If \code{argvals} is not specified and \code{yfdPar} is a \code{fd} object,
  it defaults to 101 equally-spaced points on the range of \code{yfdPar}.
}
\value{
  A list with components
  \itemize{
    \item{pval}{the observed p-value of the permutation test.}
    \item{qval}{the \code{q}th quantile of the null distribution.}
    \item{Fobs}{the observed maximal F-statistic.}
    \item{Fnull}{
      a vector of length \code{nperm} giving the observed values of the
      permutation distribution.
    }
    \item{Fvals}{the pointwise values of the observed F-statistic.}
    \item{Fnullvals}{
      the pointwise values of of the permutation observations.
    }
    \item{pvals.pts}{pointwise p-values of the F-statistic.}
    \item{qvals.pts}{
      pointwise \code{q}th quantiles of the null distribution
    }
    \item{fRegressList}{
      the result of \code{fRegress} on the observed data
    }
    \item{argvals}{
      argument values for evaluating the F-statistic if \code{yfdPar} is
      a functional data object.
    }
  }
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\section{Side Effects}{
  a plot of the functional observations
} \seealso{
  \code{\link{fRegress}}
  \code{\link{Fstat.fd}}
%  \code{\link{tstat.fd}}
}
%\examples{}
\keyword{smooth}
