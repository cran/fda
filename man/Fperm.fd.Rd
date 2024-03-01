\name{Fperm.fd}
\alias{Fperm.fd}
\title{
  Permutation F-test for functional linear regression.
}
\description{
  \code{Fperm.fd} creates a null distribution for a test of no effect in
  functional linear regression.  It makes generic use of \code{fRegress}
  and permutes the \code{yfdPar} input.
}
\usage{
Fperm.fd(yfdPar, xfdlist, betalist, wt=NULL, nperm=200,
         argvals=NULL, q=0.05, plotres=TRUE, ...)
}
\arguments{
  \item{yfdPar}{
    the dependent variable object.  It may be an object of three
    possible classes:

    \describe{
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
    a list of length equal to the number of independent
    variables. Members of this list are the independent variables.  They
    be objects of either of these two classes:
    \describe{
      \item{vector:}{a vector if the independent dependent variable is scalar.}
      \item{fd:}{a functional data object if the dependent variable is 
        functional.}
    }
    In either case, the object must have the same number of replications
    as the dependent variable object.  That is, if it is a scalar, it
    must be of the same length as the dependent variable, and if it is
    functional, it must have the same number of replications as the
    dependent variable.
  }
  \item{betalist}{
    a list of length equal to the number of independent variables.
    Members of this list define the regression functions to be
    estimated.  They are functional parameter objects.  Note that even
    if corresponding independent variable is scalar, its regression
    coefficient will be functional if the dependent variable is
    functional.  Each of these functional parameter objects defines a
    single functional data object, that is, with only one replication.
  }
  \item{wt}{
    weights for weighted least squares, defaults to all 1.
  }
  \item{nperm}{
    number of permutations to use in creating the null distribution.
  }
  \item{argvals}{
    If \code{yfdPar} is a \code{fd} object, the points at which to
    evaluate the point-wise F-statistic.
  }
  \item{q}{
    Critical upper-tail quantile of the null distribution to compare to
    the observed F-statistic.
  }
  \item{plotres}{
    Argument to plot a visual display of the null distribution
    displaying the \code{q}th quantile and observed F-statistic.
  }
  \item{...}{
    Additional plotting arguments that can be used with \code{plot}.
  }
}
\details{
  An F-statistic is calculated as the ratio of residual variance to
  predicted variance. The observed F-statistic is returned along with
  the permutation distribution.

  If \code{yfdPar} is a \code{fd} object, the maximal value of the
  pointwise F-statistic is calculated. The pointwise F-statistics are
  also returned.

  The default of setting \code{q = 0.95} is, by now, fairly
  standard. The default \code{nperm = 200} may be small, depending on
  the amount of computing time available.

  If \code{argvals} is not specified and \code{yfdPar} is a \code{fd}
  object, it defaults to 101 equally-spaced points on the range of
  \code{yfdPar}.
}
\value{
  A list with the following components:

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
    argument values for evaluating the F-statistic if \code{yfdPar} is a
    functional data object.
  }
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\section{Side Effects}{
  a plot of the functional observations
} 
\seealso{
  \code{\link{fRegress}},
  \code{\link{Fstat.fd}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
##
## 1.  yfdPar = vector
##
annualprec <- log10(apply(
    CanadianWeather$dailyAv[,,"Precipitation.mm"], 2,sum))

#  set up a smaller basis using only 40 Fourier basis functions
#  to save some computation time

smallnbasis <- 40
smallbasis  <- create.fourier.basis(c(0, 365), smallnbasis)
tempfd      <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
                            smallbasis)$fd
constantfd <- fd(matrix(1,1,35), create.constant.basis(c(0, 365)))

xfdlist <- vector("list",2)
xfdlist[[1]] <- constantfd
xfdlist[[2]] <- tempfd[1:35]

betalist   <- vector("list",2)
#  set up the first regression function as a constant
betabasis1 <- create.constant.basis(c(0, 365))
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)
betalist[[1]] <- betafdPar1

nbetabasis  <- 35
betabasis2  <- create.fourier.basis(c(0, 365), nbetabasis)
betafd2     <- fd(matrix(0,nbetabasis,1), betabasis2)

lambda          <- 10^12.5
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
betafdPar2      <- fdPar(betafd2, harmaccelLfd365, lambda)
betalist[[2]]   <- betafdPar2

# Should use the default nperm = 200
# but use 10 to save test time for illustration

F.res2 = Fperm.fd(annualprec, xfdlist, betalist, nperm=100)

##
## 2.  yfdpar = Functional data object (class fd)
##
# The very simplest example is the equivalent of the permutation
# t-test on the growth data.

# First set up a basis system to hold the smooths

if (!CRAN()) {

Knots  <- growth$age
norder <- 6
nbasis <- length(Knots) + norder - 2
hgtbasis <- create.bspline.basis(range(Knots), nbasis, norder, Knots)

# Now smooth with a fourth-derivative penalty and a very small smoothing
# parameter

Lfdobj    <- 4
lambda    <- 1e-2
growfd    <- fd(matrix(0,nbasis,1),hgtbasis)
growfdPar <- fdPar(growfd, Lfdobj, lambda)

hgtfd <- smooth.basis(growth$age,
                      cbind(growth$hgtm,growth$hgtf),growfdPar)$fd

# Now set up factors for fRegress:

cbasis = create.constant.basis(range(Knots))

maleind = c(rep(1,ncol(growth$hgtm)),rep(0,ncol(growth$hgtf)))

constfd = fd( matrix(1,1,length(maleind)),cbasis)
maleindfd = fd( matrix(maleind,1,length(maleind)),cbasis)

xfdlist = list(constfd,maleindfd)

# The fdPar object for the coefficients and call Fperm.fd

betalist = list(fdPar(hgtfd,2,1e-6),fdPar(hgtfd,2,1e-6))

# Should use nperm = 200 or so,
# but use 10 to save test time

Fres = Fperm.fd(hgtfd,xfdlist,betalist,nperm=100)

par(oldpar)
}

}
\keyword{smooth}
