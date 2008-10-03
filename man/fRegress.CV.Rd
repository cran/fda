\name{fRegress.CV}
\alias{fRegress.CV}
\title{
  Computes Cross-validated Error Sum of Integrated Squared Errors for a
  Functional Regression Model
}
\description{
  For a functional regression model, a cross-validated error sum of
  squares is computed.  For a functional dependent variable this is the
  sum of integrated squared errors. For a scalar response, this function
  has been superceded by the OCV and gcv elements returned by
  \code{fRegress}. This function aids the choice of smoothing parameters
  in this model using the cross-validated error sum of squares
  criterion.
}
\usage{
fRegress.CV(y, xfdlist, betalist, CVobs=NULL, ...)
}
\arguments{
  \item{y}{
    the dependent variable object.  It may be an object of three
    possible classes:

    \itemize{
      \item{scalar}{ a vector if the dependent variable is scalar.}
      \item{fd}{
	a functional data object if the dependent variable is
	functional.
      }
      \item{fdPar}{
	a functional parameter object if the dependent variable is
	functional
      }
    }
  }
  \item{xfdlist}{
    a list whose members are functional parameter objects specifying
    functional independent variables.  Some of these may also be vectors
    specifying scalar independent variables.
  }
  \item{betalist}{
    a list containing functional parameter objects specifying the
    regression functions and their level of smoothing.
  }
  \item{CVobs}{
    a vector giving the indices of observations to be deleted one at a
    time in the cross validation. Defaults to all of them.
    Observations not in \code{CVobs} will be included in all
    bootstrap / cross validation samples.
  }
  \item{\dots}{
    optional arguments not used by \code{fRegress.CV} but needed for
    superficial compatibability with \code{fRegress} methods.
  }
}
\value{
  A list containing
  \item{SSE.CV}{ The sum of squared errors, or integrated squared errors}
  \item{errfd.cv}{
    Either a vector or a functional data object giving the
    cross-validated errors.
  }
}
\seealso{
\code{\link{fRegress}},
\code{\link{fRegress.stderr}}
}
\examples{

#### A scalar response example.

annualprec <- log10(apply(CanadianWeather$dailyAv[,,"Precipitation.mm"],
                          2,sum))
# The simplest 'fRegress' call is singular with more bases
# than observations, so we use a small basis for this example
smallbasis  <- create.fourier.basis(c(0, 365), 25)
# There are other ways to handle this,
# but we will not discuss them here
tempfd <- smooth.basis(day.5,
          CanadianWeather$dailyAv[,,"Temperature.C"], smallbasis)$fd

# Set up a list of predictors

xfdlist <- list(const=rep(1, 35), tempfd=tempfd)

# And some functional parameter objects to define the coefficient functions.
# The intercept must be constant for a scalar response
betabasis1 <- create.constant.basis(c(0, 365))
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)

betafd2     <- with(tempfd, fd(basisobj=basis, fdnames=fdnames))
# convert to an fdPar object
betafdPar2  <- fdPar(betafd2)

betalist <- list(const=betafdPar1, tempfd=betafdPar2)

# Now call fRegress.

precip.Temp <- fRegress(annualprec, xfdlist, betalist)

precip.cv = fRegress.CV(annualprec,xfdlist,betalist)


## Are the temperature profiles different accross regions, a functional response
# model

region.contrasts <- model.matrix(~factor(CanadianWeather$region))
rgnContr3 <- region.contrasts
dim(rgnContr3) <- c(1, 35, 4)
dimnames(rgnContr3) <- list('', CanadianWeather$place, c('const',
   paste('region', c('Atlantic', 'Continental', 'Pacific'), sep='.')) )

const365 <- create.constant.basis(c(0, 365))
region.fd.Atlantic <- fd(matrix(rgnContr3[,,2], 1), const365)
%str(region.fd.Atlantic)
region.fd.Continental <- fd(matrix(rgnContr3[,,3], 1), const365)
region.fd.Pacific <- fd(matrix(rgnContr3[,,4], 1), const365)
region.fdlist <- list(const=rep(1, 35),
     region.Atlantic=region.fd.Atlantic,
     region.Continental=region.fd.Continental,
     region.Pacific=region.fd.Pacific)
%str(TempRgn.mdl$betalist)

beta1 <- with(tempfd, fd(basisobj=basis, fdnames=fdnames))
beta0 <- fdPar(beta1)
betalist <- list(const=beta0, region.Atlantic=beta0,
             region.Continental=beta0, region.Pacific=beta0)

TempRgn <- fRegress(tempfd, region.fdlist, betalist)

temprgn.ocv <- fRegress.CV(tempfd, region.fdlist, betalist)
}
% docclass is function
\keyword{smooth}
