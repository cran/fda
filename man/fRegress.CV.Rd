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
#fRegress.CV(y, xfdlist, betalist, wt=NULL, CVobs=1:N,
#            returnMatrix=FALSE, ...)

#NOTE:  The following is required by CRAN rules that
# function names like "as.numeric" must follow the documentation
# standards for S3 generics, even when they are not.
# Please ignore the following line:
\method{fRegress}{CV}(y, xfdlist, betalist, wt=NULL, CVobs=1:N,
            returnMatrix=FALSE, ...)
}
\arguments{
  \item{y}{
    the dependent variable object.
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
  \item{wt}{
    weights for weighted least squares.  Defaults to all 1's.
  }
  \item{CVobs}{
    Indices of observations to be deleted.  Defaults to 1:N.
  }
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
  \item{\dots}{
    optional arguments not used by \code{fRegress.CV} but needed for
    superficial compatibability with \code{fRegress} methods.
  }
}
\value{
  A list containing
  \item{SSE.CV}{ The sum of squared errors, or integrated squared errors}
  \item{errfd.cv}{ Either a vector or a functional data object giving the
    cross-validated errors }
}
\seealso{
\code{\link{fRegress}},
\code{\link{fRegress.stderr}}
}
\examples{
#See the analyses of the Canadian daily weather data.
}
% docclass is function
\keyword{smooth}
