\name{fRegress.CV}
\alias{fRegress.CV}
\title{
  Computes Cross-validated Error Sum of Squares for a
Functional Regression Model
}
\description{
For a functional regression model with a scalar dependent variable,
a cross-validated error sum of squares is computed.  This function
aids the choice of smoothing parameters in this model using the
cross-validated error sum of squares criterion.
}
\usage{
fRegress.CV(yvec, xfdlist, betalist)
}
\arguments{
\item{yvec}{
a vector of dependent variable values.
}
\item{xfdlist}{
a list whose members are functional parameter objects
specifying functional independent variables.  Some
of these may also be vectors specifying scalar independent
variables.
}
\item{betalist}{
a list containing functional parameter objects specifying the
regression functions and their level of smoothing.
}
}
\value{
the sum of squared errors in predicting \code{yvec}.
}
\seealso{
\code{\link{fRegress}}, 
\code{\link{fRegress.stderr}}
}
\examples{
#See the analyses of the daily weather data.
}
% docclass is function
\keyword{smooth}
