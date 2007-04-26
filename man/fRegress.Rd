\name{fRegress}
\alias{fRegress}
\title{
  A Functional Regression Analysis of the Concurrent Type
}
\description{
This function carries out a functional regression analysis,
where either the dependent variable or one or more
independent variables are functional.  Non-functional variables
may be included on either side of the equation.  In a concurrent
functional linear model all function variables are all evaluated at
a common time or argument value $t$.   That is, the fit is
defined in terms of the behavior of all variables at a fixed time,
or in terms of "now" behavior.
}
\usage{
fRegress(yfdPar, xfdlist, betalist, wt=rep(1,N))
}
\arguments{
\item{yfdPar}{
the dependent variable object.  It may be an object of
three possible classes:
\itemize{
   \item a vector if the dependent variable is scalar.
   \item a functional data object if the dependent variable is functional.
   \item a functional parameter object if the dependent variable is functional,    and if it is necessary to smooth the prediction of the dependent variable.
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
  weights for weighted least squares 
}
}
\value{
a named list of length 6 with these members:

\item{yfdPar}{
yhe first argument in the call to \code{fRegress}.
}
\item{xfdlist}{
the second argument in the call to \code{fRegress}.
}
\item{betalist}{
the third argument in the call to \code{fRegress}.
}
\item{betaestlist}{
a list of length equal to the number of independent
variables and with members having the same functional parameter
structure as the corresponding
members of \code{betalist}.  These are the estimated regression
coefficient functions.
}
\item{yhatfdobj}{
a functional data object if the dependent variable
is functional or a vector if the dependent variable is scalar.  This
is the set of predicted by the functional regression model for the
dependent variable.
}
\item{Cmatinv}{
a matrix containing the inverse of the coefficient
matrix for the linear equations that define the solution to the
regression problem.  This matrix is required for function 
\code{\link{fRegress.stderr}} that estimates confidence regions for the
regression coefficient function estimates.
}
}
\details{

In the computation of regression function estimates, all
independent variables are treated as if they are functional.  If
argument \code{xfdlist} contains one or more vectors, these
are converted to functional data objects having the constant basis
with coefficients equal to the elements of the vector.

 Needless to say, if all the variables in the model are scalar,
use this function, but rather either \code{ls} or \code{lsfit}.
}
\seealso{
\code{\link{fRegress.stderr}}, 
\code{\link{fRegress.CV}}, 
\code{\link{linmod}}
}
\examples{
#See the Canadian daily weather data analyses in the file
# \code{daily.R} for 
#examples of all the cases covered by \code{fRegress}.
}
% docclass is function
\keyword{smooth}
