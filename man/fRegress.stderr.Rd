\name{fRegress.stderr}
\alias{fRegress.stderr}
\title{
  Compute Standard errors of Coefficient Functions
Estimated by Functional Regression
Analysis
}
\description{
Function \code{fRegress} carries out a functional regression analysis
of the concurrent kind, and estimates a regression coefficient function
corresponding to each independent variable, whether it is scalar or
functional.  This function uses the list that is output by \code{fRegress}
to provide standard error functions for each regression function.  These
standard error functions are pointwise, meaning that sampling standard
deviation functions only are computed, and not sampling covariances.
}
\usage{
fRegress.stderr(fRegressList, y2cMap, SigmaE)
}
\arguments{
\item{fRegressList}{
the named list of length six that is returned from a call to function
\code{fRegress}.
}
\item{y2cMap}{
a matrix that contains the linear transformation that
takes the raw data values into the coefficients defining a smooth
functional data object. Typically, this matrix is returned from a
call to function \code{smooth.basis} that generates the
dependent variable objects.  If the dependent variable is scalar,
this matrix is an identity matrix of order equal to the length
of the vector.
}
\item{SigmaE}{
either a matrix or a bivariate functional data object
according to whether the dependent variable is scalar or functional,
respectively.
This object has a number of replications equal to
the length of the dependent variable object.  It contains an estimate
of the variance-covariance matrix or function for the residuals.
}
}
\value{
a named list of length 3 containing:

\item{betastderrlist}{
a list object of length
the number of independent variables. Each member contains a
functional parameter object
for the standard error of a regression function.
}
\item{bvar}{
a symmetric matrix containing sampling variances and
covariances for the matrix of regression coefficients
for the regression functions.  These are stored
column-wise in defining BVARIANCE.
}
\item{c2bMap}{
a matrix containing the mapping from response variable
coefficients to coefficients for regression coefficients.
}
}
\seealso{
\code{\link{fRegress}}, 
\code{\link{fRegress.CV}}
}
\examples{
#See the weather data analyses in the file \code{daily.ssc} for
#examples of the use of function \code{fRegress.stderr}.
}
% docclass is function
\keyword{smooth}