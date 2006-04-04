\name{linmod}
\alias{linmod}
\title{
  Fit Fully Functional Linear Model
}
\description{
A functional dependent variable is approximated by a single
functional covariate, and the
covariate can affect the dependent variable for all
values of its argument.  The regression function is a bivariate function.
}
\usage{
linmod(xfdobj, yfdobj, wtvec=rep(1,nrep),
       xLfdobj=int2Lfd(2), yLfdobj=int2Lfd(2),
       xlambda=0, ylambda=0)
}
\arguments{
\item{xfdobj}{
a functional data object for the covariate
}
\item{yfdobj}{
a functional data object for the dependent variable
}
\item{wtvec}{
a vector of weights for each observation.
}
\item{xLfdobj}{
either a nonnegative integer or a linear differential operator
object.  This operator is applied to the regression function's
first argument.
}
\item{yLfdobj}{
either a nonnegative integer or a linear differential operator
object.  This operator is applied to the regression function's
second argument.
}
\item{xlambda}{
a smoothing parameter for the first argument of the regression
function.
}
\item{ylambda}{
a smoothing parameter for the second argument of the regression
function.
}
}
\value{
a named list of length 3 with the following entries:

\item{alphafd}{
the intercept functional data object.
}
\item{regfd}{
a bivariate functional data object for the regression function.
}
\item{yhatfd}{
a functional data object for the approximation to the dependent variable
defined by the linear model, if the dependent variable is functional.
Otherwise the matrix of approximate values.
}
}
\seealso{
fRegress
}
\examples{
#See the prediction of precipitation using temperature as
#the independent variable in the analysis of the daily weather
#data.
}
% docclass is function
\keyword{smooth}
