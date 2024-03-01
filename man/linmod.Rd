\name{linmod}
\alias{linmod}
\title{
  Fit Fully Functional Linear Model
}
\description{
A functional dependent variable \eqn{y_i(t)} is approximated by a single
functional covariate \eqn{x_i(s)} plus an intercept function \eqn{\alpha(t)},
and the covariate can affect the dependent variable for all
values of its argument. The equation for the model is

  \deqn{y_i(t) = \beta_0(t) + \int \beta_1(s,t) x_i(s) ds + e_i(t)}

for \eqn{i = 1,...,N}. The regression function \eqn{\beta_1(s,t)} is a
bivariate function. The final term \eqn{e_i(t)} is a residual, lack of
fit or error term.   There is no need for values \eqn{s} and \eqn{t} to
be on the same continuum.

}
\usage{
linmod(xfdobj, yfdobj, betaList, wtvec=NULL)
}
\arguments{
\item{xfdobj}{
a functional data object for the covariate
}
\item{yfdobj}{
a functional data object for the dependent variable
}
\item{betaList}{
a list object of length 2.  The first element is a functional parameter
object specifying a basis and a roughness penalty for the intercept term.
The second element is a bivariate functional parameter object for the
bivariate regression function.}
\item{wtvec}{
a vector of weights for each observation.  Its default value is NULL,
in which case the weights are assumed to be 1.
}
}
\value{
a named list of length 3 with the following entries:

\item{beta0estfd}{
the intercept functional data object.
}
\item{beta1estbifd}{
a bivariate functional data object for the regression function.
}
\item{yhatfdobj}{
a functional data object for the approximation to the dependent variable
defined by the linear model, if the dependent variable is functional.
Otherwise the matrix of approximate values.
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
\seealso{
  \code{\link{bifdPar}},
  \code{\link{fRegress}}
}
\source{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009)
    \emph{Functional Data Analysis in R and Matlab}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York
}
\examples{
#See the prediction of precipitation using temperature as
#the independent variable in the analysis of the daily weather
#data, and the analysis of the Swedish mortality data.
}
% docclass is function
\keyword{smooth}
