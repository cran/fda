\name{linmod.fd}
\alias{linmod.fd}
\title{Fit One of Three Types of Functional Linear Model. }
\description{
The linear model function described below fits the three types of linear models described in Chapters 9, 10, and 11. At this point the function can only handle a single functional independent variable. The regularization features described in Chapters 10 and 11 are incorporated. 
To fit a linear model. Three cases are considered: 

functional dependent variable and multivariate independent variable, and the functional dependent variable can be multivariate or vector-valued. 
multivariate dependent variable and functional independent variable, 
functional dependent variable and a single functional independent variable. 
}
\usage{
linmod.fd(xfd, yfd, wtvec=rep(1, nrep), xLfd=2, yLfd=2, xlambda=0, ylambda=0, zmatrnk=p)
}
\arguments{
  \item{xfd }{If the independent variable is multivariate, a design matrix. If the independent variable is functional, a fd object. }
  \item{yfd }{If the dependent variable is multivariate, a design matrix. If the dependent variable is functional, a fd object. }
  \item{wtvec }{A vector of weights for the replications. By default these are 1's. }
  \item{xLfd }{For the independent variable, the order derivative to be penalized if an integer, or a linear differential operator if a functional data object. }
  \item{yLfd }{For the dependent variable, the order derivative to be penalized if an integer, or a linear differential operator if a functional data object. }
  \item{xlambda }{A smoothing parameter for the independent variable. This is zero by default. }
  \item{ylambda }{A smoothing parameter for the dependent variable. This is zero by default. }
  \item{zmatrnk }{The actual rank of independent variable matrix for the function DV/multivariate IV case. This is permitted to be less than the number of columns of this matrix. }
}
\details{
}
\value{
Returns a list with the following entries: 
"alpha" The intercept, either a scalar or a functional data object as the model requires. If the dependent variable is multivariate, there is one for each column of the matrix input as argument yfd. Otherwise this is a single constant. 
"regfd" A functional data object for the regression function. The structure of this object depends on the three cases: 
functional DV/multivariate IV: a univariate functional data object with a replication corresponding to each column dimension in the matrix input for argument xfd, and a function (3rd dimension in the coefficient matrix) for each function in the DV, 
multivariate DV/functional IV: a univariate functional data object having a replication correponding to each column of the matrix input as argument yfd. 
functional DV/functional IV: a bivariate functional object of the bifd class. 
"hatfd" A functional data object for the approximation to the dependent variable defined by the linear model, if the dependent variable is functional. Otherwise the matrix of approximate values. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{}
\keyword{}
