\name{smooth.monotone}
\alias{smooth.monotone}
\title{Monotone Smoothing of Data }
\description{
When the discrete data that are observed reflect a smooth strictly increasing or strictly decreasing function, it is often desirable to smooth the data with a strictly monotone function, even though the data themselves may not be monotone due to observational error. An example is when data are collected on the size of a growing organism over time. This function computes such a smoothing function, but, unlike other smoothing functions, for only for one curve at a time. 

The smoothing function minimizes a weighted error sum of squares criterion. This minimization requires iteration, and therefore is more computationally intensive than normal smoothing. 
}
\usage{
smooth.monotone(x, y, wt=rep(1, nobs), Wfdobj, zmat=matrix(1, nobs, 1),
Lfdobj=1, lambda=0, conv=0.0001, iterlim=20, active=c(FALSE, rep(TRUE, ncvec -1)), dbglev=1)
}
\arguments{
  \item{x }{A vector of argument values associated with the data. }
  \item{y }{A vector of observations to be smoothed. }
  \item{Wfdobj }{A functional data object that determines the final monotone smooth. It must be univariate and define a single functional observation. }
  \item{wt }{A vector of positive weights for the observations. The default is all weights one. }
  \item{zmat }{A design matrix that permits an additive constant to possibly depend on covariate variables. The default is a column of one's. }
  \item{Lfdobj }{Either a nonnegative integer or a linear differential operator object that defines a rougness penalty to be applied to functional data object Wfdobj. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
  \item{lambda }{A nonnegative smoothing parameter that controls the degree of roughness in Wfdobj, and consequently in the monotone smoothing function. }
  \item{conv }{A criterion for convergence of the iterations. }
  \item{iterlim }{A limit on the number of iterations. }
  \item{active }{A logical vector of length equal to the number of coefficients defining Wfdobj. If an entry is T, the corresponding coefficient is estimated, and if F, it is held at the value defining the argument Wfdobj. Normally the first coefficient is set to 0 and not estimated, since it is assumed that W(0) = 0. }
  \item{dbglev }{Either 0, 1, or 2. This controls the amount information printed out on each iteration, with 0 implying no output, 1 intermediate output level, and 2 full output. }
}
\details{
The monotone smoothing function has the following structure: 

f(x) = b0 + b1 D-1 [exp W(x)] 
Here the operator D-1 means taking the indefinite integral of the function that follows, namely exp W(x). The smoothing function f(x) is, by this structure, determined by three objects that need to be estimated from the data: 

W(x)
A functional data object that is first exponentiated and then the result integrated. This is the heart of the monotone smooth. The closer W(x) is to zero, the closer the monotone smooth becomes a straight line. The closer W(x) becomes a constant, the more the monotone smoother becomes an exponential function. It is assumed that W(0) = 0. 
b0
An intercept term that determines the value of the smoothing function at x = 0. 
b1
A regression coefficient that determines the slope of the smoothing function at x = 0. 
In addition, it is possible to have the intercept b0 depend in turn on the values of one or more covariates through the design matrix Zmatas follows: 

b0 = Zc 
In this case, the single intercept coefficient is replaced by the regression coefficients in vector c multipying the design matrix. 

}
\value{
A list containing: 
Wfdobj 
A functional data object defining function W(x) that that optimizes the fit to the data of the monotone function that it defines. 
beta 
The optimal regression coefficient values. 
Flist 
A list containing results for the final converged solution: 
f 
The optimal function value being minimized. 
grad 
The gradient vector at the optimal solution. 
norm 
The norm of the gradient vector at the optimal solution. 
iternum 
The number of iterations. 
iterhist 
A iternum+1 by 5 matrix containing the iteration history. 
}
\references{}
\note{}
\author{}
\seealso{}
\examples{
}
\keyword{}
