\name{posfd}
\alias{posfd}
\title{Smooth Data with a Positive Function using a Roughness Penalty }
\description{
This function smooths data from a single curve or function using a roughness penalty. Unlike function smooth. basis, however, the function is constrained to be positive. The function is defined by a functional data object Wfdobj, which defines the logarithm of the positive function. The roughness of Wfdobj, and hence of the positive smoothing function, is definable in a wide variety of ways using either derivatives or a linear differential operator. 
}
\usage{
posfd(y, argvals, Wfdobj, Lfd=3, lambda=0, conv=0.0001, iterlim=20,
      dbglev=1)
}
\arguments{
  \item{y }{A vector of observations on a single curve. }
  \item{argvals }{A vector of argument values corresponding to y. }
  \item{Wfdobj }{A functional data object defining a single univariate function. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
  \item{lambda }{A nonnegative value controlling the amount of roughness in the data. }
  \item{conv }{A criterion for convergence of the iterations. }
  \item{iterlim }{A limit on the number of iterations. }
  \item{dbglev }{Either 0, 1, or 2. This controls the amount information printed out on each iteration, with 0 implying no output, 1 intermediate output level, and 2 full output. }
}
\details{
The computation of the positive function requires iterative optimization of the fitting criterion, and is therefore more computationally intensive than unconstrained smoothing. For this reason, only one curve at a time is allowed to be smoothed. 

The computational problem of estimating a positive function is very similar to that for estimating a density function, except that there is no normalizing constraint. For this reason, this function has much in common with function densityfd. 

To plot the positive function or to evaluate it, evaluate Wfdobj, and then exponentiate the resulting vector. 
}
\value{
A list containing: 
Wfdobj 
A functional data object defining function W(x) that that optimizes the fit to the data of the monotone function that it defines. 
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
