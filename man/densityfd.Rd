\name{densityfd}
\alias{densityfd}
\title{Compute a Probability Density Function }
\description{
Like the regular S-PLUS function density, this function computes a probability density function for a sample of values of a random variable. However, in this case the density function is defined by a functional data object Wfdobj along with a normalizing constant C. The density function p(x) has the form 
p(x) = C exp[W(x)] where function W(x) is defined by the functional data object Wfdobj. 
}
\usage{
densityfd(x, Wfdobj, Lfdobj=3, lambda=0, conv=0.0001, iterlim=20,
          active=2:nbasis, dbglev=1)
}
\arguments{
  \item{x}{A vector of variable values for which the density is required. }
  \item{Wfdobj }{A functional data object defining a single univariate function. }
  \item{Lfdobj }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is evaluated rather than the functions themselves. }
  \item{lambda }{A nonnegative value controlling the amount of roughness in the data. }
  \item{conv }{A criterion for convergence of the iterations. }
  \item{iterlim }{A limit on the number of iterations. }
  \item{active }{A logical vector of length equal to the number of coefficients defining Wfdobj. If an entry is T, the corresponding coefficient is estimated, and if F, it is held at the value defining the argument Wfdobj. Normally the first coefficient is set to 0 and not estimated, since it is assumed that W(0) = 0. }
  \item{dbglev}{Either 0, 1, or 2. This controls the amount information printed out on each iteration, with 0 implying no output, 1 intermediate output level, and 2 full output. }
}
\details{
The goal of the function is provide a smooth density function estimate that approaches some target density by an amount that is controlled by the linear differential operator Lfdobj and the penalty parameter. For example, if the second derivative of W(t) is penalized heavily, this will force the function to approach a straight line, which in turn will force the density function itself to be nearly normal or Gaussian. Similarly, to each textbook density function there corresponds a W(t), and to each of these in turn their corresponds a linear differential operator that will, when apply to W(t), produce zero as a result. 

To plot the density function or to evaluate it, evaluate Wfdobj, exponentiate the resulting vector, and then divide by the normalizing constant C. 
}
\value{
A list containing: 
Wfdobj 
A functional data object defining function W(x) that that optimizes the fit to the data of the monotone function that it defines. 
C 
The normalizing constant. 
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
\references{
See Chapter 5 in J. O. Ramsay and B. W. Silverman (2002) Applied Functional Data Analysis for an example of density estimation and further references. 
}
\note{}
\author{}
\seealso{
posfd, smooth.monotone, eval.monfd 
}
\examples{}
\keyword{}
