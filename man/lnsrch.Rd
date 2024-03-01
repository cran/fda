\name{lnsrch}
\alias{lnsrch}
\title{
Search along a line for a minimum within an optimisation algorithm.
}
\description{
This is a version of the function in Numerical Recipes.  
It is initialized with a function value and gradient, and it does a 
series of quadratic searches until a convergence criterion is reached.
This version includes code for display the progress of iteration for
debugging purposes.
}
\usage{
lnsrch(xold, fold, g, p, func, dataList, stpmax, itermax=20, TOLX=1e-10, dbglev=0)
}
\arguments{
  \item{xold}{The current parameter vector value.}
  \item{fold}{The current function value.}
  \item{g}{The current gradient vector.}
  \item{p}{The current search direction vector.}
  \item{func}{The name of the function being optimized.}
  \item{dataList}{A list object containing objects specifying the 
  function to be minimized.}
  \item{stpmax}{The maximum step size.}
  \item{itermax}{The maximum number of iterations. Default 20}
  \item{TOLX}{Convergence criterion.}
  \item{dbglev}{The level of output.  0: no output; 1: function and slope
  at each iteration.  2 and above: also results within iteration.}
}
\value{
A named list containing:
    \item{x:}{The parameter vector at the final minimum.}
    \item{check:}{A boolean value indicating that problems were encountered.}
}
\author{Juan Li and James Ramsay}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.

  Press, W. H., Taukolsky, S. A., Vetterline, W. T. and Flannery, B. P. (2020) 
  \emph{Numerical Recipes, Third Edition}, Cambridge.
}
\seealso{
  \code{\link{smooth.surp}}
}
