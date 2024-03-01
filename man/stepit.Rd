\name{stepit}
\alias{stepit}
\title{
Compute a new step size for a linear search within an optimization algorithm.}
\description{
Using cubic interpolation, a new step size is computed for minimizing a 
function value.}
\usage{
  stepit(linemat, ips, dblwrd, MAXSTEP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linemat}{
  A 3 by 5 matrix containing step sizes, function values and slopes at
  multiple values.}
  \item{ips}{If 1, previous slope was positive.}
  \item{dblwrd}{Vector of length 2:  dblwrd[1] TRUE means step halved
  dblwrd[2] TRUE means step doubled.}
  \item{MAXSTEP}{maximum allowed size of a step.}
}
\value{A positive step size.}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
