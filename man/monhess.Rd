\name{monhess}
\alias{monhess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate the Hessian matrix of a monotone function
}
\description{
Evaluates the hessian or second derivative matrix of a monotone function 
with respect to the coefficients definingthe log-first derivative $W(t)$ 
at each of a set of argument values.
}
\usage{
monhess(x, Wfd, basislist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numerical vector at which function and derivative are
                 evaluated.}
  \item{Wfd}{A functional data object.}
  \item{basislist}{A list containing values of basis functions.}
}
\value{
A three dimensional array with first dimension corresponding to argument values 
and second and third dimensions to number of basis functions defining $W$.
}
\references{
  Ramsay, James O., Hooker, G. and Graves, S. (2009), \emph{Functional
    Data Analysis with R and Matlab}, Springer, New York.
    
  Ramsay, James O., and Silverman, 
  Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, 
  Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\author{
  J. O. Ramsay
}
\seealso{
  \code{\link{monfn}},
  \code{\link{mongrad}},
  \code{\link{landmarkreg}},
  \code{\link{smooth.morph}}
  }
  
