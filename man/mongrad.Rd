\name{mongrad}
\alias{mongrad}
\title{
Evaluate the gradient of a monotone function
}
\description{
Evaluates the gradient of a monotone function with respect to the coefficients defining
the log-first derivative $W(t)$ at each of a set of argument values.}
\usage{
mongrad(x, Wfdobj, basislist=vector("list",JMAX), 
                    returnMatrix=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numerical vector at which function and derivative are
                 evaluated.}
  \item{Wfdobj}{A functional data object.}
  \item{basislist}{A list containing values of basis functions.}
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.}
}
\value{
A matrix with as many rows as argument values and as many columns as basis functions
defining $W$.
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
  \code{\link{monhess}},
  \code{\link{landmarkreg}},
  \code{\link{smooth.morph}}
}
