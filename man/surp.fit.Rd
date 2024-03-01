
\name{surp.fit}
\alias{surp.fit}
\title{
  Evaluate the fit of surprisal curves to binned psychometric data.
}
\description{
  Evaluate the error sum of squares, its gradient and its hessian for the fit of 
  surprisal curves to binned psychometric data.  The function value is optimized
  by function \code{smooth.surp} in package \code{TestGardener.}
}
\usage{
surp.fit(x, surpList)
}
\arguments{
  \item{x}{The parameter vector, which is the vectorized form of the K by M-1   coefficient matrix for the functional data object.}
  \item{surpList}{
  A named list object containing objects essential to evaluating the fitting 
  criterion.  See \code{smooth.surp.R} for the composition of this list.
  }
}
\value{
A named list object for the returned objects with these names:
    \item{PENSSE:}{The error sum of squares associated with parameter value \code{x}.}
    \item{DPENSSE:}{A column vector containing gradient of the error sum of squares.}
    \item{D2PENSSE:}{A square matrix of hessian values.}
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{smooth.surp}}
}
