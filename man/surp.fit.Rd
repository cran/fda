
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
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
\emph{Journal of Educational and Behavioral Statistics}, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  \emph{Psych}, 2, 347-360.

http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{smooth.surp}}
}
