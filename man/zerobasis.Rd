\name{zerobasis}
\alias{zerobasis}
\title{
Orthonormal Matrix with Columns Summing to Zero}
\description{
The matrix is k by k-1, columns are orthonormal and sum to zero. This matrix
is used when quantities are unchanged to adding a constant, such as in
multinomial data modelling.}
\usage{
zerobasis(k)
}
\arguments{
  \item{k}{A positive integer greater than one.}
}
\details{
There are many ways of generating these matrix. This functions uses
the discrete Fourier series less the constant vector.}
\value{
  A \code{k} by \code{k-1} matrix with orthonormal columns that sum to zero.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{smooth.surp}},
  \code{\link{eval.surp}}
}
