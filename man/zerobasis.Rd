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
\seealso{
  \code{\link{smooth.surp}},
  \code{\link{eval.surp}}
}
