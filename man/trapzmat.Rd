\name{trapzmat}
\alias{trapzmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximate the functional inner product of two functional data objects using the 
trapezpoidal rule over a fine mesh of value.
}
\description{
The first two arguments are matrices of values of two functional data objects
of the same fine mesh of argument values.  The mesh is assumed to be equally spaced.}
\usage{
trapzmat(X,Y,delta=1,wt=rep(1,n))
}
\arguments{
  \item{X}{The first matrix of functional data object values over a fine mesh.}
  \item{Y}{The second matrix of functional data object values over the same fine mesh.}
  \item{delta}{The difference between adjacent argument values, assumed to 
  be a constant.}
  \item{wt}{An optional vector of weights for the products of pairs of argument 
  values, otherwise assumed to be all ones.}
}
\value{
A matrix of approximate inner products.  The number of rows the number of columns of 
\code{X} and the number of columns is the number of columns of \code{Y}.
}
\seealso{
  \code{\link{inprod}},
  \code{\link{inprod.bspline}}
}
