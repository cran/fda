\name{project.basis}
\alias{project.basis}
\title{
  Approximate Functional Data Using a Basis
}
\description{
A vector or matrix of discrete data is projected into the space
spanned by the values of a set of basis functions.  This amounts to
a least squares regression of the data on to the values of the basis
functions.  A small penalty can be applied to deal with situations in
which the number of basis functions exceeds the number of basis points.
This function is used with function \code{data2fd}, and is not
normally used directly in a functional data analysis.
}
\usage{
project.basis(y, argvals, basisobj, penalize=FALSE, returnMatrix=FALSE)
}
\arguments{
\item{y}{
a vector or matrix of discrete data.
}
\item{argvals}{
a vector containing the argument values correspond to the
values in \code{y}.
}
\item{basisobj}{
a basis object.
}
\item{penalize}{
a logical variable.  If TRUE, a small roughness penalty is applied
to ensure that the linear equations defining the least squares
solution are linearly independent or nonsingular.
}
\item{returnMatrix}{
logical:  If TRUE,  a two-dimensional is returned using a
special class from the Matrix package.
}
}
\value{
the matrix of coefficients defining the least squares approximation.
This matrix has as many rows are there are basis functions, as many
columns as there are curves, and if the data are multivariate, as many
layers as there are functions.
}
\seealso{
\code{\link{Data2fd}}
}
% docclass is function
\keyword{smooth}
