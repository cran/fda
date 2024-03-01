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
This function is not normally used in a functional data analysis to
smooth data, since function \code{smooth.basis} is provided for that job.
}
\usage{
project.basis(y, argvals, basisobj, penalize=FALSE)
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
}
\value{
the matrix of coefficients defining the least squares approximation.
This matrix has as many rows are there are basis functions, as many
columns as there are curves, and if the data are multivariate, as many
layers as there are functions.
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
  \code{\link{smooth.basis}}
}
% docclass is function
\keyword{smooth}
