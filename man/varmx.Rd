\name{varmx}
\alias{varmx}
\title{
  Rotate a Matrix of Component Loadings using the VARIMAX Criterion
}
\description{
The matrix being rotated contains the values of the component
functional data objects computed in either a principal
components analysis or a canonical correlation analysis.
The values are computed over a fine mesh of argument values.
}
\usage{
varmx(amat, normalize=FALSE)
}
\arguments{
\item{amat}{
the matrix to be rotated.  The number of rows is
equal to the number of argument values \code{nx} used
in a fine mesh.  The number of columns is the number of
components to be rotated.
}
\item{normalize}{
either \code{TRUE} or \code{FALSE}.  If \code{TRUE}, the columns of
\code{amat} are normalized prior to computing the rotation
matrix.  However, this is seldom needed for functional data.
}
}
\value{
a square rotation matrix of order equal to the number
of components that are rotated.  A rotation matrix
$T$ has that property that $T'T = TT' = I$.
}
\details{
The VARIMAX criterion is the variance of the squared component values.
As this criterion is maximized with respect to a rotation of the
space spanned by the columns of the matrix, the squared loadings
tend more and more to be either near 0 or near 1, and this tends to
help with the process of labelling or interpreting the rotated matrix.
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
  \code{\link{varmx.pca.fd}},
  \code{\link{varmx.cca.fd}}
}
% docclass is function
\keyword{smooth}
