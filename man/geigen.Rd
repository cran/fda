\name{geigen}
\alias{geigen}
\title{
  Generalized eigenanalysis
}
\description{
  Find matrices L and M to maximize

  tr(L'AM) / sqrt(tr(L'BL) tr(M'CM'))

  where A = a p x q matrix, B = p x p symmetric, positive definite
  matrix, B = q x q symmetric positive definite matrix, L = p x s
  matrix, and M = q x s matrix, where s = the number of non-zero
  generalized eigenvalues of A.
}
\usage{
geigen(Amat, Bmat, Cmat)
}
\arguments{
  \item{Amat}{ a numeric matrix }
  \item{Bmat}{
    a symmetric, positive definite matrix with dimension = number of
    rows of A
  }
  \item{Cmat}{
    a symmetric, positive definite matrix with dimension = number of
    columns of A
  }
}
\value{
  list(values, Lmat, Mmat)
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
  \code{\link{eigen}}
}
\examples{
A <- matrix(1:6, 2)
B <- matrix(c(2, 1, 1, 2), 2)
C <- diag(1:3)
ABC <- geigen(A, B, C)
}
% docclass is function
\keyword{array}
