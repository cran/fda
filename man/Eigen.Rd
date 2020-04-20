\name{Eigen}
\alias{Eigen}
\title{
  Eigenanalysis preserving dimnames 
}
\description{
  Compute eigenvalues and vectors, assigning names to the eigenvalues
  and dimnames to the eigenvectors.  
}
\usage{
Eigen(x, symmetric, only.values = FALSE, valuenames)
}
\arguments{
  \item{x}{
    a square matrix whose spectral decomposition is to be computed.  
  }
  \item{symmetric}{
    logical:  If TRUE, the matrix is assumed to be symmetric (or 
    Hermitian if complex) and only its lower triangle (diagonal
    included) is used.  If 'symmetric' is not specified, the
    matrix is inspected for symmetry.
  }
  \item{only.values}{
    if 'TRUE', only the eigenvalues are computed and returned, otherwise
    both eigenvalues and eigenvectors are returned. 
  }
  \item{valuenames}{
    character vector of length nrow(x) or a character string that can be
    extended to that length by appening 1:nrow(x).

    The default depends on symmetric and whether
    \code{\link{rownames}} == \code{\link{colnames}}:  If
    \code{\link{rownames}} == \code{\link{colnames}} and
    symmetric = TRUE (either specified or determined by
    inspection), the default is "paste('ev', 1:nrow(x), sep='')".
    Otherwise, the default is colnames(x) unless this is NULL. 
  }
}
\details{
  1.  Check 'symmetric'  

  2.  ev <- eigen(x, symmetric, only.values = FALSE, EISPACK = FALSE);
  see \code{\link{eigen}} for more details.  

  3.  rNames = rownames(x);  if this is NULL, rNames = if(symmetric)
  paste('x', 1:nrow(x), sep='') else paste('xcol', 1:nrow(x)).  

  4.  Parse 'valuenames', assign to names(ev[['values']]).  

  5.  dimnames(ev[['vectors']]) <- list(rNames, valuenames) 

  NOTE:  This naming convention is fairly obvious if 'x' is symmetric.
  Otherwise, dimensional analysis suggests problems with almost any
  naming convention.  To see this, consider the following simple
  example:

  \deqn{
    X <- matrix(1:4, 2, dimnames=list(LETTERS[1:2], letters[3:4]))
  }
  \tabular{rrr}{
    \tab c \tab d \cr
    A \tab 1 \tab 3 \cr
    B \tab 2 \tab 4 \cr
  }
  \deqn{
    X.inv <- solve(X)
  }
  \tabular{rrr}{
    \tab A \tab B \cr
    c \tab -2 \tab 1.5 \cr
    d \tab 1 \tab -0.5 \cr
  }

  One way of interpreting this is to assume that colnames are really
  reciprocals of the units.  Thus, in this example, X[1,1] is in units
  of 'A/c' and X.inv[1,1] is in units of 'c/A'.  This would make any
  matrix with the same row and column names potentially dimensionless.
  Since eigenvalues are essentially the diagonal of a diagonal matrix,
  this would mean that eigenvalues are dimensionless, and their names
  are merely placeholders.   
}
\value{
  a list with components values and (if only.values = FALSE)
  vectors, as described in \code{\link{eigen}}.  
}
%\references{}
\author{Spencer Graves}
\seealso{
  \code{\link{eigen}},
  \code{\link{svd}}
  \code{\link{qr}}
  \code{\link{chol}}
}
\examples{
X <- matrix(1:4, 2, dimnames=list(LETTERS[1:2], letters[3:4]))
eigen(X)
Eigen(X)
Eigen(X, valuenames='eigval')

Y <- matrix(1:4, 2, dimnames=list(letters[5:6], letters[5:6]))
Eigen(Y)

Eigen(Y, symmetric=TRUE)
# only the lower triangle is used;
# the upper triangle is ignored.  
}
\keyword{array}
