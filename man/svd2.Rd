\name{svd2}
\alias{svd2}
\title{
  singular value decomposition with automatic error handling
}
\description{
  The 'svd' function in R 2.5.1 occasionally throws an error
  with a cryptic message.  In some such cases, changing the
  LINPACK argument has worked.
}
\usage{
  svd2(x, nu = min(n, p), nv = min(n, p), LINPACK = FALSE)
}
\arguments{
\item{x, nu, nv, LINPACK}{
  as for the 'svd' function in the 'base' package.
  }
}

\value{
  a list with components d, u, and v, as described in the help file for
  'svd' in the 'base' package.
}
\details{
  In R 2.5.1, the 'svd' function sometimes stops with a cryptic error
  message for a matrix x for which a second call to 'svd' with !LINPACK
  will produce an answer.

  When such conditions occur, write a warning and save the matrix with
  the result to a file with a name like 'svd.LINPACK.error.matrix*.rda'
  in the working directory.

  Except for these rare pathologies, 'svd2' should work the same as
  'svd'.

  NOTE:  This function was written to get past a pathology without
  isolating a replicable example.  The problem may have long been
  fixed.  It is maintained in case it has not.
}
\seealso{
  \code{\link{svd}}
}
\examples{
# example from svd
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
X <- hilbert(9)[, 1:6]
s <- svd(X)

s2 <- svd2(X)
\dontshow{stopifnot(}
all.equal(s, s2)
\dontshow{)}
}
\keyword{array}
