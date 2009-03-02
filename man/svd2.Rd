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
  will produce an answer.  When such conditions occur, assign 'x' with
  attributes 'nu', 'nv', and 'LINPACK' to '.svd.LINPACK.error.matrix'
  in 'env = .GlobalEnv'.

  Except for these rare pathologies, 'svd2' should work the same as
  'svd'.
}
\keyword{array}
\seealso{
\code{\link{svd}},
}

