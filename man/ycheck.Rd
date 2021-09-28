\name{ycheck}
\alias{ycheck}
\title{
Check Data Arrays for Smoothing Functions
}
\description{
Check a numeric vector, a matrix or an array of data values for consistency
with the number of argument values.
}
\usage{
ycheck(y, n)
}
\arguments{
  \item{y}{A numeric vector, matrix, or 3D array, with first dimension size 
  \code{n}.}
  \item{n}{The number of argument values for smoothing.}
}
\value{
Argument \code{y} is returned if it checks out, converted to a one-column
matrix if the input was a numeric vector.}
\seealso{
  \code{\link{smooth.basis1}},
  \code{\link{smooth.monotone}},
  \code{\link{smooth.pos}},
  \code{\link{smooth.monotone}},
  \code{\link{smooth.surp}}
}
