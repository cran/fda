\name{getbasismatrix}
\alias{getbasismatrix}
\title{
  Values of Basis Functions or their Derivatives
}
\description{
Evaluate a set of basis functions or their derivatives at
a set of argument values.
}
\usage{
getbasismatrix(evalarg, basisobj, nderiv=0, returnMatrix=FALSE)
}
\arguments{
  \item{evalarg}{
    a vector of arguments values.
  }
  \item{basisobj}{
    a basis object.
  }
  \item{nderiv}{
    a nonnegative integer specifying the derivative to be evaluated.
  }
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
}
\value{
  a matrix of basis function or derivative values.  Rows correspond
  to argument values and columns to basis functions.
}
\seealso{
  \code{\link{eval.fd}}
}
\examples{
# Minimal example:  a B-spline of order 1, i.e., a step function
# with 0 interior knots:
bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)
m <- getbasismatrix(seq(0, 1, .2), bspl1.1)

\dontshow{stopifnot(}
all.equal(m, matrix(rep(1, 6), 6))
\dontshow{)}

# Use library(Matrix)
bspl2.2 <- create.bspline.basis(norder=2, breaks=seq(0, 1, .5))
M2 <- getbasismatrix(seq(0, 1, .2), bspl2.2, returnMatrix=TRUE)

M2. <- Matrix(c(1,.6,.2,0,0,0, 0,.4,.8,.8,.4,0, 0,0,0,.2,.6,1), 6)

\dontshow{stopifnot(}
all.equal(M2, M2.)
\dontshow{)}

}
% docclass is function
\keyword{smooth}
