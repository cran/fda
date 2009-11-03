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
##
## Minimal example:  a B-spline of order 1, i.e., a step function
## with 0 interior knots:
##
bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)
m <- getbasismatrix(seq(0, 1, .2), bspl1.1)

# check
m. <- matrix(rep(1, 6), 6,
    dimnames=list(NULL, 'bspl') )
\dontshow{stopifnot(}
all.equal(m, m.)
\dontshow{)}

##
## Use library(Matrix)
##
bspl2.2 <- create.bspline.basis(norder=2, breaks=seq(0, 1, .5))
M2 <- getbasismatrix(seq(0, 1, .2), bspl2.2, returnMatrix=TRUE)

M2. <- Matrix(c(1,.6,.2,0,0,0, 0,.4,.8,.8,.4,0, 0,0,0,.2,.6,1), 6)

(M2-M2.)

\dontshow{stopifnot(}
all(abs(M2-M2.)<1e-15)
\dontshow{)}

##
## Date and POSIXct
##
# Date
July4.1776 <- as.Date('1776-07-04')
Apr30.1789 <- as.Date('1789-04-30')
AmRev <- c(July4.1776, Apr30.1789)
BspRevolution <- create.bspline.basis(AmRev)
AmRevYears <- seq(July4.1776, Apr30.1789, length.out=14)
AmRevMatrix <- getbasismatrix(AmRevYears, BspRevolution)
matplot(AmRevYears, AmRevMatrix, type='b')
# Image is correct, but
# matplot does not recogize the Date class of x

# POSIXct
AmRev.ct <- as.POSIXct1970(c('1776-07-04', '1789-04-30'))
BspRev.ct <- create.bspline.basis(AmRev.ct)
AmRevYrs.ct <- seq(AmRev.ct[1], AmRev.ct[2], length.out=14)
AmRevMat.ct <- getbasismatrix(AmRevYrs.ct, BspRev.ct)
matplot(AmRevYrs.ct, AmRevMat.ct, type='b')
# Image is correct, but
# matplot does not recogize the Date class of x

}
% docclass is function
\keyword{smooth}
