\name{eval.basis}
\alias{eval.basis}
\alias{predict.basisfd}
\title{
  Values of Basis Functions or their Derivatives
}
\description{
  A set of basis functions are evaluated at a vector of argument values.
  If a linear differential object is provided, the  values are the
  result of applying the the operator to each basis function.
}
\usage{
eval.basis(evalarg, basisobj, Lfdobj=0, returnMatrix=FALSE)
\method{predict}{basisfd}(object, newdata=NULL, Lfdobj=0,
                          returnMatrix=FALSE, ...)
}
\arguments{
  \item{evalarg, newdata}{
    a vector of argument values at which the basis functiona is to be
    evaluated.
  }
  \item{basisobj}{
    a basis object defining basis functions whose values
    are to be computed.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential.
    operator object.
  }
  \item{object}{
    an object of class \code{basisfd}
  }
  \item{\dots}{
    optional arguments for \code{predict}, not currently used
  }
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
}
\details{
  If a linear differential operator object is supplied, the basis must
  be such that the highest order derivative can be computed. If a
  B-spline basis is used, for example, its order must be one larger than
  the highest order of derivative required.
}
\value{
  a matrix of basis function values with rows corresponding
  to argument values and columns to basis functions.

  \code{predict.basisfd} is a convenience wrapper for
  \code{eval.basis}.
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York
}
\seealso{
  \code{\link{getbasismatrix}},
  \code{\link{eval.fd}},
  \code{\link{plot.basisfd}}
}
\examples{
##
## 1.  B-splines
##
# The simplest basis currently available:
# a single step function
bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)
eval.bspl1.1 <- eval.basis(seq(0, 1, .2), bspl1.1)

# check
eval.bspl1.1. <- matrix(rep(1, 6), 6,
       dimnames=list(NULL, 'bspl') )
\dontshow{stopifnot(}
all.equal(eval.bspl1.1, eval.bspl1.1.)
\dontshow{)}

# The second simplest basis:
# 2 step functions, [0, .5], [.5, 1]
bspl1.2 <- create.bspline.basis(norder=1, breaks=c(0,.5, 1))
eval.bspl1.2 <- eval.basis(seq(0, 1, .2), bspl1.2)

# Second order B-splines (degree 1:  linear splines)
bspl2.3 <- create.bspline.basis(norder=2, breaks=c(0,.5, 1))
eval.bspl2.3 <- eval.basis(seq(0, 1, .1), bspl2.3)
# 3 bases:  order 2 = degree 1 = linear
# (1) line from (0,1) down to (0.5, 0), 0 after
# (2) line from (0,0) up to (0.5, 1), then down to (1,0)
# (3) 0 to (0.5, 0) then up to (1,1).

##
## 2.  Fourier
##
# The false Fourier series with 1 basis function
falseFourierBasis <- create.fourier.basis(nbasis=1)
eval.fFB <- eval.basis(seq(0, 1, .2), falseFourierBasis)

# Simplest real Fourier basis with 3 basis functions
fourier3 <- create.fourier.basis()
eval.fourier3 <- eval.basis(seq(0, 1, .2), fourier3)

# 3 basis functions on [0, 365]
fourier3.365 <- create.fourier.basis(c(0, 365))
eval.F3.365 <- eval.basis(day.5, fourier3.365)
oldpar <- par(no.readonly=TRUE)
matplot(eval.F3.365, type="l")

# The next simplest Fourier basis (5  basis functions)
fourier5 <- create.fourier.basis(nbasis=5)
eval.F5 <- eval.basis(seq(0, 1, .1), fourier5)
matplot(eval.F5, type="l")

# A more complicated example
dayrng <- c(0, 365)

nbasis <- 51
norder <- 6

weatherBasis <- create.fourier.basis(dayrng, nbasis)
basisMat <- eval.basis(day.5, weatherBasis)

matplot(basisMat[, 1:5], type="l")

##
## 3.  predict.basisfd
##
basisMat. <- predict(weatherBasis, day.5)
\dontshow{stopifnot(}
all.equal(basisMat, basisMat.)
\dontshow{)}

##
## 4.  Date and POSIXct
##
# Date
July4.1776 <- as.Date('1776-07-04')
Apr30.1789 <- as.Date('1789-04-30')
AmRev <- c(July4.1776, Apr30.1789)
BspRevolution <- create.bspline.basis(AmRev)
AmRevYears <- seq(July4.1776, Apr30.1789, length.out=14)
AmRevBases <- predict(BspRevolution, AmRevYears)
matplot(AmRevYears, AmRevBases, type='b')
# Image is correct, but
# matplot does not recogize the Date class of x

# POSIXct
AmRev.ct <- as.POSIXct1970(c('1776-07-04', '1789-04-30'))
BspRev.ct <- create.bspline.basis(AmRev.ct)
AmRevYrs.ct <- seq(AmRev.ct[1], AmRev.ct[2], length.out=14)
AmRevBas.ct <- predict(BspRev.ct, AmRevYrs.ct)
matplot(AmRevYrs.ct, AmRevBas.ct, type='b')
# Image is correct, but
# matplot does not recognize the POSIXct class of x
par(oldpar)
}
% docclass is function
\keyword{smooth}
