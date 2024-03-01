\name{bsplinepen}
\alias{bsplinepen}
\title{
  B-Spline Penalty Matrix
}
\description{
  Computes the matrix defining the roughness penalty for functions
  expressed in terms of a B-spline basis.
}
\usage{
  bsplinepen(basisobj, Lfdobj=2, rng=basisobj$rangeval, returnMatrix=FALSE)
}
\arguments{
  \item{basisobj}{
  a B-spline basis object.
  }
  \item{Lfdobj}{
  either a nonnegative integer or a linear differential operator object.
  }
  \item{rng}{
  a vector of length 2 defining range over which the basis penalty is to
  be computed.
  }
  \item{returnMatrix}{
    logical:  If TRUE,  a two-dimensional is returned using a
    special class from the Matrix package.
  }
}
\value{
a symmetric matrix of order equal to the number of basis functions
defined by the B-spline basis object.  Each element is the inner product
of two B-spline basis functions after applying the derivative or linear
differential operator defined by \code{Lfdobj}.
}
\details{
A roughness penalty for a function $x(t)$ is defined by
integrating the square of either the derivative of $x(t)$ or,
more generally, the result of applying a linear differential operator
$L$ to it.  The most common roughness penalty is the integral of
the square of the second derivative, and
this is the default. To apply this roughness penalty, the matrix of
inner products of the basis functions (possibly after applying the
linear differential operator to them) defining this function
is necessary. This function just calls the roughness penalty evaluation
function specific to the basis involved.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
  \emph{Functional data analysis with R and Matlab}, Springer, New
  York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\examples{
##
## bsplinepen with only one basis function
##
bspl1.1 <- create.bspline.basis(nbasis=1, norder=1)
pen1.1 <- bsplinepen(bspl1.1, 0)

##
## bspline pen for a cubic spline with knots at seq(0, 1, .1)
##
basisobj <- create.bspline.basis(c(0,1),13)
#  compute the 13 by 13 matrix of inner products of second derivatives
penmat <- bsplinepen(basisobj)

##
## with rng of class Date or POSIXct
##
# Date
invasion1 <- as.Date('1775-09-04')
invasion2 <- as.Date('1812-07-12')
earlyUS.Canada <- c(invasion1, invasion2)
BspInvade1 <- create.bspline.basis(earlyUS.Canada)
Binvadmat <- bsplinepen(BspInvade1)

# POSIXct
AmRev.ct <- as.POSIXct1970(c('1776-07-04', '1789-04-30'))
BspRev1.ct <- create.bspline.basis(AmRev.ct)
Brevmat <- bsplinepen(BspRev1.ct)

}
% docclass is function
\keyword{smooth}
