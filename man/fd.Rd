\name{fd}
\alias{fd}
\alias{[.fd}
\title{
  Define a Functional Data Object
}
\description{
  This is the constructor function for objects of the \code{fd} class.
  Each function that sets up an object of this class must call this
  function.  This includes functions  \code{\link{smooth.basis}}, \code{\link{density.fd}}, and so forth
  that estimate functional data objects that smooth or otherwise
  represent data.  Ordinarily, users of the functional data analysis
  software will not need to call this function directly, but these notes
  are valuable to understanding the components of a \code{list} of class
  \code{fd}.
}
\usage{
fd(coef=NULL, basisobj=NULL, fdnames=NULL)
}
\arguments{
  \item{coef}{
    a vector, matrix, or three-dimensional array of coefficients.

    The first dimension (or elements of a vector) corresponds to basis
    functions.

    A second dimension corresponds to the number of functional
    observations, curves or replicates.  If \code{coef} is a vector, it
    represents only a single functional observation.

    If \code{coef} is an array, the third dimension corresponds to
    variables for multivariate functional data objects.

    A functional data object is "univariate" if \code{coef} is a vector
    or matrix and "multivariate" if it is a three-dimensional array.

    if(is.null(coef)) coef <- rep(0, basisobj[['nbasis']])
  }
  \item{basisobj}{
    a functional basis object defining the basis

    \code{
      if(is.null(basisobj)){
	if(is.null(coef)) basisobj <- basisfd()
	else {
	  rc <- range(coef)
	  if(diff(rc)==0) rc <- rc+0:1
	  nb <- max(4, nrow(coef))
	  basisobj <- create.bspline.basis(rc, nbasis = nb)
	}
      }
    }
  }
  \item{fdnames}{
    A list of length 3, each member being a string vector containing
    labels for the levels of the corresponding dimension of the discrete
    data.  The first dimension is for argument values, and is given the
    default name "time", the second is for replications, and is given
    the default name "reps", and the third is for functions, and is
    given the default name "values".
  }
}
\value{
  A functional data object (i.e., having class \code{fd}), which is a
  list with components named \code{coefs}, \code{basis}, and
  \code{fdnames}.
}
\details{
  To check that an object is of this class, use function
  \code{is.fd}.

  Normally only developers of new functional data analysis
  functions will actually need to use this function.
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York
}
\seealso{
  \code{\link{smooth.basis}}
  \code{\link{smooth.fdPar}}
  \code{\link{smooth.basisPar}}
  \code{\link{density.fd}}
  \code{\link{create.bspline.basis}}
  \code{\link{arithmetic.fd}}
}
\examples{
##
## default
##
fd()
oldpar <- par(no.readonly=TRUE)
##
## The simplest b-spline basis:  order 1, degree 0, zero interior knots:
##       a single step function
##
bspl1.1    <- create.bspline.basis(norder=1, breaks=0:1)
fd.bspl1.1 <- fd(0, basisobj=bspl1.1)

fd.bspl1.1a <- fd(basisobj=bspl1.1)
\dontshow{ stopifnot( }
all.equal(fd.bspl1.1, fd.bspl1.1a)
\dontshow{ ) }
# TRUE

# the following three lines shown an error in a non-cran check:
# if(!CRAN()) {
#   fd.bspl1.1b <- fd(0)
# }

##
## Cubic spline:  4  basis functions
##
bspl4 <- create.bspline.basis(nbasis=4)
plot(bspl4)
parab4.5 <- fd(c(3, -1, -1, 3)/3, bspl4)
# = 4*(x-.5)^2
plot(parab4.5)

##
## Fourier basis
##
f3 <- fd(c(0,0,1), create.fourier.basis())
plot(f3)
# range over +/-sqrt(2), because
# integral from 0 to 1 of cos^2 = 1/2
# so multiply by sqrt(2) to get
# its square to integrate to 1.

##
## subset of an fd object
##

gaitbasis3 <- create.fourier.basis(nbasis=5)
gaittime = (1:20)/21
gaitfd3    <- smooth.basis(gaittime, gait, gaitbasis3)$fd
gaitfd3[1]
par(oldpar)
}
% docclass is function
\keyword{smooth}
\keyword{internal}
