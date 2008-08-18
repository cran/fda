\name{fd}
\alias{fd}
\title{
  Define a Functional Data Object
}
\description{
  This is the constructor function for objects of the \code{fd} class.
  Each function that sets up an object of this class must call this
  function.  This includes functions \code{data2fd},
  \code{smooth.basis}, \code{density.fd}, and so forth that estimate
  functional data objects that smooth or otherwise represent data.
  Ordinarily, user of the functional data analysis software will not
  need to call this function directly, but these notes are valuable to
  understanding what the "slots" or "members" of the \code{fd} class
  are.  
}
\usage{
fd(coef=NULL, basisobj=NULL, fdnames=defaultnames)
}
\arguments{
  \item{coef}{
    a vector, matrix, or three-dimensional array of coefficients.  The
    first dimension of matrices and arrays, as well as that of a vector,
    corresponds to basis functions.  The second dimension corresponds to
    the number of functions or curves, or to replicates.  The third
    dimension if present corresponds to variables for multivariate
    functional data objects.  A functional data object is referred to as
    "univariate" if this argument is a vector or a matrix, and
    "multivariate" if it is a three-dimensional array.

    if(is.null(coef)) coef <- rep(0, basisobj[['nbasis']]) 
  }
  \item{basisobj}{
    a functional basis object defining the basis

    if(is.null(basisobj)){
      rc <- range(coef)
      if(diff(rc)==0) rc <- rc+0:1
      nb <- max(4, nrow(coef)) 
      basisobj <- create.bspline.basis(rc, nbasis = nb)
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
  \code{\link{data2fd}} 
  \code{\link{smooth.basis}}
  \code{\link{density.fd}}
  \code{\link{create.bspline.basis}}
}
\examples{
##
## The simplest b-spline basis:  order 1, degree 0, zero interior knots:  
##       a single step function 
##
bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)
fd.bspl1.1 <- fd(0, basisobj=bspl1.1)

fd.bspl1.1a <- fd(basisobj=bspl1.1)
\dontshow{ stopifnot( }
all.equal(fd.bspl1.1, fd.bspl1.1a)
\dontshow{ ) }
# TRUE

\dontrun{
fd.bspl1.1b <- fd(0)
Error in fd(0) : 
  Number of coefficients does not match number of basis functions.

... because fd by default wants to create a cubic spline 
}
##
## Cubic spline:  4  basis functions 
##
bspl4 <- create.bspline.basis(nbasis=4)
plot(bspl4) 
parab4.5 <- fd(c(3, -1, -1, 3)/3, bspl4)
# = 4*(x-.5)^2
plot(parab4.5) 

}
% docclass is function
\keyword{smooth}
\keyword{internal}
