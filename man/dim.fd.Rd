\name{dim.fd}
\alias{dim.fd}
\title{
  Dimension of a functional data object 
}
\description{
  This returns the dimension of a functional data object, defined as the
  number of fuctional observations, optionally with the number of
  response variables.  
}
\usage{
dim.fd(x)
}
\arguments{
  \item{x}{
    a functional data object for which the dimenions are desired.  
  }
}
\value{
  a numeric vector of length 1 or 2.  
}
\details{
  If \code{x} is NOT a functional data object, \code{dim.fd} returns
  dim(x).  
}
\seealso{
  \code{\link{dim}} 
  \code{\link{fd}}
}
\examples{
##
## The simplest b-spline basis:  order 1, degree 0, zero interior knots:  
##       a single step function 
##
bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)
fd.bspl1.1 <- fd(0, basisobj=bspl1.1)

dim.fd(fd.bspl1.1) 

##
## handwrit
##

\dontrun{
# This took 132 seconds on an IBM T60 with 1GB RAM:  
HandWr <- smooth.basisPar(1:nrow(handwrit), handwrit)

dim.fd(HandWr$fd)
}

}
% docclass is function
\keyword{smooth}
