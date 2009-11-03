\name{fd2list}
\alias{fd2list}
\title{
  Convert a univariate functional data object to s list
}
\description{
  Convert a univariate functional data object to a list for input to
  \code{\link{Lfd}}.
}
\usage{
fd2list(fdobj)
}
\arguments{
  \item{fdobj}{
    a univariate functional data object.
  }
}
\value{
  a list as required for the second argument of \code{\link{Lfd}}.
}
\seealso{
\code{\link{Lfd}}
}
\examples{
Lbasis  = create.constant.basis(c(0,365));  #  create a constant basis
Lcoef   = matrix(c(0,(2*pi/365)^2,0),1,3)   #  set up three coefficients
wfdobj  = fd(Lcoef,Lbasis)      # define an FD object for weight functions
wfdlist = fd2list(wfdobj)       # convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdlist)  #  define the operator object
}
% docclass is function
\keyword{smooth}
