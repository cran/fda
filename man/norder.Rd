\name{norder}
\alias{norder}
\alias{norder.fd}
\alias{norder.basisfd}
\alias{norder.default}
\alias{norder.bspline}
\title{
  Order of a B-spline
}
\description{
  norder = number of basis functions minus the number of interior
  knots.
}
\usage{
norder(x, ...)
\method{norder}{fd}(x, ...)
\method{norder}{basisfd}(x, ...)
\method{norder}{default}(x, ...)

#norder.bspline(x, ...)

#NOTE:  The following is required by CRAN rules that
# function names like "as.numeric" must follow the documentation
# standards for S3 generics, even when they are not.
# Please ignore the following line:
\method{norder}{bspline}(x, ...)
}
\arguments{
  \item{x}{
    Either a basisfd object or an object containing a basisfd object as
    a component.
  }
  \item{\dots}{optional arguments currently unused}
}
\details{
  norder throws an error of basisfd[['type']] != 'bspline'.
}
\value{
  An integer giving the order of the B-spline.
}
%\references{}
\author{
  Spencer Graves
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{create.bspline.basis}}
}
\examples{
  bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)

  stopifnot(norder(bspl1.1)==1)

  stopifnot(norder(fd(0, basisobj=bspl1.1))==1)

  stopifnot(norder(fd(rep(0,4)))==4)

  stopifnot(norder(list(fd(rep(0,4))))==4)
  \dontrun{
    norder(list(list(fd(rep(0,4)))))
    Error in norder.default(list(list(fd(rep(0, 4))))) :
    input is not a 'basisfd' object and does not have a 'basisfd'
    component.
  }

  stopifnot(norder(create.bspline.basis(norder=1, breaks=c(0,.5, 1))) == 1)

  stopifnot(norder(create.bspline.basis(norder=2, breaks=c(0,.5, 1))) == 2)

  # Default B-spline basis:  Cubic spline:  degree 3, order 4,
  # 21 breaks, 19 interior knots.
  stopifnot(norder(create.bspline.basis()) == 4)

  # these five lines throw an error of for a nocran check 
  # if (!CRAN()) {
  #   norder(create.fourier.basis(c(0,12) ))
    # Error in norder.bspline(x) :
    # object x is of type = fourier;  'norder' is only defined for type = 'bspline'
  # }

}
\keyword{smooth}
