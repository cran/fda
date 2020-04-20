\name{smooth.fdPar}
\alias{smooth.fdPar}
\title{
  Smooth a functional data object using a directly specified roughness
  penalty
}
\description{
  Smooth data already converted to a functional data object, fdobj,
  using directly specified criteria.
}
\usage{
smooth.fdPar(fdobj, Lfdobj=NULL, lambda=1e-4,
             estimate=TRUE, penmat=NULL)
}
\arguments{
  \item{fdobj}{
    a functional data object to be smoothed.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object.

    If \code{NULL}, Lfdobj depends on fdobj[['basis']][['type']]:

    \itemize{
      \item{bspline}{
	Lfdobj <- int2Lfd(max(0, norder-2)), where norder =
	norder(fdobj).
      }
      \item{fourier}{
	Lfdobj = a harmonic acceleration operator:

	\code{Lfdobj <- vec2Lfd(c(0,(2*pi/diff(rng))^2,0), rng)}

	where rng = fdobj[['basis']][['rangeval']].
      }
      \item{anything else}{Lfdobj <- int2Lfd(0)}
    }
  }
  \item{lambda}{
    a nonnegative real number specifying the amount of smoothing
    to be applied to the estimated functional parameter.
  }
  \item{estimate}{
    a logical value:  if \code{TRUE}, the functional parameter is
    estimated, otherwise, it is held fixed.
  }
  \item{penmat}{
    a roughness penalty matrix.  Including this can eliminate the need
    to compute this matrix over and over again in some types of
    calculations.
  }
}
\details{
  1.  fdPar

  2.  smooth.fd
}
\value{
  a functional data object.
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\seealso{
\code{\link{smooth.fd}},
\code{\link{fdPar}},
\code{\link{smooth.basis}},
\code{\link{smooth.pos}},
\code{\link{smooth.morph}}
}
\examples{
	#  see smooth.basis
}
% docclass is function
\keyword{smooth}
