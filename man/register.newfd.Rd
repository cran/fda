\name{register.newfd}
\alias{register.newfd}
\title{
  Register Functional Data Objects with Pre-Computed Warping Functions
}
\description{ This function registers a new functional data object to pre-computed
warping functions.
}
\usage{
register.newfd(yfd, Wfd,type=c('direct','monotone','periodic'))
}
\arguments{
  \item{yfd}{
    a multivariate functional data object defining the functions to be
    registered with \code{Wfd}.
  }
  \item{Wfd}{
    a functional data object defining the registration functions to be used
    to register \code{yfd}. This can be the result of either \code{landmarkreg}
    or \code{register.fd}.
  }
  \item{type}{
    indicates the type of registration function.
    \itemize{
      \item{direct}{
	assumes \code{Wfd} is a direct definition of the registration
	functions. This is produced by \code{landmarkreg}.
      }
      \item{monotone}{
	assumes that \code{Wfd} defines a monotone functional data
	objected, up to shifting and scaling to make endpoints
	agree. This is produced by \code{register.fd}.
      }
      \item{periodic}{
	does shift registration for periodic functions. This is output
	from \code{register.fd} if \code{periodic=TRUE}.
      }
    }
  }
}
\value{
  a functional data object defining the registered curves.
}
\details{
  Only shift registration is considered for the periodic case.
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York, ch. 6 & 7.
}
\seealso{
  \code{\link{landmarkreg}},
  \code{\link{register.fd}}
}
\examples{

# Register the lip data with landmark registration, then register the first
# derivatives with the pre-computed warping functions.

# Lip data:
lipfd <- smooth.basisPar(liptime, lip, 6, Lfdobj=int2Lfd(4),
                         lambda=1e-12)$fd
names(lipfd$fdnames) <- c("time(seconds)", "replications", "mm")

# Landmark Registration:

lipmeanmarks <- apply(lipmarks,2,mean)

wnbasis <- 6
wnorder <- 4
wbreaks <- c(0,lipmeanmarks,0.35)

warpbasis <- create.bspline.basis(nbasis=wnbasis, norder=wnorder,
                                  breaks=wbreaks);
WfdPar    <- fdPar(fd(basisobj=warpbasis), 2, 1e-4)

lipreglist <- landmarkreg(lipfd, as.matrix(lipmarks), lipmeanmarks, WfdPar)

# And warp:

Dlipfd = deriv.fd(lipfd,Lfdobj=1)

Dlipregfd = register.newfd(Dlipfd,lipreglist$warpfd,type='direct')
}
\keyword{smooth}
