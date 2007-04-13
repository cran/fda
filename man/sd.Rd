\name{sd.fd}
\alias{sd.fd}
\alias{std.fd}
\alias{stdev.fd}
\alias{stddev.fd}
\title{
  Standard Deviation of Functional Data
}
\description{
  Evaluate the standard deviation of a set of functions in a functional
  data object.
}
\usage{
sd.fd(fdobj)
std.fd(fdobj)
stdev.fd(fdobj)
stddev.fd(fdobj)
}
\arguments{
  \item{fdobj}{
    a functional data object.
  }
}
\details{
  The multiple aliases are provided for compatibility with previous
  versions and with other languages.  The name for the standard
  deviation function in R is 'sd'.  Matlab uses 'std'.  S-Plus and
  Microsoft Excal use 'stdev'.  'stddev' was used in a previous version
  of the 'fda' package and is retained for compatibility.  
}
\value{
  a functional data object with a single replication
  that contains the standard deviation of the one or several functions in
  the object \code{fdobj}.
}
\seealso{
  \code{\link{mean.fd}}, 
  \code{\link{sum.fd}}, 
  \code{\link{center.fd}}
}
\examples{
liptime  <- seq(0,1,.02)
liprange <- c(0,1)

#  -------------  create the fd object -----------------
#       use 31 order 6 splines so we can look at acceleration

nbasis <- 51
norder <- 6
lipbasis <- create.bspline.basis(liprange, nbasis, norder)
lipbasis <- create.bspline.basis(liprange, nbasis, norder)

#  ------------  apply some light smoothing to this object  -------

Lfdobj   <- int2Lfd(4)
lambda   <- 1e-12
lipfdPar <- fdPar(lipbasis, Lfdobj, lambda)

lipfd <- smooth.basis(liptime, lip, lipfdPar)$fd
names(lipfd$fdnames) = c("Normalized time", "Replications", "mm")

lipstdfd <- sd.fd(lipfd)
plot(lipstdfd)

all.equal(lipstdfd, std.fd(lipfd))
all.equal(lipstdfd, stdev.fd(lipfd))
all.equal(lipstdfd, stddev.fd(lipfd))

}
% docclass is function
\keyword{smooth}
