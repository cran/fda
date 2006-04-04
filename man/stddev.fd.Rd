\name{stddev.fd}
\alias{stddev.fd}
\title{
  Standard Deviation of Functional Data
}
\description{
Evaluate the standard deviation of a set of functions in a functional
data object.
}
\usage{
stddev.fd(fdobj)
}
\arguments{
\item{fdobj}{
a functional data object.
}
}
\value{
a functional data object with a single replication
that contains the standard deviation of the one or several functions in
the object \code{fdobj}.
}
\seealso{
std.fd, mean.fd, sum.fd, center.fd
}
\examples{
liptime  <- seq(0,1,.02)
liprange <- c(0,1)

#  -------------  create the fd object -----------------
#       use 31 order 6 splines so we can look at acceleration

nbasis <- 51
norder <- 6
lipbasis <- create.bspline.basis(liprange, nbasis, norder)

#  ------------  apply some light smoothing to this object  -------

Lfdobj   <- int2Lfd(4)
lambda   <- 1e-12
lipfdPar <- fdPar(lipbasis, Lfdobj, lambda)

lipfd <- smooth.basis(liptime, lip, lipfdPar)$fd
names(lipfd$fdnames) = c("Normalized time", "Replications", "mm")

lipstdfd <- stddev.fd(lipfd)
plot(lipstdfd)

}
% docclass is function
\keyword{smooth}
