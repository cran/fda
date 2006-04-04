\name{varmx.pca.fd}
\alias{varmx.pca.fd}
\title{
  Rotation of Functional Principal Components with VARIMAX
Criterion
}
\description{
Principal components are often easier to interpret if they are
rotated.  Among the many possible ways in which this rotation can be
defined, the VARIMAX criterion seems to give satisfactory results most
of the time.
}
\usage{
varmx.pca.fd(pcafd, nharm=scoresd[2], nx=501)
}
\arguments{
\item{pcafd}{
an object of class \code{pca.fd} that is produced by function
\code{pca.fd}.
}
\item{nharm}{
the number of harmonics or principal components to be
rotated.
}
\item{nx}{
the number of argument values in a fine mesh
used to define the harmonics to be
rotated.
}
}
\value{
a rotated principal components analysis object of class
\code{pca.fd}.
}
\seealso{
\code{\link{varmx}}, 
\code{\link{varmx.cca.fd}}
}
% docclass is function
\keyword{smooth}
