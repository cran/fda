\name{varmx.cca.fd}
\alias{varmx.cca.fd}
\title{
  Rotation of Functional Canonical Components with VARIMAX
}
\description{
Results of canonical correlation analysis are often easier to interpret if
they are rotated.  Among the many possible ways in which this rotation can be
defined, the VARIMAX criterion seems to give satisfactory results most
of the time.
}
\usage{
varmx.cca.fd(ccafd, nx=201)
}
\arguments{
\item{ccafd}{
an object of class "cca.fd" that is produced by function
\code{cca.fd}.
}
\item{nx}{
the number of points in a fine mesh of points that is
required to approximate canonical variable functional
data objects.
}
}
\value{
  a rotated version of argument \code{cca.fd}.
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
  \code{\link{varmx}}, 
  \code{\link{varmx.pca.fd}}
}
% docclass is function
\keyword{smooth}
