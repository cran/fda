\name{cca.fd}
\alias{cca.fd}
\title{
  Functional Canonical Correlation Analysis
}
\description{
Carry out a functional canonical correlation analysis with
regularization or roughness penalties on the estimated
canonical variables.
}
\usage{
cca.fd(fdobj1, fdobj2=fdobj1, ncan = 2,
       ccafdPar1=fdPar(basisobj1, 2, 1e-10),
       ccafdPar2=ccafdPar1, centerfns=TRUE)
}
\arguments{
\item{fdobj1}{
a functional data object.
}
\item{fdobj2}{
a functional data object.  By default this is \code{ fdobj1 }, in
which case the first argument must be a bivariate functional data
object.
}
\item{ncan}{
the number of canonical variables and weight functions to be
computed.  The default is 2.
}
\item{ccafdPar1}{
a functional parameter object defining the first set of canonical
weight functions.  The object may contain specifications for a
roughness penalty. The default is defined using the same basis
as that used for \code{ fdobj1 } with a slight penalty on its
second derivative.
}
\item{ccafdPar2}{
a functional parameter object defining the second set of canonical
weight functions.  The object may contain specifications for a
roughness penalty. The default is \code{ ccafdParobj1 }.
}
\item{centerfns}{
if TRUE, the functions are centered prior to analysis. This is the
default.
}
}
\value{
an object of class \code{cca.fd} with the 5 slots:

\item{ccwtfd1}{
a functional data object for the first
canonical variate weight function
}
\item{ccwtfd2}{
a functional data object for the second
canonical variate weight function
}
\item{cancorr}{
a vector of canonical correlations
}
\item{ccavar1}{
a matrix of scores on the first canonical variable.
}
\item{ccavar2}{
a matrix of scores on the second canonical variable.
}
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
\code{\link{plot.cca.fd}},
\code{\link{varmx.cca.fd}},
\code{\link{pca.fd}}
}
\examples{
#  Canonical correlation analysis of knee-hip curves

gaittime  <- (1:20)/21
gaitrange <- c(0,1)
gaitbasis <- create.fourier.basis(gaitrange,21)
lambda    <- 10^(-11.5)
harmaccelLfd <- vec2Lfd(c(0, 0, (2*pi)^2, 0))

gaitfd    <- fda::fd(matrix(0,gaitbasis$nbasis,1), gaitbasis)
gaitfdPar <- fda::fdPar(gaitfd, harmaccelLfd, lambda)
gaitfd    <- fda::smooth.basis(gaittime, gait, gaitfdPar)$fd
ccafdPar  <- fda::fdPar(gaitfd, harmaccelLfd, 1e-8)
ccafd0    <- cca.fd(gaitfd[,1], gaitfd[,2], ncan=3, ccafdPar, ccafdPar)
#  display the canonical correlations
round(ccafd0$ccacorr[1:6],3)
#  compute a VARIMAX rotation of the canonical variables
ccafd <- varmx.cca.fd(ccafd0)
#  plot the canonical weight functions
oldpar <- par(no.readonly= TRUE)
plot.cca.fd(ccafd)
par(oldpar)

}
\keyword{smooth}
