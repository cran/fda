\name{plot.cca.fd}
\alias{plot.cca.fd}
\title{
  Plot Functional Canonical Correlation Weight Functions
}
\description{
  A canonical correlation analysis produces a series of pairs of functional 
  data objects which, when used as weighting functions, successively maximize
  the corresponding canonical correlation between two functional data objects.
  Like functional principal component weight functions, successive weight
  within either side fo the pair are required to be orthogonal to all previous
  weight functions.  Consequently, each successive canonical correlation will
  no larger than its predecessor, and more likely substantially smaller.
  This function plots an object of class \code{cca.fd} that results from the
  use of function \code{cca.fd}.  Each pair of weight functions is plotted
  after a left mouse click indicating that you are ready for the next plot.
}
\usage{
\method{plot}{cca.fd}(x, cexval = 1, ...)
}
\arguments{
  \item{x}{
    an object of class \code{cca.fd} produced by an invocation of function
    \code{cca.fd.R}.
  }
  \item{cexval}{
    A number used to determine label sizes in the plots.
  }
  \item{\dots}{
    other arguments for 'plot'.
  }
}
\details{
  Produces a plot of a pair of weight functions corresponding to each
  canonical correlation between two functional data objects.
}
\value{
  invisible(NULL)
}
\seealso{
  \code{\link{cca.fd}},
  \code{\link{pda.fd}}
  \code{\link{plot.pca.fd}}
}
\examples{

#  Canonical correlation analysis of knee-hip curves

gaittime  <- (1:20)/21
gaitrange <- c(0,1)
gaitbasis <- create.fourier.basis(gaitrange,21)
lambda    <- 10^(-11.5)
harmaccelLfd <- vec2Lfd(c(0, 0, (2*pi)^2, 0))
gaitfdPar <- fdPar(gaitbasis, harmaccelLfd, lambda)
gaitfd    <- smooth.basis(gaittime, gait, gaitfdPar)$fd
ccafdPar  <- fdPar(gaitfd, harmaccelLfd, 1e-8)
ccafd0    <- cca.fd(gaitfd[,1], gaitfd[,2], ncan=3, ccafdPar, ccafdPar)
#  display the canonical correlations
round(ccafd0$ccacorr[1:6],3)
#  plot the unrotated canonical weight functions
plot.cca.fd(ccafd0)
#  compute a VARIMAX rotation of the canonical variables
ccafd <- varmx.cca.fd(ccafd0)
#  plot the rotated canonical weight functions
plot.cca.fd(ccafd)

}
% docclass is function
\keyword{smooth}
