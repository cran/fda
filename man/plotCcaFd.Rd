\name{plotCcaFd}
\alias{plotCcaFd}
\title{Plot Output from Canonical Correlation Function CCA.FD }
\description{
Function cca.fd produces an object of class ccafd. This function plots the information in this object. 
}
\usage{
plotCcaFd(ccafd, overplt=FALSE, jcan=0, flip=FALSE, ...)
}
\arguments{
  \item{ccafd }{An object of class cca.fd that is produced by function cca.fd. }
  \item{overplt }{A logical variable. If overplt is TRUE, then each pair of weight functions is plotted in a single plot. The line types and colours of the "x" and "y" curves respectively are specified as in function plotFd. }
  \item{jcan }{If jcan=0, then all the pairs of variates are plotted. Otherwise only the variates in vector jcan are plotted (eg if jcan=1, only the leading variate is plotted, if jcan=c(1,3) only the first and third.) }
  \item{flip }{A logical vector of length ncan. If flip[j] is TRUE then the jth pair of weight functions are multiplied by -1. If flip is a scalar it is replicated to the necessary length. }
  \item{... }{Additional arguments that can be used in all plotting programs. }
}
\details{
}
\value{
}
\references{}
\note{}
\author{}
\seealso{}
\examples{
}
\keyword{}
