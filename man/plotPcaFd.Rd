\name{plotPcaFd}
\alias{plotPcaFd}
\title{Plot Output from Principal Components Analysis Function PCA.FD }
\description{
Function pca.fd produces an object of class pcafd. This function plots the information in this object. 
}
\usage{
plotPcaFd(pcafd, nx=128, pointplot=TRUE, harm=0, expand=0, cycle=FALSE, ...)
}
\arguments{
  \item{pcafd }{An object of class pcafd. }
  \item{nx }{The number of argument values used to plot the harmonics. }
  \item{pointplot }{A logical variable. If pointplot is TRUE, then the harmonics are plotted as + and -, otherwise lines are used. }
  \item{harm }{The number of harmonics to plot }
  \item{expand }{If expand = 0 then the effect of adding +/- 2 standard deviations of each principal component are given. Otherwise the factor expand is used. }
  \item{cycle }{A logical variable. If cycle is TRUE and there are two variables, then a cycle plot will be drawn. If the number of variables is anything else, cycle will be ignored. }
  \item{... }{Additional arguments that can be used in all plotting programs. }
}
\details{
}
\value{
A plot of the results of a principal components analysis. 
}
\references{}
\note{}
\author{}
\seealso{
pca.fd 
}
\examples{
}
\keyword{}
