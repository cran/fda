\name{plotFd}
\alias{plotFd}
\title{Plot a Functional Data Object }
\description{
Functional data observations, or a derivative of them, are plotted. These may be either plotted simultaneously, as matplot does for multivariate data, or one by one with a mouse click to move from one plot to another. The function also accepts the other plot specification arguments that the regular plot does. 
}
\usage{
plotFd(fd, Lfd=0, matplt=TRUE, href=TRUE, nx=101, xlab=xlabel, ylab=ylabel, xlim=rangex, ylim=rangey, ...)
}
\arguments{
  \item{fd}{A functional data object. }
  \item{Lfd }{Either a nonnegative integer or a linear differential operator object. If present, the derivative or the value of applying the operator is plotted rather than the functions themselves. }
  \item{matplt }{A logical variable. If TRUE, the observations are all plotted simultaneously. If not, they are plotted one by one, with a mouse click to move to the next plot. }
  \item{href }{A logical variable. If TRUE, a horizontal dotted line is plotted at zero. }
  \item{nx }{The number of points to use to define the plot. The default is usually enough, but for a highly variable function more may be required. }
  \item{xlab }{A label for the abscissa. }
  \item{ylab }{A label for the ordinate. }
  \item{xlim }{Axis limits for the abscissa. }
  \item{ylim }{Axis limits for the ordinate. }
  \item{... }{Additional arguments that can be used in all plotting programs. }
}
\details{
}
\value{
A plot of the functional data. Note that for multivariate data, a suitable array must first be defined using the par function. 
}
\references{}
\note{}
\author{}
\seealso{
linesFd, plotfit.fd 
}
\examples{}
\keyword{}
