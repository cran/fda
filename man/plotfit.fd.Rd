\name{plotfit.fd}
\alias{plotfit.fd}
\title{Plot a Functional Data Object Along with Data }
\description{
Functional observations are plotted one by one along with the data used to estimate them. 
}
\usage{
plotfit.fd(y, argvals, fd, rng=rangeval, index=1:nrep,
           nfine=101, residual=FALSE, sortwrd=FALSE)
}
\arguments{
  \item{y }{The discrete observations used to estimate the functional data object. }
  \item{argvals }{The argument values corresponding to y. }
  \item{fd }{The functional data object to be plotted. }
  \item{rng}{A vector of two strictly increasing values specifying a range over which the plots are to be shown. This range must be contained within the range used to define the functional data object. }
  \item{index }{A set of indices of functions if only a subset of the observations are to be plotted. }
  \item{nfine }{The number of argument values used to define the plot of the functional data object. }
  \item{residual }{A logical variable. If TRUE, the residuals are plotted rather than the data and functional data object. }
  \item{sortwrd }{A logical variable. If TRUE, the observations are sorted for plotting by the size of the sum of squared residuals. }
}
\details{
}
\value{
A set of plots, where a mouse click advances from one plot to the next. 
}
\references{}
\note{}
\author{}
\seealso{
plotFd, linesFd 
}
\examples{
}
\keyword{}
