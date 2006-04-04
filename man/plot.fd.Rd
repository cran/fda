\name{plot.fd}
\alias{plot.fd}
\title{
  Plot a Functional Data Object
}
\description{
Functional data observations, or a derivative of them, are plotted.
These may be either plotted
simultaneously, as \code{matplot} does for multivariate data, or one
by one with a mouse click to move from one plot to another.  The
function also accepts the other plot specification arguments that the
regular \code{plot} does.
}
\usage{
plot.fd(x, Lfdobj=0, href=TRUE, nx=201,
        titles=NULL, xlab=xlabel, ylab=ylabel,
        xlim=rangex, ylim=rangey, ask=FALSE, ...)
}
\arguments{
\item{x}{
a functional data object to be plotted.
}
\item{Lfdobj}{
either a nonnegative integer or a linear differential operator object.
If present, the derivative or the value of applying the operator is
plotted rather than the functions themselves.
}
\item{href}{
a logical variable: I
if \code{TRUE}, the observations are all plotted simultaneously.
If not, they are plotted one by one, with a mouse click to move to the
next plot.
}
\item{nx}{
the number of points to use to define the plot.  The default is
usually enough, but for a highly variable function more may be
required.
}
\item{titles}{
a vector of strings for identifying curves
}
\item{xlab}{
a label for the horizontal axis.
}
\item{ylab}{
a label for the vertical axis.
}
\item{xlim}{
a vector of length 2 containing axis limits for the horizontal axis.
}
\item{ylim}{
a vector of length 2 containing axis limits for the vertical axis.
}
\item{ask}{
a logical value:
if \code{TRUE}, each curve is shown separately, and the plot
advances with a mouse click
}
\item{\dots }{
additional plotting arguments that can be used with function
\code{plot}
}
}
\value{
none
}
\section{Side Effects}{
a plot of the functional observations
}
\details{
Note that for multivariate data, a
suitable array must first be defined using the \code{par} function.
}
\seealso{
lines.fd, plotfit.fd
}
\examples{
daytime   <- (1:365)-0.5
dayrange  <- c(0,365)
dayperiod <- 365
nbasis     <- 65
dayrange  <- c(0,365)
daybasis65 <- create.fourier.basis(dayrange, nbasis, dayperiod)
harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), dayrange)
harmfdPar     <- fdPar(daybasis65, harmaccelLfd, 1e5)
daytempfd <- data2fd(daily$tempav, daytime, daybasis65,
                     argnames=list("Day", "Station", "Deg C"))

#  plot all the temperature functions for the monthly weather data
#plot(daytempfd, main="Temperature Functions")
#  plot each temperature function in turn, advancing by a click
#plot(daytempfd, matplt=FALSE)

}
% docclass is function
\keyword{smooth}
