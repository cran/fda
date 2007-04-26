\name{plot.fd}
\alias{plot.fd}
\title{
  Plot a Functional Data Object
}
\description{
  Functional data observations, or a derivative of them, are plotted.
  These may be either plotted simultaneously, as \code{matplot} does for
  multivariate data, or one by one with a mouse click to move from one
  plot to another.  The function also accepts the other plot
  specification arguments that the regular \code{plot} does.
}
\usage{
plot.fd(x, y, Lfdobj=0, href=TRUE,
        titles=NULL, xlim=rangex, ylim=rangey,
        xlab=xlabel, ylab=ylabel, ask=FALSE, nx=201, ...)
}
\arguments{
  \item{x}{
    functional data object(s) to be plotted.
  }
  \item{y}{
    sequence of points at which to evaluate the functions 'x' and plot
    on the horizontal axis.  Defaults to seq(rangex[1], rangex[2],
    length = nx).  
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object. If present, the derivative or the value of applying the
    operator is plotted rather than the functions themselves.
  }
  \item{href}{
    a logical variable:  If \code{TRUE}, add a horizontal reference line
    at 0.  
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
    a logical value:  If \code{TRUE}, each curve is shown separately, and
    the plot advances with a mouse click
  }
  \item{nx}{
    the number of points to use to define the plot.  The default is
    usually enough, but for a highly variable function more may be
    required.
  }
  \item{\dots }{
    additional plotting arguments that can be used with function
    \code{plot}
  }
}
\value{
  invisible(NULL)
}
\section{Side Effects}{
  a plot of the functional observations
}
\details{
  Note that for multivariate data, a
  suitable array must first be defined using the \code{par} function. 
}
\seealso{
  \code{\link{lines.fd}}, 
  \code{\link{plotfit.fd}}
}
\examples{

#daytime   <- (1:365)-0.5
#dayrange  <- c(0,365)
#dayperiod <- 365
#nbasis     <- 65
#dayrange  <- c(0,365)

daybasis65 <- create.fourier.basis(c(0, 365), 65)
harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))

harmfdPar     <- fdPar(daybasis65, harmaccelLfd, lambda=1e5)

daytempfd <- with(CanadianWeather, data2fd(
        dailyAv[,,"Temperature.C"], day.5, 
        daybasis65,argnames=list("Day", "Station", "Deg C")) )

#  plot all the temperature functions for the monthly weather data
plot(daytempfd, main="Temperature Functions")

\dontrun{
# To plot one at a time:  
# The following pauses to request page changes.
# (Without 'dontrun', the package build process
# might encounter problems with the par(ask=TRUE)
# feature.)
plot(daytempfd, main="Temperature Functions", ask=TRUE)
}

}
% docclass is function
\keyword{smooth}
\keyword{hplot}
