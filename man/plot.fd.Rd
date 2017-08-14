\name{plot.fd}
\alias{plot.fd}
\alias{plot.fdPar}
\alias{plot.fdSmooth}
\title{
  Plot a Functional Data Object
}
\description{
  Functional data observations, or a derivative of them, are plotted.
  These may be either plotted simultaneously, as \code{matplot} does
  for multivariate data, or one by one with a mouse click to move
  from one plot to another.  The function also accepts the other plot
  specification arguments that the regular \code{plot} does.  Calling
  \code{plot} with an \code{fdSmooth} or an \code{fdPar}
  object plots its \code{fd} component.
}
\usage{
\method{plot}{fd}(x, y, Lfdobj=0, href=TRUE, titles=NULL,
                    xlim=NULL, ylim=NULL, xlab=NULL,
                    ylab=NULL, ask=FALSE, nx=NULL, axes=NULL, ...)
\method{plot}{fdPar}(x, y, Lfdobj=0, href=TRUE, titles=NULL,
                    xlim=NULL, ylim=NULL, xlab=NULL,
                    ylab=NULL, ask=FALSE, nx=NULL, axes=NULL, ...)
\method{plot}{fdSmooth}(x, y, Lfdobj=0, href=TRUE, titles=NULL,
                    xlim=NULL, ylim=NULL, xlab=NULL,
                    ylab=NULL, ask=FALSE, nx=NULL, axes=NULL, ...)
%plot.fd(x, y, Lfdobj=0, href=TRUE, titles=NULL, xlim=rangex, ylim=rangey,
%        xlab=xlabel, ylab=ylabel, ask=FALSE, nx=201, ...)
}
\arguments{
  \item{x}{
    functional data object(s) to be plotted.
  }
  \item{y}{
    sequence of points at which to evaluate the functions 'x' and
    plot on the horizontal axis.  Defaults to seq(rangex[1],
    rangex[2], length = nx).

    NOTE:  This will be the values on the horizontal axis, NOT the
    vertical axis.
  }
  \item{Lfdobj}{
    either a nonnegative integer or a linear differential operator
    object. If present, the derivative or the value of applying the
    operator is plotted rather than the functions themselves.
  }
  \item{href}{
    a logical variable:  If \code{TRUE}, add a horizontal reference
    line at 0.
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
    a vector of length 2 containing axis limits for the horizontal 
    axis.
  }
  \item{ylim}{
    a vector of length 2 containing axis limits for the vertical 
    axis.
  }
  \item{ask}{
    a logical value:  If \code{TRUE}, each curve is shown separately,
    and the plot advances with a mouse click
  }
  \item{nx}{
    the number of points to use to define the plot.  The default is
    usually enough, but for a highly variable function more may be
    required.
  }
  \item{axes}{
    Either a logical or a list or \code{NULL}.

    \itemize{
      \item{logical}{
      	whether axes should be drawn on the plot
      }
      \item{list}{
	      a list used to create custom \code{axes} used to create axes
	      via \code{do.call(x$axes[[1]], x$axes[-1])}.  The primary
	      example of this uses \code{list("axesIntervals", ...)}, e.g.,
	      with \code{Fourier} bases to create \code{CanadianWeather} 
	      plots
      }
    }
  }
  \item{\dots }{
    additional plotting arguments that can be used with function
    \code{plot}
  }
}
\details{
  Note that for multivariate data, a suitable array must first be
  defined using the \code{par} function.
}
\value{
  'done'
  % This used to return 'invisible(NULL)'.
  % However, with R 2.7.0 under XEmacs with ESS,
  % it sometimes failed to plot.  I changed the return value,
  % and the problem disappeared.
}
\section{Side Effects}{
  a plot of the functional observations
}
\seealso{
  \code{\link{lines.fd}}, \code{\link{plotfit.fd}}
}
\examples{
##
## plot.fd
##

daybasis65 <- create.fourier.basis(c(0, 365), 65,
                    axes=list("axesIntervals"))
%harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
%harmfdPar     <- fdPar(daybasis65, harmaccelLfd, lambda=1e5)

daytempfd <- with(CanadianWeather, Data2fd(day.5,
        dailyAv[,,"Temperature.C"], daybasis65))

#  plot all the temperature functions for the monthly weather data
plot(daytempfd, main="Temperature Functions")

\dontrun{
# To plot one at a time:
# The following pauses to request page changes.
plot(daytempfd, ask=TRUE)
}

# test plot(fd with NA)
daytempfd$coefs[,1] <- NA
plot(daytempfd, main="Temperature Functions")

##
## plot.fdSmooth
##
b3.4 <- create.bspline.basis(norder=3, breaks=c(0, .5, 1))
# 4 bases, order 3 = degree 2 =
# continuous, bounded, locally quadratic
fdPar3 <- fdPar(b3.4, lambda=1)

# Penalize excessive slope Lfdobj=1;
# (Can not smooth on second derivative Lfdobj=2 at it is discontinuous.)
fd3.4s0 <- smooth.basis(0:1, 0:1, fdPar3)

# using plot.fd directly
plot(fd3.4s0$fd)

# same plot via plot.fdSmooth
plot(fd3.4s0)
##
## with Date and POSIXct argvals
##
# Date
invasion1 <- as.Date('1775-09-04')
invasion2 <- as.Date('1812-07-12')
earlyUS.Canada <- c(invasion1, invasion2)
BspInvasion <- create.bspline.basis(earlyUS.Canada)

earlyUSyears <- seq(invasion1, invasion2, length.out=7)
(earlyUScubic <- (as.numeric(earlyUSyears-invasion1)/365.24)^3)
fitCubic <- smooth.basis(earlyUSyears, earlyUScubic, BspInvasion)
plot(fitCubic)

# POSIXct
AmRev.ct <- as.POSIXct1970(c('1776-07-04', '1789-04-30'))
BspRev.ct <- create.bspline.basis(AmRev.ct)
AmRevYrs.ct <- seq(AmRev.ct[1], AmRev.ct[2], length.out=14)
(AmRevLin.ct <- as.numeric(AmRevYrs.ct-AmRev.ct[2]))
fitLin.ct <- smooth.basis(AmRevYrs.ct, AmRevLin.ct, BspRev.ct)
plot(fitLin.ct)

}
% docclass is function
\keyword{smooth}
\keyword{hplot}
