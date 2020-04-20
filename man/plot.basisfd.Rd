\name{plot.basisfd}
\alias{plot.basisfd}
\title{
  Plot a Basis Object
}
\description{
  Plots all the basis functions.
}
\usage{
\method{plot}{basisfd}(x, knots=TRUE, axes=NULL, ...)
}
\arguments{
  \item{x}{
    a basis object
  }
  \item{knots}{
    logical:  If TRUE and x[['type']] == 'bspline', the knot locations
    are plotted using vertical dotted, red lines.  Ignored otherwise.
  }
  \item{axes}{
    Either a logical or a list or \code{NULL}.

    \itemize{
      \item{logical}{
	whether axes should be drawn on the plot
      }
      \item{list}{
	a list used to create custom \code{axes} used to create axes via
	\code{x$axes[[1]]} and  \code{x$axes[-1]}.  The primary example of
	this uses \code{list("axesIntervals", ...)}, e.g., with
	\code{Fourier} bases to create \code{CanadianWeather} plots
      }
    }
  }
  \item{\dots }{
    additional plotting parameters passed to \code{matplot}.
  }
}
\value{
  none
}
\section{Side Effects}{
  a plot of the basis functions
}
\seealso{
\code{\link{plot.fd}}
}
\examples{
##
## 1.  b-spline
##
# set up the b-spline basis for the lip data, using 23 basis functions,
#   order 4 (cubic), and equally spaced knots.
#  There will be 23 - 4 = 19 interior knots at 0.05, ..., 0.95
lipbasis <- create.bspline.basis(c(0,1), 23)
# plot the basis functions
plot(lipbasis)

##
## 2.  Fourier basis
##
yearbasis3 <- create.fourier.basis(c(0,365),
                    axes=list("axesIntervals") )
#  plot the basis
plot(yearbasis3)

##
## 3.  With Date and POSIXct rangeval
##
# Date
July4.1776 <- as.Date('1776-07-04')
Apr30.1789 <- as.Date('1789-04-30')
AmRev <- c(July4.1776, Apr30.1789)
BspRevolution <- create.bspline.basis(AmRev)
plot(BspRevolution)

# POSIXct
July4.1776ct <- as.POSIXct1970('1776-07-04')
Apr30.1789ct <- as.POSIXct1970('1789-04-30')
AmRev.ct <- c(July4.1776ct, Apr30.1789ct)
BspRev.ct <- create.bspline.basis(AmRev.ct)
plot(BspRev.ct)

}
% docclass is function
\keyword{smooth}
