\name{axisIntervals}
\alias{axisIntervals}
\alias{axesIntervals}
\title{
  Mark Intervals on a Plot Axis
}
\description{
  Adds an axis (axisintervals) or two axes (axesIntervals)
  to the current plot with tick marks delimiting interval
  described by labels
}
\usage{
axisIntervals(side=1, atTick1=fda::monthBegin.5, atTick2=fda::monthEnd.5,
      atLabels=fda::monthMid, labels=month.abb, cex.axis=0.9, ...)
axesIntervals(side=1:2, atTick1=fda::monthBegin.5, atTick2=fda::monthEnd.5,
      atLabels=fda::monthMid, labels=month.abb, cex.axis=0.9, las=1, ...)
}
\arguments{
  \item{side}{
    an integer specifying which side of the plot the axis is to
    be drawn on.  The axis is placed as follows: 1=below, 2=left,
    3=above and 4=right.
  }
  \item{atTick1}{
    the points at which tick-marks marking the starting points of the
    intervals are to be drawn.  This defaults to 'monthBegin.5' to mark
    monthly periods for an annual cycle.  These are constructed by
    calling axis(side, at=atTick1, labels=FALSE, ...).  For more detail
    on this, see 'axis'.
  }
  \item{atTick2}{
    the points at which tick-marks marking the ends of the
    intervals are to be drawn.  This defaults to 'monthEnd.5' to mark
    monthly periods for an annual cycle.  These are constructed by
    calling axis(side, at=atTick2, labels=FALSE, ...).  Use atTick2=NA
    to rely only on atTick1.  For more detail
    on this, see 'axis'.
  }
  \item{atLabels}{
    the points at which 'labels' should be typed.  These are constructed
    by calling axis(side, at=atLabels, tick=FALSE, ...).  For more detail
    on this, see 'axis'.
  }
  \item{labels}{
    Labels to be typed at locations 'atLabels'.  This is accomplished by
    calling axis(side, at=atLabels, labels=labels, tick=FALSE, ...).
    For more detail on this, see 'axis'.
  }
  \item{cex.axis}{
    Character expansion (magnification) used for axis annotations
    ('labels' in this function call) relative
    to the current setting of 'cex'.  For more detail, see 'par'.
  }
  \item{las}{
    line axis style;  see \code{par}.
  }
  \item{\dots }{
    additional arguments passed to \code{axis}.
  }
}

\value{
  The value from the third (labels) call to 'axis'.  This function is
  usually invoked for its side effect, which is to add an axis to an
  already existing plot.

  \code{axesIntervals} calls \code{axisIntervals(side[1], ...)} then
  \code{axis(side[2], ...)}.
}
\section{Side Effects}{
  An axis is added to the current plot.
}
\author{ Spencer Graves  }
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
  \emph{Functional data analysis with R and Matlab}, Springer, New
  York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{axis}},
  \code{\link{par}}
  \code{\link{monthBegin.5}}
  \code{\link{monthEnd.5}}
  \code{\link{monthMid}}
  \code{\link{month.abb}}
  \code{\link{monthLetters}}
}
\examples{
oldpar <- par(no.readonly= TRUE)
daybasis65 <- create.fourier.basis(c(0, 365), 65)

daytempfd <- with(CanadianWeather, smooth.basis(
       day.5,  dailyAv[,,"Temperature.C"], 
       daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd )

with(CanadianWeather, plotfit.fd(
      dailyAv[,,"Temperature.C"], argvals=day.5,
          daytempfd, index=1, titles=place, axes=FALSE) )
# Label the horizontal axis with the month names
axisIntervals(1)
axis(2)
# Depending on the physical size of the plot,
# axis labels may not all print.
# In that case, there are 2 options:
# (1) reduce 'cex.lab'.
# (2) Use different labels as illustrated by adding
#     such an axis to the top of this plot

with(CanadianWeather, plotfit.fd(
      dailyAv[,,"Temperature.C"], argvals=day.5,
          daytempfd, index=1, titles=place, axes=FALSE) )
# Label the horizontal axis with the month names
axesIntervals()

axisIntervals(3, labels=monthLetters, cex.lab=1.2, line=-0.5)
# 'line' argument here is passed to 'axis' via '...'
par(oldpar)
}
% docclass is function
\keyword{smooth}
\keyword{hplot}
