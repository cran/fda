\name{plotfit}
%\alias{plotfit}
\alias{plotfit.fd}
\alias{plotfit.fdSmooth}
\title{
  Plot a Functional Data Object With Data
}
\description{
  Plot either functional data observations 'x' with a fit 'fdobj' or
  residuals from the fit.

  This function is useful for assessing how well a functional data
  object fits the actual discrete data.

  The default is to make one plot per functional observation with fit
  if residual is FALSE and superimposed lines if residual==TRUE.

  With multiple plots, the system waits to confirm a desire to move to
  the next page unless ask==FALSE.
}
\usage{
plotfit.fd(y, argvals, fdobj, rng = NULL, index = NULL,
      nfine = 101, residual = FALSE, sortwrd = FALSE, titles=NULL,
      ylim=NULL, ask=TRUE, type=c("p", "l")[1+residual],
      xlab=NULL, ylab=NULL, sub=NULL, col=1:9, lty=1:9, lwd=1,
      cex.pch=1, axes=NULL, ...)
plotfit.fdSmooth(y, argvals, fdSm, rng = NULL, index = NULL,
      nfine = 101, residual = FALSE, sortwrd = FALSE, titles=NULL,
      ylim=NULL, ask=TRUE, type=c("p", "l")[1+residual],
      xlab=NULL, ylab=NULL, sub=NULL, col=1:9, lty=1:9, lwd=1,
      cex.pch=1, axes=NULL, ...)
}
\arguments{
  \item{y}{
    a vector, matrix or array containing the discrete observations used
    to estimate the functional data object.
  }
  \item{argvals}{
    a vector containing the argument values corresponding to the first
    dimension of \code{y}.
  }
  \item{fdobj}{
    a functional data object estimated from the data.
  }
  \item{fdSm}{ an object of class \code{fdSmooth} }
  \item{rng}{
    a vector of length 2 specifying the limits for the horizontal axis.
    This must be a subset of fdobj[['basis']][['rangeval']], which is
    the default.
  }
  \item{index}{
    a set of indices of functions if only a subset of the observations
    are to be plotted.  Subsetting can also be achieved by subsetting
    \code{y};  see \code{details}, below.
  }
  \item{nfine}{
    the number of argument values used to define the plot of the
    functional data object.  This may need to be increased if the
    functions have a great deal of fine detail.
  }
  \item{residual}{
    a logical variable:  if \code{TRUE}, the residuals are plotted
    rather than the data and functional data object.
  }
  \item{sortwrd}{
    a logical variable:  if \code{TRUE}, the observations (i.e., second
    dimension of \code{y}) are sorted for plotting by the size of the
    sum of squared residuals.
  }
  \item{titles}{
    a vector containing strings that are titles for each observation.
  }
  \item{ylim}{
    a numeric vector of length 2 giving the y axis limits;  see 'par'.
  }
  \item{ask}{
    If TRUE and if 'y' has more levels than the max length of col, lty,
    lwd and cex.pch, the user must confirm page change before the next
    plot will be created.
  }
  \item{type}{
    type of plot desired, as described with \code{\link{plot}}.  If
    residual == FALSE, 'type' controls the representation for 'x', which
    will typically be 'p' to plot points but not lines;  'fdobj' will
    always plot as a line.  If residual == TRUE, the default type ==
    "l";  an alternative is "b" for both.
  }
  \item{xlab}{x axis label. }
  \item{ylab}{
    Character vector of y axis labels.  If(residual), ylab defaults to
    'Residuals', else to varnames derived from names(fdnames[[3]] or
    fdnames[[3]] or dimnames(y)[[3]].
  }
  \item{sub}{
    subtitle under the x axis label.  Defaults to the RMS residual from
    the smooth.
  }
  \item{col, lty, lwd, cex.pch}{
    Numeric or character vectors specifying the color (col), line type
    (lty), line width (lwd) and size of plotted character symbols
    (cex.pch) of the data representation on the plot.

    If ask==TRUE, the length of the longest of these determines the
    number of levels of the array 'x' in each plot before asking the
    user to acknowledge a desire to change to the next page.  Each of
    these is replicated to that length, so col[i] is used for x[,i] (if
    x is 2 dimensional), with line type and width controlled by lty[i]
    and lwd[i], respectively.

    If ask==FALSE, these are all replicated to length = the number of
    plots to be superimposed.

    For more information on alternative values for these paramters, see
    'col', 'lty', 'lwd', or 'cex' with \code{\link{par}}.
  }
  \item{axes}{
    Either a logical or a list or \code{NULL}.

    \itemize{
      \item{logical}{
	     whether axes should be drawn on the plot
      }
      \item{list}{
	      a list used to create custom \code{axes} used to create axes via
	      \code{x$axes[[1]]} and \code{}x$axes[-1])}.  The primary example of
	this uses \code{list("axesIntervals", ...)}, e.g., with
	\code{Fourier} bases to create \code{CanadianWeather} plots
    }
  }
  \item{\dots }{
    additional arguments such as axis labels that may be used with other
    \code{plot} functions.
  }
}
\details{
  \code{plotfit} plots discrete data along with a functional data object
  for fitting the data.  It is designed to be used after something like
  \code{smooth.fd}, \code{smooth.basis} or
  \code{smooth.fdPar} to check the fit of the data offered by the
  \code{fd} object.

  \code{plotfit.fdSmooth} calls \code{plotfit} for its 'fd' component.

  The plot can be restricted to a subset of observations (i.e., second
  dimension of \code{y}) or variables (i.e., third dimension of
  \code{y}) by providing \code{y} with the dimnames for its second and
  third dimensions matching a subset of the dimnames of fdobj[['coef']]
  (for \code{plotfit.fd} or fdSm[['fdobj']][['coef']] for
  \code{plotfit.fdSmooth}).  If only one observation or variable is to
  be plotted, \code{y} must include 'drop = TRUE', as, e.g., y[, 2, 3,
  drop=TRUE].  If \code{y} or fdobj[['coef']] does not have dimnames on
  its second or third dimension, subsetting is achieved by taking the
  first few colums so the second or third dimensions match.  This is
  achieved using checkDims3(y, fdobj[['coef']], defaultNames =
  fdobj[['fdnames']]]).
}
\value{
  A matrix of mean square deviations from predicted is returned
  invisibly.  If fdobj[["coefs"]] is a 3-dimensional array, this is a
  matrix of dimensions equal to the last two dimensions of
  fdobj[["coefs"]].  This will typically be the case when x is also a
  3-dimensional array with the last two dimensions matching those of
  fdobj[["coefs"]].  The second dimension is typically replications and
  the third different variables.

  If x and fobj[["coefs"]] are vectors or 2-dimensional arrays, they are
  padded to three dimensions, and then MSE is computed as a matrix with
  the second dimension = 1;  if x and fobj[["coefs"]] are vectors, the
  first dimension of the returned matrix will also be 1.
}
\section{Side Effects}{
  a plot of the the data 'x' with the function or the deviations between
  observed and predicted, depending on whether residual is FALSE or
  TRUE.
}

\seealso{
  \code{\link{plot}},
  \code{\link{plot.fd}},
  \code{\link{lines.fd}},
  \code{\link{plot.fdSmooth}},
  \code{\link{lines.fdSmooth}},
  \code{\link{par}},
  \code{\link{smooth.fd}},
  \code{\link{smooth.basis}},
  \code{\link{smooth.basisPar}},
  \code{\link{checkDims3}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
#  set up a Fourier basis for smoothing temperature data
daybasis65 <- create.fourier.basis(c(0, 365), 65,
                    axes=list("axesIntervals"))
#  smooth the average temperature data using function smooth.basis
Daytempfd <- with(CanadianWeather, smooth.basis(day.5,
         dailyAv[,,"Temperature.C"], daybasis65)$fd )
daytempfd <- with(CanadianWeather, smooth.basis(day.5,
         dailyAv[,,"Temperature.C"],
         daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd )
#  Plot the temperature data along with the fit to the data for the first
#  station, St. John's Newfoundland
# If you want only the fitted functions, use plot(daytempfd)
# To plot only a single fit vs. observations, use argument index
# to request which one you want.
with(CanadianWeather, plotfit.fd(dailyAv[, , "Temperature.C",
     drop=FALSE], argvals= day.5, daytempfd, index=1, titles=place) )
# Default ylab = daytempfd[['fdnames']]

with(CanadianWeather, plotfit.fd(dailyAv[, , "Temperature.C"],
     argvals= day.5, Daytempfd, index=1, titles=place) )

plot(daytempfd)

if (!CRAN()) {
# plot all the weather stations, one by one after a click on the plot
# in response to a request.
# This example is within the "dontrun" environment to prevent the
# the R package checking process from pausing: without 'dontrun', the package
# build process might encounter problems with the par(ask=TRUE) feature.
with(CanadianWeather, plotfit.fd(dailyAv[,, "Temperature.C"], day.5,
     daytempfd, ask=TRUE) )
}
#  Now plot results for two weather stations.
op <- par(mfrow=c(2,1), xpd=NA, bty="n")
# xpd=NA:  clip lines to the device region,
#       not the plot or figure region
# bty="n":  Do not draw boxes around the plots.
ylim <- range(CanadianWeather$dailyAv[,,"Temperature.C"])
# Force the two plots to have the same scale
with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
          daytempfd, index=2, titles=place, ylim=ylim) )
with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5,
          daytempfd, index=35, titles=place, ylim=ylim) )
\dontrun{
# Plot residuals with interactive display of stations one by one
par(op)
with(CanadianWeather, plotfit.fd(dailyAv[, , "Temperature.C"],
          day.5, daytempfd, residual=TRUE) )
}
#  The gait data are bivariate, and this code illustrates how plotfit.fd
#  deals with the plotting of two variables at the same time
#  First define normalized times and their range
gaittime  <- (0:19) + 0.5
gaitrange <- c(0,20)
#  Define the harmonic acceleration differential operator
harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
#  Set up basis for representing gait data.
gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)
#  Smooth the data
gaitfd <- smooth.basisPar(gaittime, gait, gaitbasis, 
                          Lfdobj=harmaccelLfd, lambda=1e-2)$fd
#  Assign names to the data
names(gaitfd$fdnames) <- c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] <- c("Hip", "Knee")
\dontrun{
#  plot each pair of curves interactively, two plots per page, the top
#  for hip angle, and the bottom for knee angle
plotfit.fd(gait, gaittime, gaitfd)
#  Plot the residuals, sorting cases by residual sum of squares summed over
#  both hip and knee angles.
#  The first series of 39 plots are for hip angle, two plots per page,
#  and the second 39 are for knee angle.  The plots are sorted by the
#  size of the total residual sum of squares, but RMS residual values
#  for specific angles are not all going to be in order.
plotfit.fd(gait, gaittime, gaitfd, residual=TRUE, sort=TRUE)
}
par(oldpar)
}
\keyword{smooth}
\keyword{hplot}
