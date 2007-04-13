\name{plotfit.fd}
\alias{plotfit.fd}
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
plotfit.fd(x, argvals, fdobj, rng = rangeval,
           index = 1:nrep, nfine = 101, residual = FALSE,
           sortwrd = FALSE, titles=NULL, ylim=NULL,
           ask=FALSE, type=c("p", "l")[1+residual],
           col=1:9, lty=1:9, lwd=1, cex.pch=1, 
           ...)
}
\arguments{
  \item{x}{
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
  \item{rng}{
    a vector of length 2 specifying a reduced range contained within the
    range defined by the functional data object or the arguments.  This
    is used to examine the fit over a limited segment of the whole range
    of it.  
  }
  \item{index}{
    a set of indices of functions if only a subset of the observations
    are to be plotted. 
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
    a logical variable:  if \code{TRUE}, the observations are sorted for
    plotting by the size of the sum of squared residuals.  
  }
  \item{titles}{
    a vector containing strings that are titles for each observation.
  }
  \item{ylim}{
    a numeric vector of length 2 giving the y axis limits;  see 'par'.  
  }
  \item{ask}{
    If TRUE and if 'x' has more levels than the max length of col, lty,
    lwd and cex.pch, the user must confirm page change before the next
    plot will be created.  The number of lines / cases per plot is
    controled by the maximum of the lengths of col, lty, lwd and cex.pch
    (see below), except that if length(dim(x))>2, a separate plot is
    created for each level of the third dimension of x.  
  }
  \item{type}{
    type of plot desired, as described with \code{\link{plot}}.  If
    residual == FALSE, 'type' controls the representation for 'x', which
    will typically be 'p' to plot points but not lines;  'fdobj' will
    always plot as a line.  If residual == TRUE, the default type ==
    "l";  an alternative is "b" for both.  
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
  \item{\dots }{
    additional arguments such as axis labels that may be used with other
    \code{plot} functions. 
  }
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
\details{
  'plotfit.fd' plots discrete data along with a functional data object for
  fitting the data.  It is designed to be used after 'data2fd',
  'smooth.fd', or 'smooth.basis' to check the fit of the data offered by
  the 'fd' object. 
}

}
\seealso{
  \code{\link{plot.fd}}, 
  \code{\link{plot}}, 
  \code{\link{lines.fd}}
  \code{\link{par}}
  \code{\link{data2fd}}
  \code{\link{smooth.fd}}
  \code{\link{smooth.basis}}
}
\examples{
daybasis65 <- create.fourier.basis(c(0, 365), 65)

daytempfd <- with(CanadianWeather, data2fd(
       dailyAv[,,"Temperature.C"], day.5, 
       daybasis65, argnames=list("Day", "Station", "Deg C")) )
 
with(CanadianWeather, plotfit.fd(x=dailyAv[, , "Temperature.C"],
          argvals= day.5, daytempfd, index=1, titles=place, axes=FALSE) )
# Label the horizontal axis with the month names
axis(1, monthBegin.5, labels=FALSE)
axis(1, monthEnd.5, labels=FALSE)
axis(1, monthMid, monthLetters, tick=FALSE)
axis(2)

\dontrun{
# The following pauses to request page changes.
# (Without 'dontrun', the package build process
# might encounter problems with the par(ask=TRUE)
# feature.)
with(CanadianWeather, plotfit.fd(
          dailyAv[,, "Temperature.C"], day.5,
          daytempfd, ask=TRUE) )
}

# If you want only the fitted functions, use plot(daytempfd)

# To plot only a single fit vs. observations, use index
# to request which one you want.  

op <- par(mfrow=c(2,1), xpd=NA, bty="n")
# xpd=NA:  clip lines to the device region,
#       not the plot or figure region
# bty="n":  Do not draw boxes around the plots.  
ylim <- range(CanadianWeather$dailyAv[,,"Temperature.C"])
# Force the two plots to have the same scale 
with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5, 
          daytempfd, index=2, titles=place, ylim=ylim, axes=FALSE) )
axis(1, monthBegin.5, labels=FALSE)
axis(1, monthEnd.5, labels=FALSE)
axis(1, monthMid, monthLetters, tick=FALSE)
axis(2)

with(CanadianWeather, plotfit.fd(dailyAv[,,"Temperature.C"], day.5, 
          daytempfd, index=35, titles=place, ylim=ylim) )
axis(1, monthBegin.5, labels=FALSE)
axis(1, monthEnd.5, labels=FALSE)
axis(1, monthMid, monthLetters, tick=FALSE)
axis(2)
par(op)

# plot residuals
with(CanadianWeather, plotfit.fd(dailyAv[, , "Temperature.C"], 
          day.5, daytempfd, residual=TRUE) )
# Can't read this, so try with 2 lines per page with ask=TRUE, 
# and limiting length(col), length(lty), etc. <=2
\dontrun{
with(CanadianWeather, plotfit.fd(
          dailyAv[,,"Temperature.C"], day.5, 
          daytempfd, residual=TRUE, col=1:2, lty=1, ask=TRUE) )
}

}
% docclass is function
\keyword{smooth}
\keyword{hplot}
