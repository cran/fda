\name{plot.pca.fd}
\alias{plot.pca.fd}
\title{
  Plot Functional Principal Components
}
\description{
  Display the types of variation across a sample of functions.  Label
  with the eigenvalues that indicate the relative importance of each
  mode of variation.
}
\usage{
  \method{plot}{pca.fd}(x, nx = 128, pointplot = TRUE, harm = 0,
                        expand = 0, cycle = FALSE, ...)
}
\arguments{
  \item{x}{
    a functional data object.
  }
  \item{nx}{
    Number of points to plot or vector (if length > 1) to use as
    \code{evalarg} in evaluating and plotting the functional principal
    components.
  }
  \item{pointplot}{
    logical:  If TRUE, the harmonics / principal components are plotted
    as '+' and '-'.   Otherwise lines are used.
  }
  \item{harm}{
    Harmonics / principal components to plot.  If 0, plot all.

    If length(harm) > sum(par("mfrow")), the user advised, "Waiting to
    confirm page change..." and / or 'Click or hit ENTER for next page'
    for each page after the first.
  }
  \item{expand}{
    nonnegative real:  If expand == 0 then effect of +/- 2 standard
    deviations of each pc are given otherwise the factor expand is
    used.
  }
  \item{cycle}{
    logical:  If cycle=TRUE and there are 2 variables then a cycle plot
    will be drawn If the number of variables is anything else, cycle
    will be ignored.
  }
  \item{\dots}{
    other arguments for 'plot'.
  }
}
\details{
  Produces one plot for each principal component / harmonic to be
  plotted.
}
\value{
  invisible(NULL)
}
\seealso{
  \code{\link{cca.fd}},
  \code{\link{pda.fd}}
  \code{\link{plot.pca.fd}}
}
\examples{

#  carry out a PCA of temperature
#  penalize harmonic acceleration, use varimax rotation

daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)

harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
harmfdPar     <- fdPar(daybasis65, harmaccelLfd, lambda=1e5)
daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
                     daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd

daytemppcaobj <- pca.fd(daytempfd, nharm=4, harmfdPar)
#  plot harmonics, asking before each new page after the first:
plot.pca.fd(daytemppcaobj)

# plot 4 on 1 page
op <- par(mfrow=c(2,2))
plot.pca.fd(daytemppcaobj, cex.main=0.9)
par(op)

}
% docclass is function
\keyword{smooth}
