\name{pinch}
\alias{pinch}
\alias{pinchraw}
\alias{pinchtime}
\title{
  pinch force data
}
\description{
  151 measurements of pinch force during 20 replications with time from
  start of measurement.
}
\usage{
pinch
pinchraw
pinchtime
}
\format{

  \describe{
    \item{pinch, pinchraw}{
      Matrices of dimension c(151, 20) = 20 replications of measuring
      pinch force every 2 milliseconds for 300 milliseconds.  The
      original data included 300 observations.  \code{pinchraw} consists
      of the first 151 of the 300 observations.  \code{pinch} selected
      151 observations so the maximum of each curve occurred at 0.076
      seconds.
    }
    \item{pinchtime}{
      time in seconds from the start = seq(0, 0.3, 151) = every 2
      milliseconds.
    }
  }
}
\details{
  Measurements every 2 milliseconds.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\examples{
matplot (pinchtime, pinchraw, type="l", lty=1, cex=2,
         col=1, lwd=1,  xlab = "Seconds", ylab="Force (N)")
abline(h=2, lty=2)

matplot (pinchtime, pinch, type="l", lty=1, cex=2,
         col=1, lwd=1,  xlab = "Seconds", ylab="Force (N)")
abline(h=2, v=0.075, lty=2)
}
\keyword{datasets}
