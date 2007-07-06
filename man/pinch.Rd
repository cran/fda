\name{pinch}
\alias{pinch}
\alias{pinchtime}
\title{
  pinch force data
}
\description{
  151 measurements of pinch force during 20 replications, registered, with
  time from start of measurement.   
}
\usage{
pinch
pinchtime
}
\format{

\describe{
  \item{pinch}{
    Matrix of dimension c(151, 20) = 20 measurements of pinch force
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
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York, p. 13, Figure 1.11,
  pp. 22-23, Figure 2.2, and p. 144, Figure 7.13.  
}
\examples{
  plot(pinchtime, pinch[, 1], type="b")
}
\keyword{datasets}
