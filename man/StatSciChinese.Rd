\name{StatSciChinese}
\alias{StatSciChinese}
\docType{data}
\title{Statistical Science in Chinese}
\description{
  (x, y, z) coordinates of the location of the tip of a pen during fifty
  replications of writing 'Statistical Science' in simplified Chinese at
  10 millisecond intervals
}
\usage{data(StatSciChinese)}
\format{
  a 3-dimensional array of dimensions (601, 50, 3) containing 601
  observations of (x, y, z) coordinates of the tip of a pen at 2.5
  millisecond intervals for each of 50 repitions of writing 'Statistical
  Science' in simplified Chinese.
}
\details{
  Xiaochun Li wrote 'Statistical Science' in simplified Chinese 50
  times.  An infra-red detecting tablet was attached to the tip of the
  pen, and a wall-mounted set of three cameras recorded its position 400
  times per second with an error level of about 0.5 millimeters.  Each
  sample required about 6 seconds to produce, and for simplicity, time
  was normalized to this interval for all 50 records.  The script
  requires 50 strokes, with an average time of 120 milliseconds per
  stroke.  These raw data were shifted and rotated so the numbers more
  accruately reflected \code{x} and \code{y} coordinates relative to the
  drawn characters plus vertical distance from the paper.
}
%\source{}
\references{
  Ramsay, James O. (2000) Functional Components of Variation in
  Handwriting, \emph{Journal of the American Statistical Association},
  95, 9-15.
}
\examples{
data(StatSciChinese)

i <- 3
StatSci1 <- StatSciChinese[, i, ]
# Where does the pen leave the paper?
plot(StatSci1[, 3], type='l')
thresh <- quantile(StatSci1[, 3], .8)
abline(h=thresh)

sel1 <- (StatSci1[, 3] < thresh)
StatSci1[!sel1, 1:2] <- NA
plot(StatSci1[, 1:2], type='l')

mark <- seq(1, 601, 12)
points(StatSci1[mark, 1], StatSci1[mark, 2])
}
\keyword{datasets}
