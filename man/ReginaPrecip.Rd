\name{ReginaPrecip}
\alias{ReginaPrecip}
\docType{data}
\title{
  Regina Daily Precipitation
}
\description{
  Temperature in millimeters in June in Regina, Saskatchewan, Canada,
  1960 - 1993, omitting 16 missing values.
}
\usage{data(ReginaPrecip)}
\format{
  A numeric vector of length 1006.
}
%\details{}
%\source{}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer B. (2009)
  \emph{Functional Data Analysis with R and Matlab}, Springer, New York
  (esp. Section 5.4.3)
}
\seealso{
  \code{\link{CanadianWeather}}
  \code{\link{MontrealTemp}}
}
\examples{
data(ReginaPrecip)
hist(ReginaPrecip)

}
\keyword{datasets}
