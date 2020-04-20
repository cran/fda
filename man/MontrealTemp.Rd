\name{MontrealTemp}
\alias{MontrealTemp}
\docType{data}
\title{
  Montreal Daily Temperature
}
\description{
  Temperature in degrees Celsius in Montreal each day from 1961 through
  1994
}
\usage{data(MontrealTemp)}
\format{
  A numeric array with dimnames = list(1961:1994, names(dayOfYear)).
}
%\details{}
%\source{}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer B. (2009)
  \emph{Functional Data Analysis with R and Matlab}, Springer, New York
  (esp. Section 4.3)
}
\seealso{
  \code{\link{CanadianWeather}}
  \code{\link{monthAccessories}}
}
\examples{
data(MontrealTemp)

JanuaryThaw <- t(MontrealTemp[, 16:47])

}
\keyword{datasets}
