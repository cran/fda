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
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
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
