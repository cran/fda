\name{infantGrowth}
\alias{infantGrowth}
\docType{data}
\title{
  Tibia Length for One Baby
}
\description{
  Measurement of the length of the tibia for the first 40 days of life
  for one infant.
}
\usage{data(infantGrowth)}
\format{
  A matrix with three columns:

  \itemize{
    \item{day}{age in days}
    \item{tibiaLength}{
      The average of five measurements of tibia length in millimeters
    }
    \item{sd.length}{
      The standard deviation of five measurements of tibia length in
      millimeters
    }
  }
}
%\details{}
\source{
  Hermanussen, M., Thiel, C., von Bueren, E., de los Angeles Rol. de
  Lama, M., Perez Romero, A., Arizonavarreta Ruiz, C., Burmeister, J.,
  Tresguerras, J. A. F. (1998) Micro and macro perspectives in auxology:
  Findings and considerations upon the variability of short term and
  individual growth and the stability of population derived parameters,
  \emph{Annals of Human Biology}, 25:  359-395.
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
}
\examples{
data(infantGrowth)
plot(tibiaLength~day, infantGrowth)
}
\keyword{datasets}
