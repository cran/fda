\name{seabird}
\alias{seabird}
\docType{data}
\title{Sea Bird Counts}
\description{
  Numbers of sightings of different species of seabirds by year 1986 -
  2005 at E. Sitkalidak, Uganik, Uyak, and W. Sitkalidak by people
  affiliated with the Kodiak National Wildlife Refuge, Alaska.
}
\usage{data(seabird)}
\format{
  A data frame with 3793 observations on the following 22 variables.
  \describe{
    \item{BAGO, BLSC, COME, COMU, CORM, HADU, HOGR, LOON, MAMU, OLDS, PIGU, RBME, RNGR, SUSC, WWSC}{
      integer count of the numbers of sightings of each species by
      transect by year
    }
    \item{Year}{integer year, 1986 - 2005}
    \item{Site}{
      Integer codes for \code{Bay}:  10 = Uyak, 20 = Uganik, 60 =
      E. Sitkalidak, 70 = W. Sitkalidak
    }
    \item{Transect}{
      Integer code (101 - 749) for the specific plot of ground observed
    }
    \item{Temp}{a numeric vector}
    \item{ObservCond}{
      a factor with levels \code{Average}, \code{Excellent},
      \code{Fair}, \code{Good}, and \code{Ideal}.  
    }
    \item{Bay}{
      a factor with levels \code{E. Sitkalidak} \code{Uganik}
      \code{Uyak} \code{W. Sitkalidak}
    }
    \item{ObservCondFactor3}{
      a factor with levels \code{ExcellentIdeal}, 
      \code{FairAverage}, and \code{Good}.  These combine levels 
      from \code{ObservCond}.
    }
  }
}
\details{
  Data provided by the Kodiak National Wildlife Refuge
}
\source{
  Zwiefelhofer, D., Reynolds, J. H., and Keim, M. (2008) Population
  trends and annual density estimates for select wintering seabird
  species on Kodiak Island, Alaska, \emph{U.S. Fish and Wildlife
    Service, Kodiak National Wildlife Refuge}, Technical Report,
  no. 08-00x
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer B. (2009)
  \emph{Functional Data Analysis with R and Matlab}, Springer, New York
  (esp. ch. 10)
}
\examples{
data(seabird)
## maybe str(seabird) ; plot(seabird) ...
}
\keyword{datasets}
