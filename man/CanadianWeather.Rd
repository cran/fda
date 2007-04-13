\name{CanadianWeather}
\alias{CanadianWeather}
\non_function{}
\title{
  Canadian average annual weather cycle
}
\description{
  Daily temperature and precipitation at 35 different locations
  in Canada averaged over 1960 to 1994.  
}
\format{
  A list with the following components:
  \describe{
    \item{dailyAv}{
      a three dimensional array c(365, 35, 3) summarizing data collected
      at 35 different weather stations in Canada on the following:

      [[1]] Temperature.C:  average daily temperature for each day of the year

      [[2]] Precipitation.mm:  average daily rainfall for each day of
      the year rounded to 0.1 mm.  

      [[3]] log10precip:  base 10 logarithm of Precipitation.mm after
      first replacing 27 zeros by 0.05 mm (Ramsay and Silverman 2006,
      p. 248).   
    }    
    \item{place}{
      Names of the 35 different weather stations in Canada whose data
      are summarized in 'dailyAv'.  
    }
    \item{monthlyTemp}{
      A matrix of dimensions (12, 35) giving the average temperature
      in degrees celcius for each month of the year.  
    }
    \item{monthlyPrecip}{
      A matrix of dimensions (12, 35) giving the average daily 
      precipitation in milimeters for each month of the year.   
    }
    \item{geogindex}{
      Order the weather stations from East to West to North 
    }
  }
}
\details{
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York. 
}
%\seealso{\code{\link{eval.penalty}}}
\examples{
# Expand the left margin to allow space for place names 
op <- par(mar=c(5, 4, 4, 5)+.1)
# Plot
stations <- c("Pr. Rupert", "Montreal", "Edmonton", "Resolute")
matplot(day.5, CanadianWeather$dailyAv[, stations, "Temperature.C"],
        type="l", axes=FALSE, xlab="", ylab="Mean Temperature (deg C)") 
axis(2, las=1)
# Label the horizontal axis with the month names
axis(1, monthBegin.5, labels=FALSE)
axis(1, monthEnd.5, labels=FALSE)
axis(1, monthMid, monthLetters, tick=FALSE)
# Add the monthly averages 
matpoints(monthMid, CanadianWeather$monthlyTemp[, stations])
# Add the names of the weather stations
mtext(stations, side=4,
      at=CanadianWeather$dailyAv[365, stations, "Temperature.C"],
     las=1)
# clean up 
par(op)
}
% docclass is function
\keyword{datasets}

