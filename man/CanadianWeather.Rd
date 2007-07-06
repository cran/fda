\name{CanadianWeather}
\docType{data}
\alias{CanadianWeather}
\alias{daily}
\title{
  Canadian average annual weather cycle
}
\description{
  Daily temperature and precipitation at 35 different locations
  in Canada averaged over 1960 to 1994.  
}
\usage{
CanadianWeather
daily
}
\format{
  'CanadianWeather' and 'daily' are lists containing essentially the
  same data.  'CanadianWeather' may be preferred for most purposes;
  'daily' is included primarily for compatability with scripts written
  before the other format became available and for compatability with
  the Matlab 'fda' code.    

  \item{CanadianWeather}{
    A list with the following components:
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
      are summarized in 'dailyAv'.  These names vary between 6 and 11
      characters in length.  By contrast, daily[["place"]] which are all
      11 characters, with names having fewer characters being extended
      with trailing blanks.  
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

  \item{daily}{
    A list with the following components:
    \item{place}{
      Names of the 35 different weather stations in Canada whose data
      are summarized in 'dailyAv'.  These names are all 11 characters,
      with shorter names being extended with trailing blanks.  This is
      different from CanadianWeather[["place"]], where trailing blanks
      have been dropped.    
    }
    \item{tempav}{
      a matrix of dimensions (365, 35) giving the average temperature in
      degrees celcius for each day of the year.  This is essentially the
      same as CanadianWeather[["dailyAv"]][,, "Temperature.C"].  
    }
    \item{precipav}{
      a matrix of dimensions (365, 35) giving the average temperature in
      degrees celcius for each day of the year.  This is essentially the
      same as CanadianWeather[["dailyAv"]][,, "Precipitation.mm"].  
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

