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
  'daily' is included primarily for compatibility with scripts written
  before the other format became available and for compatibility with
  the Matlab 'fda' code.
  \describe{
    \item{CanadianWeather}{
    A list with the following components:
    \describe{
	    \item{dailyAv}{
	      a three dimensional array c(365, 35, 3) summarizing data
	      collected at 35 different weather stations in Canada on the
	      following:

	      [,,1] = [,, 'Temperature.C']:  average daily temperature for
	      each day of the year

	      [,,2] = [,, 'Precipitation.mm']:  average daily rainfall for
	      each day of the year rounded to 0.1 mm.

	      [,,3] = [,, 'log10precip']:  base 10 logarithm of
	      Precipitation.mm after first replacing 27 zeros by 0.05 mm
	      (Ramsay and Silverman 2006, p. 248).
	    }
	    \item{place}{
	      Names of the 35 different weather stations in Canada whose data
	      are summarized in 'dailyAv'.  These names vary between 6 and 11
	      characters in length.  By contrast, daily[["place"]] which are
	      all 11 characters, with names having fewer characters being
	      extended with trailing blanks.
	    }
	    \item{province}{
	      names of the Canadian province containing each place
	    }
	    \item{coordinates}{
	      a numeric matrix giving 'N.latitude' and 'W.longitude' for each
	      place.
	    }
	    \item{region}{
	      Which of 4 climate zones contain each place:  Atlantic, Pacific,
	      Continental, Arctic.
	    }
	    \item{monthlyTemp}{
	      A matrix of dimensions (12, 35) giving the average temperature
	      in degrees celcius for each month of the year.
	    }
	    \item{monthlyPrecip}{
	      A matrix of dimensions (12, 35) giving the average daily
	      precipitation in millimeters for each month of the year.
	    }
	    \item{geogindex}{
	      Order the weather stations from East to West to North
	    }
    }
  }
  \item{daily}{
    A list with the following components:
    \describe{
	    \item{place}{
	      Names of the 35 different weather stations in Canada whose
	      data are summarized in 'dailyAv'.  These names are all 11
	      characters, with shorter names being extended with trailing
	      blanks.  This is different from CanadianWeather[["place"]],
	      where trailing blanks have been dropped.
	    }
	    \item{tempav}{
	      a matrix of dimensions (365, 35) giving the average
	      temperature in degrees celcius for each day of the year.  This
	      is essentially the same as CanadianWeather[["dailyAv"]][,,
	      "Temperature.C"].
	    }
	    \item{precipav}{
	      a matrix of dimensions (365, 35) giving the average
	      temperature in degrees celcius for each day of the year.  This
	      is essentially the same as CanadianWeather[["dailyAv"]][,,
	      "Precipitation.mm"].
	    }
    }
  }
}
}

\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{monthAccessories}}
  \code{\link{MontrealTemp}}
}
\examples{
##
## 1.  Plot (latitude & longitude) of stations by region
##
with(CanadianWeather, plot(-coordinates[, 2], coordinates[, 1], type='n',
                           xlab="West Longitude", ylab="North Latitude",
                           axes=FALSE) )
Wlon <- pretty(CanadianWeather$coordinates[, 2])
axis(1, -Wlon, Wlon)
axis(2)

rgns <- 1:4
names(rgns) <- c('Arctic', 'Atlantic', 'Continental', 'Pacific')
Rgns <- rgns[CanadianWeather$region]
with(CanadianWeather, points(-coordinates[, 2], coordinates[, 1],
                             col=Rgns, pch=Rgns) )
legend('topright', legend=names(rgns), col=rgns, pch=rgns)

##
## 2.  Plot dailyAv[, 'Temperature.C'] for 4 stations
##
data(CanadianWeather)
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

