\name{ElectricDemand}
\alias{ElectricDemand}
\title{
  Predicting electricity demand in Adelaide from temperature
}
\description{
  The data sets used in this demonstration analysis consist of half-hourly electricity demands from Sunday to Saturday in Adelaide, Australia, between July 6, 1976 and March 31, 2007.  Also provided in the same format and times are the temperatures in degrees Celsius at the Adelaide airport.  The shapes of the demand curves for each day resemble those for temperature, and the goal is to see how well a concurrent functional regress model fit by function \code{fRegress} can fit the demand curves.
}
\format{
   There is a data object for each day of the week, and in each object, the member y, such as \code{mondaydemand$y}, is a 48 by 508 matrix of electricity demand values, one for each of 48 half-hourly points, and for each of 508 weeks.   
}
\details{
  The demonstration is designed to show how to use the main functional regression function \code{fRegress} in order to fit a functional dependent variable (electricity demand) from an obviously important input or covariate functional variable (temperature).  In the texts cited, this model is referred to as a concurrent functional regression because the model assumes that changes in the input functional variable directly and immediately change the dependent functional variable.  
  
  Also illustrated is the estimation of pointwise confidence intervals for the regression coefficient functions and for the fitting functions that approximate electricity demand.  These steps involve functons \code{fRegress.stderr} and \code{predict.fRegress}.
  
  The data are supplied through the CRAN system by package \code{fds}, and further information as well as many other datasets are to be found there.
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York. 
    
  Ramsay, James O., Hooker, Giles and Graves, Spencer (2009)  Functional
    Data Analysis with R and Matlab, Springer, New York.
}
\keyword{datasets}
