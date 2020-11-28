\name{predict.fRegress}
\alias{predict.fRegress}
\title{
  Predict method for Functional Regression
}
\description{
  Model predictions for object of class \code{fRegress}.
}
\usage{
\method{predict}{fRegress}(object, newdata=NULL, se.fit = FALSE,
     interval = c("none", "confidence", "prediction"),
     level = 0.95, ...)
}
\arguments{
  \item{object}{
    Object of class inheriting from \code{fRegress}
  }
  \item{newdata}{
    Either NULL or a list matching object\$xfdlist.

    If(is.null(newdata)) predictions <- object\$yhatfdobj

    If newdata is a list, predictions = the sum of either 
    
    newdata[i] * betaestfdlist[i] if object\$yfdobj has class \code{fd} 
    
    or 
    
    inprod(newdata[i], betaestfdlist[i]) if class(object\$yfdobj) =
    \code{numeric}.
  }
  \item{se.fit}{
    a switch indicating if standard errors of predictions are required
    
    NOTE:  se.fit = TRUE is NOT IMPLEMENTED YET.  
  }
  \item{interval}{
    type of prediction (response or model term)
    
    NOTE:  Only "intervale = 'none'" has been implemented so far.  
  }
  \item{level}{Tolerance/confidence level}
  \item{\dots}{
    additional arguments for other methods
  }
}
\details{
  1.  Without \code{newdata}, fit <- object\$yhatfdobj.

  2.  With \code{newdata}, if(class(object\$y) == 'numeric'), fit <- sum
  over i of inprod(betaestlist[i], newdata[i]).  Otherwise, fit <- sum
  over i of betaestlist[i] * newdata[i].

  3.  If(se.fit | (interval != 'none')) compute \code{se.fit}, then
  return whatever is desired.
}
\value{
  The predictions produced by \code{predict.fRegress} are either a
  vector or a functional parameter (class \code{fdPar}) object, matching
  the class of \code{object\$y}.

  If \code{interval} is not "none", the predictions will be
  multivariate for \code{object\$y} and the requested \code{lwr} and
  \code{upr} bounds.  If \code{object\$y} is a scalar, these predictions
  are returned as a matrix;  otherwise, they are a multivariate
  functional parameter object (class \code{fdPar}).

  If \code{se.fit} is \code{TRUE}, \code{predict.fRegress} returns a
  list with the following components:

  \item{fit}{
    vector or matrix or univariate or multivariate functional parameter
    object depending on the value of \code{interval} and the class of
    \code{object\$y}.
  }
  \item{se.fit}{standard error of predicted means}
}
\author{ Spencer Graves }
\seealso{
  \code{\link{fRegress}}
  \code{\link[stats]{predict}}
}
\examples{

##
## vector response with functional explanatory variable
##

annualprec <- log10(apply(CanadianWeather$dailyAv[,,
                   "Precipitation.mm"], 2,sum))
smallbasis  <- create.fourier.basis(c(0, 365), 25)
tempfd <- smooth.basis(day.5,
   CanadianWeather$dailyAv[,,"Temperature.C"], smallbasis)$fd
precip.Temp.f <- fRegress(annualprec ~ tempfd)

precip.Temp.p <- predict(precip.Temp.f)

# plot response vs. fitted 
plot(annualprec, precip.Temp.p)

}
\keyword{ models}

