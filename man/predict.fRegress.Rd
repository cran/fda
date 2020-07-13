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
    
    NOTE:  Only "interval = 'none'" has been implemented so far.  
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
###
###
### vector response with functional explanatory variable
### vector response with vector explanatory variable
###
###
##
## example from help('lm')
##
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
fRegress.D9 <- fRegress(weight ~ group)

pred.fR.D9 <- predict(fRegress.D9)
#  Now compare with 'lm'
lm.D9 <- lm(weight ~ group)
pred.lm.D9 <- predict(lm.D9)

\dontshow{stopifnot(}
all.equal(as.vector(pred.fR.D9), as.vector(pred.lm.D9))
\dontshow{)}

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

\dontshow{
\dontrun{
# The following does not yet work, because  
# the confidence interval option has not yet been coded.  
precip.Temp.p <- predict(precip.Temp.f, interval='confidence')
\dontshow{stopifnot(}
class(precip.Temp.p) == 'matrix'
\dontshow{)}

## ***** not yet implemented *****

##
## Example using se.fit
##

precip.Temp.p <- predict(precip.Temp.f, se.fit=TRUE)

\dontshow{stopifnot(}
class(precip.Temp.p) == 'list'
\dontshow{)}

## ***** not yet implemented *****

###
###
### functional response with
###            (concurrent) functional explanatory variable
###
###
##
##  predict knee angle from hip angle;  from demo('gait', package='fda')
##
(gaittime <- as.numeric(dimnames(gait)[[1]])*20)
gaitrange <- c(0,20)
gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)
harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
gaitfd <- smooth.basisPar(gaittime, gait,
       gaitbasis, Lfdobj=harmaccelLfd, lambda=1e-2)$fd
hipfd  <- gaitfd[,1]
kneefd <- gaitfd[,2]
knee.hip.f <- fRegress(kneefd ~ hipfd)

knee.hip.pred <- predict(knee.hip.f)
plot(knee.hip.pred)
}
}
}
\keyword{ models}

