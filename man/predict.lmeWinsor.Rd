\name{predict.lmeWinsor}
\alias{predict.lmeWinsor}
\title{
  Predict method for Winsorized linear model fits with mixed effects 
}
\description{
  Model predictions for object of class 'lmeWinsor'.  
}
\usage{
\method{predict}{lmeWinsor}(object, newdata, level=Q, asList=FALSE,
      na.action=na.fail, ...)
}
\arguments{
  \item{ object }{
    Object of class inheriting from 'lmeWinsor', representing a fitted
    linear mixed-effects model.  
  }
  \item{ newdata }{
    an optional data frame to be used for obtaining the predictions.
    All variables used in the fixed and random effects models, as well
    as the grouping factors, must be present in the data frame. If
    missing, the fitted values are returned.
  }
  \item{ level }{
    an optional integer vector giving the level(s) of grouping to be
    used in obtaining the predictions. Level values increase from
    outermost to innermost grouping, with level zero corresponding to
    the population predictions. Defaults to the highest or innermost
    level of grouping. 
  }
  \item{ asList }{
    an optional logical value. If 'TRUE' and a single value is given
    in 'level', the returned object is a list with the predictions
    split by groups; else the returned value is either a vector or a
    data.frame, according to the length of 'level'. 
  }
  \item{na.action}{
    a function that indicates what should happen when 'newdata'
    contains 'NA's.  The default action ('na.fail') causes the
    function to print an error message and terminate if there are any
    incomplete observations. 
  }
  \item{\dots}{
    additional arguments for other methods
  }
}
\details{
  1.  Identify inputs and outputs as with \link{lmeWinsor}.  

  2.  If 'newdata' are provided, clip all numeric xNames to
  (object[["lower"]], object[["upper"]]). 

  3.  Call \link[nlme]{predict.lme}  

  4.  Clip the responses to the relevant components of
  (object[["lower"]], object[["upper"]]).

  5.  Done.  
}
\value{
  'predict.lmeWinsor' produces a vector of predictions or a matrix of
  predictions with limits or a list, as produced by
  \link[nlme]{predict.lme} 
}
\author{ Spencer Graves }
\seealso{
  \code{\link{lmeWinsor}}
  \code{\link[nlme]{predict.lme}}
  \code{\link{lmWinsor}}
  \code{\link[stats]{predict.lm}}
}
\examples{
requireNamespace('nlme')
fm1w <- lmeWinsor(distance ~ age, data = nlme::Orthodont,
                 random=~age|Subject)
# predict with newdata 
newDat <- data.frame(age=seq(0, 30, 2),
           Subject=factor(rep("na", 16)) )
pred1w <- predict(fm1w, newDat, level=0)

# fit with 10 percent Winsorization 
fm1w.1 <- lmeWinsor(distance ~ age, data = nlme::Orthodont,
                 random=~age|Subject, trim=0.1)
pred30 <- predict(fm1w.1)
stopifnot(all.equal(as.numeric(
              quantile(nlme::Orthodont$distance, c(.1, .9))),
          range(pred30)) )

}
\keyword{ models }

