\name{predict.lmWinsor}
\alias{predict.lmWinsor}
\title{
  Predict method for Winsorized linear model fits 
}
\description{
  Model predictions for object of class 'lmWinsor'.  
}
\usage{
\method{predict}{lmWinsor}(object, newdata, se.fit = FALSE,
     scale = NULL, df = Inf, 
     interval = c("none", "confidence", "prediction"),
     level = 0.95, type = c("response", "terms"), 
     terms = NULL, na.action = na.pass,
     pred.var = res.var/weights, weights = 1, ...)
}
\arguments{
  \item{ object }{
    Object of class inheriting from 'lmWinsor'  
  }
  \item{ newdata }{
    An optional data frame in which to look for variables with which to
    predict.  If omitted, the fitted values are used. 
  }
  \item{ se.fit}{
    a switch indicating if standard errors of predictions are required 
  }
  \item{ scale }{
    Scale parameter for std.err. calculation  
  }
  \item{df}{
    degrees of freedom for scale 
  }
  \item{interval}{
    type of prediction (response or model term) 
  }
  \item{ level }{ Tolerance/confidence level }
  \item{ type }{
    Type of prediction (response or model term);  see
    \link[stats]{predict.lm} 
  }
  \item{terms}{
    If 'type="terms"', which terms (default is all terms)
  }
  \item{na.action}{
    function determining what should be done with missing values in
    'newdata'.  The default is to predict 'NA'. 
  }
  \item{ pred.var }{
    the variance(s) for future observations to be assumed for prediction
    intervals.  See \link[stats]{predict.lm} 'Details'.
  }
  \item{ weights }{
    variance weights for prediction. This can be a numeric vector or a
    one-sided model formula. In the latter case, it is interpreted as an
    expression evaluated in 'newdata'
  }
  \item{\dots}{
    additional arguments for other methods
  }
}
\details{
  1.  Identify inputs and outputs via mdly <- mdlx <- formula(object);
  mdly[[3]] <- NULL;  mdlx[[2]] <- NULL;  xNames <- all.vars(mdlx);
  yNames <- all.vars(mdly).  Give an error if as.character(mdly[[2]]) !=
  yNames.  

  2.  If 'newdata' are provided, clip all xNames to (object[["lower"]],
  object[["upper"]]). 

  3.  Call \link[stats]{predict.lm}  

  4.  Clip the responses to the relevant components of
  (object[["lower"]], object[["upper"]]).

  5.  Done.  
}
\value{
  'predict.lmWinsor' produces a vector of predictions or a matrix of
  predictions with limits or a list, as produced by
  \link[stats]{predict.lm} 
}
\author{ Spencer Graves }
\seealso{
  \code{\link{lmWinsor}}
  \code{\link[stats]{predict.lm}}
}
\examples{
# example from 'anscombe' 
lm.1 <- lmWinsor(y1~x1, data=anscombe)

newD <- data.frame(x1=seq(1, 22, .1))
predW <- predict(lm.1, newdata=newD) 
plot(y1~x1, anscombe, xlim=c(1, 22)) 
lines(newD[["x1"]], predW, col='blue')
abline(h=lm.1[['lower']]['y1'], col='red', lty='dashed') 
abline(h=lm.1[['upper']]['y1'], col='red', lty='dashed')
abline(v=lm.1[['lower']]['x1'], col='green', lty='dashed') 
abline(v=lm.1[['upper']]['x1'], col='green', lty='dashed') 

}
\keyword{ models }

