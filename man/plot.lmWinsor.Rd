\name{plot.lmWinsor}
\alias{plot.lmWinsor}
\title{
  lmWinsor plot 
}
\description{
  plot an lmWinsor model or list of models as line(s) with the data as
  points
}
\usage{
\method{plot}{lmWinsor}(x, n=101, lty=1:9, col=1:9,
         lwd=c(2:4, rep(3, 6)), lty.y=c('dotted', 'dashed'),
         lty.x = lty.y, col.y=1:9, col.x= col.y, lwd.y = c(1.2, 1),
         lwd.x=lwd.y, ...)
}
\arguments{
  \item{x}{
    an object of class 'lmWinsor', which is either a list of objects of
    class c('lmWinsor', 'lm') or is a single object of that double
    class.  Each object of class c('lmWinsor', 'lm') is the result of a
    single 'lmWinsor' fit.  If 'x' is a list, it summarizes multiple
    fits with different limits to the same data.  
  }
  \item{n}{
    integer;  with only one explanatory variable 'xNames' in the model,
    'n' is the number of values at which to evaluate the model
    predictions.  This is ignored if the number of explanatory variable
    'xNames' in the model is different from 1.  
  }
  \item{lty, col, lwd, lty.y, lty.x, col.y, col.x, lwd.y, lwd.x}{
    'lty', 'col' and 'lwd' are each replicated to a length matching the
    number of fits summarized in 'x' and used with one line for each fit
    in the order appearing in 'x'.  The others refer to horizontal and
    vertical limit lines. 
  }
  \item{\dots}{
    optional arguments for 'plot'  
  }
}
\details{
  1.  One fit or several?  
  
  2.  How many explanatory variables are involved in the model(s) in
  'x'?  If only one, then the response variable is plotted vs. that one
  explanatory variable.  Otherwise, the response is plotted
  vs. predictions. 

  3.  Plot the data.

  4.  Plot one line for each fit with its limits.  
}
\value{
invisible(NULL) 
}
\author{ Spencer Graves }
\seealso{
  \code{\link{lmWinsor}}
  \code{\link{plot}}
}
\examples{
lm.1 <- lmWinsor(y1~x1, data=anscombe)
plot(lm.1)
plot(lm.1, xlim=c(0, 15), main="other title")

# list example
lm.1. <- lmWinsor(y1~x1, data=anscombe, trim=c(0, 0.25, .4, .5)) 
plot(lm.1.)

}
\keyword{ hplot }
