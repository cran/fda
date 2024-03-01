\name{refinery}
\alias{refinery}
\title{
  Reflux and tray level in a refinery
}
\description{
  194 observations on reflux and "tray 47 level" in a distallation
  column in an oil refinery.  
}
\format{
  A data.frame with the following components:  
  \describe{
    \item{Time}{
      observation time 0:193 
    }
    \item{Reflux}{
      reflux flow centered on the mean of the first 60 observations 
    }
    \item{Tray47}{
      tray 47 level centered on the mean of the first 60
      observations 
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
%\seealso{\code{\link{eval.penalty}}}
\examples{
    attach(refinery)
# allow space for an axis on the right 
    op <- par(mar=c(5, 4, 4, 5)+0.1)
# plot uval 
    plot(Time, Reflux, type="l", bty="n")
# add yval 
    y.u <- diff(range(Tray47))/diff(range(Reflux))
    u0 <- min(Reflux)
    y0 <- min(Tray47)

    lines(Time, u0+(Tray47-y0)/y.u, lty=3, lwd=1.5, col="red")
    y.tick <- pretty(range(Tray47))
    axis(4, at=u0+(y.tick)/y.u, labels=y.tick, col="red", lty=3,
            lwd=1.5)
# restore previous plot margins
    par(op)
    detach(refinery)
}
% docclass is function
\keyword{datasets}

