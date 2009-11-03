\name{residuals.fRegress}
\alias{residuals.fRegress}
\title{
  Residuals from a Functional Regression
}
\description{
  Residuals from a model object of class \code{fRegress}.
}
\usage{
\method{residuals}{fRegress}(object, ...)
}
\arguments{
  \item{ object }{
    Object of class inheriting from \code{lmWinsor}
  }
  \item{\dots}{
    additional arguments for other methods
  }
}
\details{
  object\$y - predict(object)
}
\value{
  The residuals produced by \code{resid.fRegress} or
  \code{residuals.fRegress} are either a vector or a functional
  parameter (class \code{fdPar}) object, matching
  the class of \code{object\$y}.
}
\author{ Spencer Graves }
\seealso{
  \code{\link{fRegress}}
  \code{\link{predict.fRegress}}
  \code{\link[stats]{residuals}}
}
\examples{
##
## example from help('lm')
##
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2,10,20, labels=c("Ctl","Trt"))
weight <- c(ctl, trt)
fRegress.D9 <- fRegress(weight ~ group)

resid.fR.D9 <- resid(fRegress.D9)
#  Now compare with 'lm'
lm.D9 <- lm(weight ~ group)
resid.lm.D9 <- resid(lm.D9)

\dontshow{stopifnot(}
all.equal(as.vector(resid.fR.D9), as.vector(resid.lm.D9))
\dontshow{)}

##
##  resid from knee angle prediciton from hip angle;
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

#knee.hip.e <- resid(knee.hip.f)
#plot(knee.hip.e)



}
\keyword{ models }

