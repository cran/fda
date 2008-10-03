\name{df.residual.fRegress}
\alias{df.residual.fRegress}
\title{
  Degress of Freedom for Residuals from a Functional Regression
}
\description{
  Effective degrees of freedom for residuals, being the trace of the
  idempotent \code{hat} matrix transforming observations into residuals
  from the fit.
}
\usage{
\method{df.residual}{fRegress}(object, ...)
}
\arguments{
  \item{ object }{
    Object of class inheriting from \code{fRegress}
  }
  \item{\dots}{
    additional arguments for other methods
  }
}
\details{
  1.  Determine N = number of observations

  2.  df.model <- object\$df

  3.  df.residual <- (N - df.model)

  4.  Add attributes
}
\value{
  The numeric value of the residual degrees-of-freedom extracted from
  \code{object} with the following attributes:

  \item{nobs}{number of observations}
  \item{df.model}{
    effective degrees of freedom for the model, being the trace of the
    idempotent linear projection operator transforming the observations
    into their predictions per the model.  This includes the intercept,
    so the 'degrees of freedom for the model' for many standard purposes
    that compare with a model with an estimated mean will be 1 less than
    this number.
  }
}
\author{ Spencer Graves }
\references{
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.
  Hastie, Trevor, Tibshirani, Robert, and Friedman, Jerome (2001)
  \emph{The Elements of Statistical Learning:  Data Mining, Inference,
  and Prediction}, Springer, New York.
}
\seealso{
  \code{\link{fRegress}}
  \code{\link[stats]{df.residual}}
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
dfr.D9 <- df.residual(fRegress.D9)

# Check the answer manually
lm.D9 <- lm(weight ~ group)
dfr.D9l <- df.residual(lm.D9)
attr(dfr.D9l, 'nobs') <- length(predict(lm.D9))
attr(dfr.D9l, 'df.model') <- 2
\dontshow{stopifnot(}
all.equal(dfr.D9, dfr.D9l)
\dontshow{)}

##
## functional response with (concurrent) functional explanatory variable
##
# *** NOT IMPLEMENTED YET FOR FUNCTIONAL RESPONSE
#     BUT WILL EVENTUALLY USE THE FOLLOWING EXAMPLE:

(gaittime <- as.numeric(dimnames(gait)[[1]])*20)
gaitrange <- c(0,20)
gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)
harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
gaitfd <- smooth.basisPar(gaittime, gait,
       gaitbasis, Lfdobj=harmaccelLfd, lambda=1e-2)$fd
hipfd  <- gaitfd[,1]
kneefd <- gaitfd[,2]

knee.hip.f <- fRegress(kneefd ~ hipfd)

#knee.hip.dfr <- df.residual(knee.hip.f)


# Check the answer
#kh.dfr <- knee.hip.f$df ???


}
\keyword{ models }

