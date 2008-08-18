\name{fda-package}
\alias{fda-package}
\alias{fda}
\docType{package}
\title{Functional Data Analysis in R}
\description{
  Functions and data sets companion to Ramsay, J. O., and Silverman,
  B. W. (2005) Functional Data Analysis, 2nd ed. and (2002) Applied
  Functional Data Analysis (Springer).  This includes finite bases
  approximations (such as splines and Fourier series) to functions fit
  to data smoothing on the integral of the squared deviations from an
  arbitrary differential operator.   
}
\details{
  \tabular{ll}{
    Package: \tab fda\cr
    Type: \tab Package\cr
    Version: \tab 2.0.0\cr
    Date: \tab 2008-05-05\cr
    License: \tab GPL-2\cr
    LazyLoad: \tab yes\cr
  }
}
\author{
  J. O. Ramsay, 

  Maintainer:  J. O. Ramsay <ramsay@psych.mcgill.ca>
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional 
    Data Analysis, 2nd ed.}, Springer, New York. 

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\examples{
##
## Simple smoothing
##
girlGrowthSm <- with(growth, smooth.basisPar(argvals=age, y=hgtf))
plot(girlGrowthSm$fd, xlab="age", ylab="height (cm)",
         main="Girls in Berkeley Growth Study" )
plot(deriv(girlGrowthSm$fd), xlab="age", ylab="growth rate (cm / year)",
         main="Girls in Berkeley Growth Study" )
plot(deriv(girlGrowthSm$fd, 2), xlab="age",
        ylab="growth acceleration (cm / year^2)",
        main="Girls in Berkeley Growth Study" )
##
## Simple basis 
##
bspl1.2 <- create.bspline.basis(norder=1, breaks=c(0,.5, 1))
plot(bspl1.2)
# 2 bases, order 1 = degree 0 = step functions:  
# (1) constant 1 between 0 and 0.5 and 0 otherwise
# (2) constant 1 between 0.5 and 1 and 0 otherwise.

fd1.2 <- Data2fd(0:1, basisobj=bspl1.2)
op <- par(mfrow=c(2,1))
plot(bspl1.2, main='bases') 
plot(fd1.2, main='fit')
par(op) 
# A step function:  0 to time=0.5, then 1 after 

}
\keyword{smooth}
