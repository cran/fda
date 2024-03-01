\name{mean.fd}
\alias{mean.fd}
\title{
  Mean of Functional Data
}
\description{
Evaluate the mean of a set of functions in a functional data object.
}
\usage{
\method{mean}{fd}(x, ...)
}
\arguments{
\item{x}{
a functional data object.
}
\item{...}{Other arguments to match the generic function for 'mean'}
}
\value{
a functional data object with a single replication
that contains the mean of the functions in the object \code{fd}.
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{stddev.fd}},
  \code{\link{var.fd}},
  \code{\link{sum.fd}},
  \code{\link{center.fd}}
  \code{\link{mean}}
}
\examples{
oldpar <- par(no.readonly=TRUE)
##
## 1.  univeriate:  lip motion
##
liptime  <- seq(0,1,.02)
liprange <- c(0,1)

#  -------------  create the fd object -----------------
#       use 31 order 6 splines so we can look at acceleration

nbasis <- 51
norder <- 6
lipbasis <- create.bspline.basis(liprange, nbasis, norder)

#  ------------  apply some light smoothing to this object  -------

lipLfdobj <- int2Lfd(4)
lipLambda <- 1e-12
lipfdPar  <- fdPar(fd(matrix(0,nbasis,1),lipbasis), lipLfdobj, lipLambda)

lipfd <- smooth.basis(liptime, lip, lipfdPar)$fd
names(lipfd$fdnames) = c("Normalized time", "Replications", "mm")

lipmeanfd <- mean.fd(lipfd)
plot(lipmeanfd)

##
## 2.  Trivariate:  CanadianWeather
##
dayrng <- c(0, 365)

nbasis <- 51
norder <- 6

weatherBasis <- create.fourier.basis(dayrng, nbasis)

weather.fd <- smooth.basis(day.5, CanadianWeather$dailyAv,
                           weatherBasis)$fd

str(weather.fd.mean <- mean.fd(weather.fd))
par(oldpar)
}
\keyword{smooth}
