\name{center.fd}
\alias{center.fd}
\title{
  Center Functional Data
}
\description{
Subtract the pointwise mean from each of the functions
in a functional data object; that is, to center them on the mean function.
}
\usage{
center.fd(fdobj)
}
\arguments{
\item{fdobj}{
a functional data object to be centered.
}
}
\value{
a functional data object whose mean is zero.
}
\seealso{
\code{\link{mean.fd}},
\code{\link{sum.fd}},
\code{\link{stddev.fd}},
\code{\link{std.fd}}
}
\examples{
daytime    <- (1:365)-0.5
daybasis   <- create.fourier.basis(c(0,365), 365)
harmLcoef  <- c(0,(2*pi/365)^2,0)
harmLfd    <- vec2Lfd(harmLcoef, c(0,365))
templambda <- 0.01
tempfdPar  <- fdPar(daybasis, harmLfd, templambda)

# do not run on CRAN because it takes too long.
if(!CRAN()){
tempfd     <- smooth.basis(daytime,
       CanadianWeather$dailyAv[,,"Temperature.C"], tempfdPar)$fd

tempctrfd  <- center.fd(tempfd)

plot(tempctrfd, xlab="Day", ylab="deg. C",
     main = "Centered temperature curves")
}
}
% docclass is function
\keyword{smooth}
