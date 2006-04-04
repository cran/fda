\name{eval.fd}
\alias{eval.fd}
\title{
  Values of a Functional Data Object
}
\description{
Evaluate a functional data object at specified argument values, or
evaluate a derivative or the result of applying a linear differential
operator to the functional object.
}
\usage{
eval.fd(evalarg, fdobj, Lfdobj=0)
}
\arguments{
\item{evalarg}{
a vector of argument values at which the functional data object is to be
evaluated.
}
\item{fdobj}{
a functional data object to be evaluated.
}
\item{Lfdobj}{
either a nonnegative integer or a linear differential operator object.
If present, the derivative or the value of applying the operator is
evaluated rather than the functions themselves.
}
}
\value{
an array of 2 or 3 dimensions containing the function
values.  The first dimension corresponds to the argument values in
\code{evalarg},
the second to replications, and the third if present to functions.
}
\seealso{
getbasismatrix, eval.bifd, eval.penalty, eval.monfd, eval.posfd
}
\examples{

#    set up the fourier basis
nbasis   <- 65
dayrange <- c(0,365)
daybasis <- create.fourier.basis(dayrange, nbasis)
#  Make temperature fd object
#  Temperature data are in 12 by 365 matrix tempav
#  See analyses of weather data.
#  Set up sampling points at mid days
daytime  <- (1:365)-0.5
#  Convert the data to a functional data object
tempfd <- data2fd(daily$tempav, daytime, daybasis)
#   set up the harmonic acceleration operator
Lbasis  <- create.constant.basis(dayrange);
Lcoef   <- matrix(c(0,(2*pi/365)^2,0),1,3)
bfdobj  <- fd(Lcoef,Lbasis)
bwtlist <- fd2list(bfdobj)
harmaccelLfd <- Lfd(3, bwtlist)
#   evaluate the value of the harmonic acceleration
#   operator at the sampling points
Ltempmat <- eval.fd(daytime, tempfd, harmaccelLfd)
#  Plot the values of this operator
matplot(daytime, Ltempmat, type="l")

}
% docclass is function
\keyword{smooth}
