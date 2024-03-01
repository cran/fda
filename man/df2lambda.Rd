\name{df2lambda}
\alias{df2lambda}
\title{
  Convert Degrees of Freedom to a Smoothing Parameter Value
}
\description{
The degree of roughness of an estimated function is controlled by a
smoothing parameter $lambda$ that directly multiplies the penalty.
However, it can be difficult to interpret or choose this value, and it
is often easier to determine the roughness by choosing a value that is
equivalent of the degrees of freedom used by the smoothing procedure.
This function converts a degrees of freedom value into a multiplier
$lambda$.
}
\usage{
df2lambda(argvals, basisobj, wtvec=rep(1, n), Lfdobj=0,
          df=nbasis)
}
\arguments{
\item{argvals}{
a vector containing argument values associated with the values to
be smoothed.
}
\item{basisobj}{
a basis function object.
}
\item{wtvec}{
a vector of weights for the data to be smoothed.
}
\item{Lfdobj}{
either a nonnegative integer or a linear differential operator object.
}
\item{df}{
the degrees of freedom to be converted.
}
}
\value{
a positive smoothing parameter value $lambda$
}
\details{
The conversion requires a one-dimensional optimization and may be
therefore computationally intensive.
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
  \code{\link{lambda2df}}, 
  \code{\link{lambda2gcv}}
}
\examples{

#  Smooth growth curves using a specified value of
#  degrees of freedom.
#  Set up the ages of height measurements for Berkeley data
age <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))
#  Range of observations
rng <- c(1,18)
#  Set up a B-spline basis of order 6 with knots at ages
knots  <- age
norder <- 6
nbasis <- length(knots) + norder - 2
hgtbasis <- create.bspline.basis(rng, nbasis, norder, knots)
#  Find the smoothing parameter equivalent to 12
#  degrees of freedom
lambda <- df2lambda(age, hgtbasis, df=12)
#  Set up a functional parameter object for estimating
#  growth curves.  The 4th derivative is penalyzed to
#  ensure a smooth 2nd derivative or acceleration.
Lfdobj <- 4
growfdPar <- fdPar(fd(matrix(0,nbasis,1),hgtbasis), Lfdobj, lambda)
#  Smooth the data.  The data for the girls are in matrix
#  hgtf.
hgtffd <- smooth.basis(age, growth$hgtf, growfdPar)$fd
#  Plot the curves
oldpar <- par(no.readonly=TRUE)
plot(hgtffd)
par(oldpar)
}
% docclass is function
\keyword{smooth}
