\name{register.fd0}
\alias{register.fd0}
\title{
  Correct for initial position error between functional data objects.
}
\description{
  Certain physical measurements differ only in their initial position.
  Correct for this by estimating \code{x0} to minimize the following:

  integral( (yfd(x-x0)-y0fd(x))^2
  from max(x1, x1+x0) to min(x2, x2+x0)) / (x2-x1-abs(x0))
}
\usage{
register.fd0(y0fd, yfd=NULL, ...)
%\method{plot}{register.fd0}(x, ...)
}
\arguments{
  \item{y0fd}{
    a functional data object defining the target for registration.

    If \code{yfd} is NULL and y0fd is a multivariate data object, then
    y0fd is assigned to yfd and y0fd is replaced by its mean.
  }
  \item{yfd}{
    a multivariate functional data object defining the functions to be
    registered to target \code{y0fd}.  If it is NULL and \code{y0fd} is
    a multivariate functional data object, yfd takes the value of
    \code{y0fd}.
  }
  \item{\dots}{
    optional arguments passed to other methods of \code{\link{plot}}
  }
}
\details{
  Mikkelsen (2003) noted that depth positioning errors in measuring
  subsurface movement of soils can make it appear that soils are moving
  when they are not.  \code{register.fd0} differs from
  \code{register.fd} in two ways:  \code{register.fd} estimates a
  nonlinear transformation starting and ending at the same place.  By
  contrast, \code{register.fd0} shifts the starting point, using a "time
  warping" function that is linear with unit slope.  This means that the
  length of the interval of validity is shrunk by abs(x0).

  This function is designed to find the most obvious local minimum of
  the normalized integrated squared differences between
  y0df(x) and yfd(x+x0), as defined above.  This objective function may
  have multiple local minima.  \code{register.fd0} uses
  \code{\link{optimize}} to find a global minimum over x0 in c(-1,
  1)*diff(y0fd\$basis\$rangeval).

  1.  Check y0fd and yfd.

  2.  Define dy2 = integrand function and ss.dy2 = objective function.

  3.  optimize(ss.dy2, ...) for each curve in yfd.

  4.  Compute regfd = registered yfd as follows:

  4.1.  Let x00 = min(0, x0) <= 0 and x01 = max(0, x0) >= 0.  Then
  regfd\$rangeval = (rangeval[1] + x01, rangeval[2] + x00)

  4.2.  For a B-spline basis, construct a new basis retaining all the
  knots interior to regfd\$rangeval.

  4.3.  Fit each curve adjusted to this new basis to obtain regfd.

  5.  Compute dregfd = minus.fd(regfd, y0fd, regfd\$basis).
}
\value{
  a list of class \code{register.fd0} with the following components:

  \item{regfd}{
    a functional data object containing the input functional data
    objects shifted by the offset = x0, restricted to their common
    range, c(rangeval[1]+max(0, offset), rangeval[2]+min(0, offset)).
  }
  \item{dregfd}{
    the difference between regfd (yfd registered) and y0fd.

    This difference cannot currently be computed with the naive "-.fd",
    because of a missmatch between regfd$basis$rangeval and
    y0fd$basis$rangeval.  This is handled by computing minus.fd(regfd,
    y0fd, regfd$basis$rangeval).
  }
  \item{offset}{
    the vector \code{x0} containing the offset for each function in
    \code{yfd}.
  }

}

\source{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York, ch. 6 & 7.
}
\seealso{
  \code{\link{smooth.monotone}},
  \code{\link{smooth.morph}}
}
\examples{
# Specify smoothing weight
lambda.gr2.3 <- .03

# Specify what to smooth, namely the rate of change of curvature
Lfdobj.growth    <- 2

# Establish a B-spline basis
nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
rng.growth <- range(growth$age)
# 1 18
wbasis.growth <- create.bspline.basis(rangeval=rng.growth,
                   nbasis=nbasis.growth, norder=norder.growth,
                   breaks=growth$age)

# Smooth consistent with the analysis of these data
# in afda-ch06.R, and register to individual smooths:
cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth  <- fd(cvec0.growth, wbasis.growth)
growfdPar2.3 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr2.3)
# Create a functional data object for all the boys
hgtmfd.all <- with(growth, smooth.basis(age, hgtm, growfdPar2.3))

nBoys <- 2
# nBoys <- dim(growth[["hgtm"]])[2]
# register.fd takes time, so use only 2 curves as an illustration
# to minimize compute time in this example;

#Alternative to subsetting later is to subset now:
#hgtmfd.all<-with(growth,smooth.basis(age, hgtm[,1:nBoys],growfdPar2.3))

# Register the growth velocity rather than the
# growth curves directly
smBv <- deriv(hgtmfd.all$fd, 1)

smB.reg.0 <- register.fd0(smBv[1:nBoys])

op <- par(mfrow=c(3,1))
plot(smBv[1:nBoys])
lines(mean.fd(smBv[1:nBoys]), lty='dotted', col='blue')

plot(smB.reg.0$regfd)
lines(mean.fd(smBv[1:nBoys]), lty='dotted', col='blue')
# Note change of scale.
# The peak of smBv[2] at ~15.5
# is pushed off the chart as x0[2]=4.47
plot(smB.reg.0$dregfd)
# looks OK
par(op)

}
% docclass is function
\keyword{smooth}
