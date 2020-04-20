\name{exponentiate.fd}
\alias{exponentiate.fd}
\alias{^.fd}
\title{
  Powers of a functional data ('fd') object
}
\description{
  Exponentiate a functional data object where feasible.
}
\usage{
\method{^}{fd}(e1, e2)
exponentiate.fd(e1, e2, tolint=.Machine$double.eps^0.75,
  basisobj=e1$basis,
  tolfd=sqrt(.Machine$double.eps)*
          sqrt(sum(e1$coefs^2)+.Machine$double.eps)^abs(e2),
  maxbasis=NULL, npoints=NULL)
}
\arguments{
  \item{e1}{
    object of class 'fd'.
  }
  \item{e2}{
    a numeric vector of length 1.
  }
  \item{basisobj}{reference basis}
  \item{tolint}{
    if abs(e2-round(e2))<tolint, we assume e2 is an integer.  This
    simplifies the algorithm.
  }
  \item{tolfd}{
    the maximum error allowed in the difference between the direct
    computation \code{eval.fd(e1)^e2} and the computed representation.
  }
  \item{maxbasis}{
    The maximum number of basis functions in growing
    \code{referencebasis} to achieve a fit within \code{tolfd}.  Default
    = 2*nbasis12+1 where nbasis12 =
    nbasis of e1^floor(e2).
  }
  \item{npoints}{
    The number of points at which to compute \code{eval.fd(e1)^e2} and
    the computed representation to evaluate the adequacy of the
    representation.  Default =
    2*maxbasis-1.  For a max Fourier
    basis, this samples the highest
    frequency at all its extrema and
    zeros.
  }
}
\details{
  If e1 has a B-spline basis, this uses the B-spline algorithm.

  Otherwise it throws an error unless it finds one of the following
  special cases:

  \itemize{
    \item{e2 = 0}{
      Return an \code{fd} object with a constant basis that is
      everywhere 1
    }
    \item{e2 is a positive integer to within tolint}{
      Multiply e1 by itself e2 times
    }
    \item{e2 is positive and e1 has a Fourier basis}{
      e120 <- e1^floor(e2)

      outBasis <- e120$basis

      rng <- outBasis$rangeval

      Time <- seq(rng[1], rng[2], npoints)

      e1.2 <- predict(e1, Time)^e2

      fd1.2 <- smooth.basis(Time, e1.2, outBasis)$

      d1.2 <- (e1.2 - predict(fd1.2, Time))

      if(all(abs(d1.2)<tolfd))return(fd1.2)

      Else if(outBasis$nbasis<maxbasis) increase the size of outBasis
      and try again.

      Else write a warining with the max(abs(d1.2)) and return fd1.2.
    }
  }
}
\value{
  A function data object approximating the desired power.
}
\seealso{
  \code{\link{arithmetic.fd}}
  \code{\link{basisfd}},
  \code{\link{basisfd.product}}
}
\examples{
##
## sin^2
##

basis3 <- create.fourier.basis(nbasis=3)
plot(basis3)
# max = sqrt(2), so
# integral of the square of each basis function (from 0 to 1) is 1
integrate(function(x)sin(2*pi*x)^2, 0, 1) # = 0.5

# sin(theta)
fdsin <- fd(c(0,sqrt(0.5),0), basis3)
plot(fdsin)

fdsin2 <- fdsin^2

# check
fdsinsin <- fdsin*fdsin
# sin^2(pi*time) = 0.5*(1-cos(2*pi*theta) basic trig identity
plot(fdsinsin) # good

\dontshow{stopifnot(}
all.equal(fdsin2, fdsinsin)
\dontshow{)}

% ##
% ## sqrt(sin2)
% ##
% plot(fdsin2)
% fdsin. <- sqrt(fdsin2)
% plot(fdsin, main='fdsin and sqrt(fdsin^2)')
% lines(fdsin., col='red')
% # fdsin is positive and negative
% # fdsin. = sqrt(fdsin^2)
% # has trouble, because it wants to be smooth
% # but theoretically has a discontinuous first derivative at 0.5

% fdsin.5.2 <- fdsin.^2
% resin <- fdsin2-fdsin.5.2
% plot(resin)
% \dontshow{stopifnot(}
% max(abs(resin$coefs))<0.01
% \dontshow{)}

% ##
% ## x^2, x = straight line f(x)=x
% ##
% bspl1 <- create.bspline.basis(norder=2)
% x <- fd(0:1, bspl1)
% plot(x)

% x2 <- x^2
% plot(x2)

% er.x <- x-sqrt(x2)
% er.x$coefs
% max(er.x$coefs)
% \dontshow{stopifnot(}
% max(abs(er.x$coefs))<10*.Machine$double.eps
% \dontshow{)}

}
% docclass is function
\keyword{smooth}

