\name{Lfd}
\alias{Lfd}
\title{
  Define a Linear Differential Operator Object
}
\description{
  A linear differential operator of order $m$ is defined,
  usually to specify a roughness penalty.
}
\usage{
Lfd(nderiv=0, bwtlist=vector("list", 0))
}
\arguments{
  \item{nderiv}{
    a nonnegative integer specifying the order $m$ of the
    highest order derivative in the operator
  }
  \item{bwtlist}{
    a list of length $m$.  Each member contains a
    functional data object that acts as a weight function for a
    derivative.  The first member weights the function, the
    second the first derivative, and so on up to order $m-1$.
  }
}
\value{
  a linear differential operator object
}
\details{
  To check that an object is of this class, use functions \code{is.Lfd}
  or \code{int2Lfd}.

  Linear differential operator objects are often used to define
 roughness penalties for smoothing towards a "hypersmooth" function that
 is annihilated by the operator.  For example, the harmonic acceleration
 operator used in the analysis of the Canadian daily weather data
 annihilates linear combinations of $1, sin(2 pi t/365)$ and $cos(2 pi
 t/365)$, and the larger the smoothing parameter, the closer the smooth
 function will be to a function of this shape.

 Function \code{pda.fd} estimates a linear differential operator object
 that comes as close as possible to annihilating a functional data
 object.

 A linear differential operator of order $m$ is a linear combination of
 the derivatives of a functional data object up to order $m$.  The
 derivatives of orders 0, 1, ..., $m-1$ can each be multiplied by a
 weight function $b(t)$ that may or may not vary with argument $t$.

 If the notation $D^j$ is taken to mean "take the derivative of order
 $j$", then a linear differental operator $L$ applied to function $x$
 has the expression

 $Lx(t) = b_0(t) x(t) + b_1(t)Dx(t) + ... + b_\{m-1\}(t) D^\{m-1\} x(t)
 + D^mx(t)$

 There are \code{print}, \code{summary}, and \code{plot} methods for
 objects of class \code{Lfd}.
}
\seealso{
  \code{\link{int2Lfd}},
  \code{\link{vec2Lfd}},
  \code{\link{fdPar}},
  \code{\link{pda.fd}}
  \code{\link{plot.Lfd}}
}
\examples{
#  Set up the harmonic acceleration operator
dayrange  <- c(0,365)
Lbasis  <- create.constant.basis(dayrange,
                  axes=list("axesIntervals"))
Lcoef   <- matrix(c(0,(2*pi/365)^2,0),1,3)
bfdobj  <- fd(Lcoef,Lbasis)
bwtlist <- fd2list(bfdobj)
harmaccelLfd <- Lfd(3, bwtlist)
}
% docclass is function
\keyword{smooth}
