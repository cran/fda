fourier <- function(x, nbasis = n, period = span, nderiv = 0)
{
  #  Computes the NDERIV derivative of the Fourier series basis
  #    for NBASIS functions with period PERIOD, these being evaluated
  #    at values in vector X
  #  Returns an N by NBASIS matrix of function values
  #  Note:  The number of basis functions always odd.  If the argument
  #   NBASIS is even, it is increased by one.

  #  last modified 8 June 99

  x      <- as.vector(x)
  n      <- length(x)
  onen   <- rep(1,n)
  xrange <- range(x)
  span   <- xrange[2] - xrange[1]
  if (nbasis <= 0) stop('NBASIS not positive')
  if (period <= 0) stop('PERIOD not positive')
  if (nderiv <  0) stop('NDERIV negative')

  if (2*(nbasis %/% 2) == nbasis) nbasis <- nbasis + 1
  basis  <- matrix(0,n,nbasis)
  omega  <- 2*pi/period
  omegax <- omega*x

  if (nderiv == 0) {
    #  The fourier series itself is required.
    basis[,1] <- 0.7071068
    j    <- seq(2,nbasis-1,2)
    k    <- j/2
    args <- outer(omegax,k)
    basis[,j]   <- sin(args)
    basis[,j+1] <- cos(args)
  } else {
    #  A derivative of the fourier series is required.
    basis[,1] <- 0.0
    if (nderiv == floor(nderiv/2)*2) {
      mval  <- nderiv/2
      ncase <- 1
    } else {
      mval <- (nderiv-1)/2
      ncase <- 2
    }
    j    <- seq(2,nbasis-1,2)
    k    <- j/2
    fac  <- outer(onen,((-1)^mval)*(k*omega)^nderiv)
    args <- outer(omegax,k)
    if (ncase == 1) {
      basis[,j]   <-  fac * sin(args)
      basis[,j+1] <-  fac * cos(args)
    } else {
      basis[,j]   <-  fac * cos(args)
      basis[,j+1] <- -fac * sin(args)
    }
  }
  basis <- basis/sqrt(period/2)
  return(basis)
}
