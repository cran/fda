fourierpen <- function(basisobj, Lfdobj=2)
{

  #  Computes the Fourier penalty matrix.
  #  Arguments:
  #  BASISOBJ ... a basis object of type "fourier"
  #  LFDOBJ   ... either the order of derivative or a
  #                linear differential operator to be penalized.
  #  Returns  the penalty matrix.

  #  Note:  The number of basis functions is always odd.  If BASISOBJ
  #  specifies an even number of basis functions, then the number of basis
  #  functions is increased by 1, and this function returns a matrix of
  #  order one larger.

  #  Last modified 21 January 2003

  if (!(inherits(basisobj, "basis.fd"))) stop(
    "First argument is not a basis.fd object.")

  nbasis <- basisobj$nbasis
  if (2*(nbasis %/% 2) == nbasis) basisobj$nbasis <- nbasis + 1

  type <- getbasistype(basisobj)
  if (type != "fourier") stop ("Wrong basis type")

  if (is.numeric(Lfdobj)) {
    if (length(Lfdobj) == 1) {
      nderiv <- Lfdobj
      if (nderiv != as.integer(nderiv)) {
        stop("Order of derivative must be an integer")
      }
      if (nderiv < 0) {
        stop("Order of derivative must be 0 or positive")
      }
    } else {
      stop("Order of derivative must be a single number")
    }
    if (nderiv < 0) stop ("Order of derivative cannot be negative")
  } else if (inherits(Lfdobj, "fd")) {
    derivcoef <- getcoef(Lfdobj)
    nderiv    <- ncol(derivcoef)
    Lfdbasis  <- getbasis(Lfdobj)
  } else {
    stop("Second argument must be an integer or a functional data object")
  }

  if (nderiv < 0) stop("NDERIV is negative")

  width  <- basisobj$rangeval[2] - basisobj$rangeval[1]
  period <- basisobj$params[1]

  if (period == width &&
      (is.numeric(Lfdobj) || getbasistype(Lfdbasis) == "const")) {

    #  Compute penalty matrix for penalizing integral over one period.

    pendiag <- pendiagfn(basisobj, nderiv)

    if (!is.numeric(Lfdobj)) {
      for (i in 1:nderiv) {
        if (derivcoef[1,i] != 0)
          pendiag <- pendiag + derivcoef[1,i]*pendiagfn(basisobj,i-1)
      }
    }

    penaltymatrix <- diag(pendiag)

  } else {

    #  Compute penalty matrix by numerical integration

    penaltymatrix <- inprod(basisobj, basisobj, Lfdobj, Lfdobj)

  }

  return( penaltymatrix )
}

#  ------------------------------------------------------------------

pendiagfn <- function(basisfd, nderiv) {

    nbasis  <- basisfd$nbasis
    period  <- basisfd$params[1]
    rangev  <- basisfd$rangeval
    omega   <- 2*pi/period
    halfper <- period/2
    twonde  <- 2*nderiv
    pendiag <- rep(0,nbasis)
    if (nderiv == 0) pendiag[1] <- period/2.0 else pendiag[1] <- 0
    j   <- seq(2,nbasis-1,2)
    fac <- halfper*(j*omega/2)^twonde
    pendiag[j]   <- fac
    pendiag[j+1] <- fac
    pendiag <- 2*pendiag/period
    return(pendiag)
}
