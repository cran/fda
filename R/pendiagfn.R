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
