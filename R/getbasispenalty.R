getbasispenalty <- function(basisfd, Lfd=NULL)
{
#  Computes the penalty matrix  associated with basis.fd object BASISFD.
#    This is defined in terms of a linear differential operator Lfd
#    The default for Lfd depends on the nature of the basis.

#  Last modified 13 December 2002

if (!(inherits(basisfd, "basis.fd"))) stop(
    "First argument is not a basis object.")

type   <- getbasistype(basisfd)
nbasis <- basisfd$nbasis

if        (type == "fourier") {
    if (is.null(Lfd)) Lfd <- 2
    penalty <- fourierpen(basisfd, Lfd)
} else if (type == "bspline") {
    norder <- basisfd$nbasis - length( basisfd$params )
    if (is.null(Lfd)) Lfd <- as.integer(norder/2)
    penalty <- bsplinepen(basisfd, Lfd)
} else if (type == "poly")    {
    if (is.null(Lfd)) Lfd <- 2
    penalty <- polynompen(basisfd, Lfd)
} else if (type == "expon")   {
    if (is.null(Lfd)) Lfd <- 2
    penalty <- exponpen(basisfd, Lfd)
} else if (type == "polyg")   {
    if (is.null(Lfd)) Lfd <- 1
    penalty <- polygpen(basisfd, Lfd)
} else if (type == "power")   {
    if (is.null(Lfd)) Lfd <- 2
    penalty <- powerpen(basisfd, Lfd)
} else if (type == "const")   {
    if (is.null(Lfd)) Lfd <- 0
    if (Lfd == 0) {
      penalty <- basisfd$rangeval[2] - basisfd$rangeval[1]
    } else {
      penalty <- 0
    }
} else {
    stop("Basis type not recognizable")
}

return(penalty)
}
