polynompen <- function(basisobj, Lfdobj=2)
{

#  Computes the polynomial penalty matrix for polynomials of the form
#      (x-ctr)^l
#  Arguments:
#  BASISOBJ ... a basis object of type "polynom"
#  LFDOBJ   ... either the order of derivative or a  nonhomogeneous 
#               linear differential operator to be penalized.
#  Returns the penalty matrix.

#  Last modified 17 January 2006

if (!(inherits(basisobj, "basis"))) stop(
    "First argument is not a basis.fd object.")

Lfdobj <- int2Lfd(Lfdobj)

type <- basisobj$type
  if (type != "polynom") stop("BASISOBJ not of type polynom")

#  Find the highest order derivative in LFD

if (inherits(Lfdobj, "Lfd")) {
    nderiv <- Lfdobj$nderiv
} else {
    stop("Second argument must be an integer or a functional data object")
}

#  Compute penalty matrix

bwtlist <- Lfdobj$bwtlist
isintLfd <- TRUE
if (nderiv > 0) {
	for (ideriv in 1:nderiv) {
		fdj <- bwtlist[[ideriv]]
		if (!is.null(fdj)) {
			if (any(fdj$coefs != 0)) {
				isintLfd <- FALSE
				break
			}
		}
	}
}

if (isintLfd) {
    nbasis   <- basisobj$nbasis
    rangex   <- basisobj$rangeval
    ctr      <- basisobj$params[1]
    basismat <- getbasismatrix(rangex, basisobj, nderiv)
    penmatl  <- outer(basismat[1,],basismat[1,])*(rangex[1] - ctr)
    penmatu  <- outer(basismat[2,],basismat[2,])*(rangex[2] - ctr)
    penaltymatrix <- matrix(0,nbasis,nbasis)
    for (i in (nderiv+1):nbasis) for (j in (nderiv+1):i) {
      	penaltymatrix[i,j] <- (penmatu[i,j] - penmatl[i,j])/(i + j - 2*nderiv - 1)
      	penaltymatrix[j,i] <- penaltymatrix[i,j]
    }
} else {
    penaltymatrix <- inprod(basisobj, basisobj, Lfdobj, Lfdobj)
}

  return( penaltymatrix )
}
