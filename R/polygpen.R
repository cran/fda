polygpen <- function(basisfd, Lfd=1)
{

#  Computes the polygonal penalty matrix.
#  Arguments:
#  BASISFD   ... a basis.fd object of type "bspline"
#  LFD ... either the order of derivative or a
#          linear differential operator to be penalized.
#          The highest derivative must be either 0 or 1.
#  Returns the penalty matrix.

#  Last modified 5 December 2001

if (!(inherits(basisfd, "basis.fd"))) stop(
    "First argument is not a basis.fd object.")

type <- getbasistype(basisfd)
if (type != "polyg") stop("BASISFD not of type polyg")

#  Find the highest order derivative in LFD

if (is.numeric(Lfd)) {
    if (length(Lfd) == 1) {
      	nderiv <- Lfd
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
} else if (inherits(Lfd, "fd")) {
   	derivcoef <- getcoef(Lfd)
   	if (length(dim(derivcoef))==3) derivcoef <- derivcoef[,,1]
   	nderiv <- dim(derivcoef)[2] - 1
   	if (nderiv < 0) {
   		stop("Order of derivative must be 0 or positive")
   	}
    nderiv <- ncol(derivcoef)
} else {
    stop("Second argument must be an integer or a functional data object")
}

if (nderiv > 1) stop(
    "Derivative greater than 1 cannot be taken for polygonal basis.")

if (is.numeric(Lfd)) {
    args    <- basisfd$params
    n       <- length(args)
    argdiff <- diff(args)
    penaltymatrix <- diag(rep(1,n))
    if (Lfd == 0) {
      	penaltymatrix[1,1] <- argdiff[  1]/3
      	penaltymatrix[n,n] <- argdiff[n-1]/3
      	indx <- 2:(n-1)
      	diag(penaltymatrix[indx  ,indx  ]) <- (argdiff[indx]+argdiff[indx-1])/3
      	indx <- 2:n
      	diag(penaltymatrix[indx  ,indx-1]) <- argdiff/6
      	diag(penaltymatrix[indx-1,indx  ]) <- argdiff/6
    } else {
      	argdiff <- 1/argdiff
      	penaltymatrix[1,1] <- argdiff[  1]
      	penaltymatrix[n,n] <- argdiff[n-1]
      	indx <- 2:(n-1)
      	diag(penaltymatrix[indx,  indx  ]) <- argdiff[ind]+argdiff[ind-1]
      	indx <- 2:n
      	diag(penaltymatrix[indx  ,indx-1]) <- -argdiff
      	diag(penaltymatrix[indx-1,indx  ]) <- -argdiff
    }
} else {
    penaltymatrix <- inprod(basisfd, basisfd, Lfd, Lfd)
}

return( penaltymatrix )
}
