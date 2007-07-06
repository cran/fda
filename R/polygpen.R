polygpen <- function(basisobj, Lfdobj=int2Lfd(1))
{

#  Computes the polygonal penalty matrix.
#  Arguments:
#  BASISOBJ ... a basis object of type "polyg"
#  LFDOBJ   ... either the order of derivative or a
#               linear differential operator to be penalized.
#          The highest derivative must be either 0 or 1.
#  Returns the penalty matrix.

#  Last modified 2007.11.28 by Spencer Graves
#  previously modified 26 October 2005

if (!(inherits(basisobj, "basisfd"))) stop(
    "First argument is not a basis object.")

Lfdobj <- int2Lfd(Lfdobj)

type <- basisobj$basis
if (type != "polyg") stop("BASISOBJ not of type polyg")

nderiv <- Lfdobj$nderiv

if (nderiv > 1) stop(
    "Derivative greater than 1 cannot be taken for polygonal basis.")

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
    args    <- basisobj$params
    n       <- length(args)
    argdiff <- diff(args)
    penaltymatrix <- diag(rep(1,n))
    if (nderiv == 0) {
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
      	diag(penaltymatrix[indx,  indx  ]) <- argdiff[indx]+argdiff[indx-1]
      	indx <- 2:n
      	diag(penaltymatrix[indx  ,indx-1]) <- -argdiff
      	diag(penaltymatrix[indx-1,indx  ]) <- -argdiff
    }
} else {
    penaltymatrix <- inprod(basisobj, basisobj, Lfdobj, Lfdobj)
}

return( penaltymatrix )
}
