"[.fd" <- function(fdobj, i=TRUE, j=TRUE, drop=TRUE) {
  #  select subsets of curves in a functional data object

  coef    <- fdobj$coefs
  fdnames <- fdobj$fdnames
  ndim    <- length(dim(coef))

  if(ndim == 2) {
    coefselect <- coef[, i,drop=FALSE]
    if(length(fdnames[[2]])>1){
	    fdnames[[2]] = fdnames[[2]][i]
    }
  } else {
    coefselect <- coef[, i, j,drop=drop]
    if(length(fdnames[[2]])>1){
	    fdnames[[2]] = fdnames[[2]][i]
    }
    if(length(fdnames[[3]])>1){
	    fdnames[[3]] = fdnames[[3]][j]
    }
  }
  fd(coefselect, fdobj$basis, fdnames) 
}
