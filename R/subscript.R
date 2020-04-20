"[.fd" <- function(fdobj, i=TRUE, j=TRUE, drop=TRUE) {
  
  #  Last moddified 16 January 2020
  
  if (!(is.fd(fdobj) || is.fdPar(fdobj)))  stop(
    "First argument is neither a functional data or a functional parameter object.")
  if (is.fdPar(fdobj)) fdobj <- fdobj$fd
  
  coef    <- as.array(fdobj$coefs)
  fdnames <- fdobj$fdnames
  coefdim <- dim(coef)
  ndim    <- length(coefdim)

  if(ndim == 2) {
    if (coefdim[2] == 1) {
      coefselect <- coef
    } else {
      coefselect <- coef[, i, drop=FALSE]
    }
    if (length(fdnames[[2]])>1) {
	    fdnames[[2]] = fdnames[[2]][i]
    }
  } else {
    if (coefdim[2] == 1) {
      coefselect <- coef
    } else {
      coefselect <- coef[, i, j,drop=drop]
    }
    if(length(fdnames[[2]])>1){
	    fdnames[[2]] = fdnames[[2]][i]
    }
    if(length(fdnames[[3]])>1){
	    fdnames[[3]] = fdnames[[3]][j]
    }
  }
  fd(coefselect, fdobj$basis, fdnames) 
}
