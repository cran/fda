fdParcheck = function(fdParobj, ncurve=NULL) {
  #  Last modified 16 November 2021 by Jim Ramsay
  if (inherits(fdParobj, "basisfd") && is.null(ncurve)) 
    stop("First argument is basisfd object and second argument is missing.")
  if (!inherits(fdParobj, "fdPar")) {
    if (inherits(fdParobj, "fd")) {
        fdParobj <- fdPar(fdParobj)
    }
    if (inherits(fdParobj, "basisfd")) {
      nbasis   <- fdParobj$nbasis
      fdParobj <- fdPar(fd(matrix(0,nbasis,ncurve),fdParobj))
    } else {
        stop(paste("'fdParobj' is not a functional parameter object,",
               "not a functional data object, and",
               "not a basis object."))
    }
  }
  return(fdParobj)
  
}
