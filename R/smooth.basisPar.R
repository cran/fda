smooth.basisPar <- function(argvals, y, fdobj=NULL, Lfdobj=NULL,
      lambda=0, estimate=TRUE, penmat=NULL,
      wtvec=NULL, fdnames=NULL ){
##
## 1.  fdobj
##
  {
    if(is.null(fdobj))
      fdobj <- create.bspline.basis(argvals)
    else {
      if(is.numeric(fdobj)){
        if(length(fdobj)==1) {
          if(round(fdobj) != fdobj)
            stop("'fdobj' is numeric but not an integer")
          fdobj <- create.bspline.basis(argvals, norder=fdobj)
        }
        else
          fdobj <- fd(fdobj)
      }
    }
  }
##
## 2.  fdPar
##
  fdP <- fdPar(fdobj, Lfdobj=Lfdobj, lambda=lambda,
               estimate=estimate, penmat=penmat)
##
## 3.  smooth.fd
##
#  n <- length(argvals)
  w <- wtvec
#  smoothB <- smooth.basis(argvals, y, fdP, wtvec=w,
#        dffactor=dffactor, fdnames=fdnames)
  smooth.basis(argvals, y, fdP, wtvec=w, fdnames=fdnames)
}
