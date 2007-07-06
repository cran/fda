smooth.basisPar <- function(argvals, y, fdobj=NULL, Lfdobj=int2Lfd(2),
      lambda=1/diff(range(argvals)), estimate=TRUE, penmat=NULL, 
      wtvec=rep(1,n), dffactor=1,
      fdnames=list(NULL, dimnames(y)[[2]], NULL) ){
##
## 1.  fdPar
##
  if(is.null(fdobj))
    fdobj <- create.bspline.basis(argvals) 
##
## 2.  fdPar
##
  fdP <- fdPar(fdobj, Lfdobj=Lfdobj, lambda=lambda,
               estimate=estimate, penmat=penmat)
##
## 3.  smooth.fd
##
  n <- length(argvals)
  w <- wtvec
#  smoothB <- smooth.basis(argvals, y, fdP, wtvec=w,
#        dffactor=dffactor, fdnames=fdnames)
  smooth.basis(argvals, y, fdP, wtvec=w,
        dffactor=dffactor, fdnames=fdnames)
}


  
