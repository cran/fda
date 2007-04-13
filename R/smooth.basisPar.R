smooth.basisPar <- function(argvals, y, fdobj=fd(),
      Lfdobj=int2Lfd(0), lambda=0, estimate=TRUE, penmat=NULL,
      wtvec=rep(1,n), dffactor=1,
      fdnames=list(NULL, dimnames(y)[[2]], NULL) ){
##
## 1.  fdPar
##
  fdP <- fdPar(fdobj, Lfdobj=Lfdobj, lambda=lambda,
               estimate=estimate, penmat=penmat)
##
## 2.  smooth.fd
##
  n <- length(argvals)
  w <- wtvec
#  smoothB <- smooth.basis(argvals, y, fdP, wtvec=w,
#        dffactor=dffactor, fdnames=fdnames)
  smooth.basis(argvals, y, fdP, wtvec=w,
        dffactor=dffactor, fdnames=fdnames)
}


  
