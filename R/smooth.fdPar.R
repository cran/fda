smooth.fdPar <- function(fdobj, Lfdobj=NULL,
         lambda=1e-4, estimate=TRUE, penmat=NULL){
##
## 1.  fdPar
##
  fdP <- fdPar(fdobj, Lfdobj=Lfdobj, lambda=lambda,
               estimate=estimate, penmat=penmat)
##
## 2.  smooth.fd
##
  return(smooth.fd(fdobj, fdP))
}
