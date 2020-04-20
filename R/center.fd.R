center.fd <- function(fdobj)
{
#  remove mean function for functional observations

#  Last modified 16 January 2020

if (!(is.fd(fdobj) || is.fdPar(fdobj))) 
  stop("First argument is neither an fd or an fdPar object.")
if (is.fdPar(fdobj)) fdobvj = fdobj$fd

coef     <- as.array(fdobj$coefs)
coefd    <- dim(coef)
ndim     <- length(coefd)
basisobj <- fdobj$basis
nbasis   <- basisobj$nbasis
if (ndim == 2) {
   coefmean <- apply(coef,1,mean)
   coef     <- sweep(coef,1,coefmean)
} else {
   nvar <- coefd[3]
   for (j in 1:nvar) {
       coefmean  <- apply(coef[,,j],1,mean)
       coef[,,j] <- sweep(coef[,,j],1,coefmean)
   }
}
fdnames      <- fdobj$fdnames
fdnames[[3]] <- paste("Centered",fdnames[[3]])
centerfdobj  <- fd(coef, basisobj, fdnames)
return(centerfdobj)	
}
