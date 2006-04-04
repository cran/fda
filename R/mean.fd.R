mean.fd <- function(x, ...)
{
	coef      <- x$coefs
  	coefd     <- dim(coef)
  	ndim      <- length(coefd)
  	basisobj  <- x$basis
  	nbasis    <- basisobj$nbasis
 	if (ndim == 2) {
    	coefmean  <- matrix(apply(coef,1,mean),nbasis,1)
    	coefnames <- list(dimnames(coef)[[1]],"Mean")
  	} else {
    	nvar <- coefd[3]
    	coefmean  <- array(0,c(coefd[1],1,nvar))
    	for (j in 1:nvar) coefmean[,1,j] <- apply(coef[,,j],1,mean)
    	coefnames <- list(dimnames(coef)[[1]], "Mean", dimnames(coef)[[3]])
  	}
  	fdnames <- x$fdnames
  	fdnames[[2]] <- "mean"
  	fdnames[[3]] <- paste("mean",fdnames[[3]])
  	meanfd <- fd(coefmean, basisobj, fdnames)

  	meanfd
}
