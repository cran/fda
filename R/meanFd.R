meanFd <- function(fd)
{
  #  Compute mean functional data object for functional observations
  #    in argument FD

  #  Last modified 6 Feb 2001

  if (!(inherits(fd, "fd"))) stop("Argument  FD not a functional data object.")

  coef   <- getcoef(fd)
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  basis  <- getbasis(fd)
  nbasis <- basis$nbasis
  if (ndim == 2) {
    coefmean  <- matrix(apply(coef,1,mean),nbasis,1)
    coefnames <- list(dimnames(coef)[[1]],"Mean")
  } else {
    nvar <- coefd[3]
    coefmean  <- array(0,c(coefd[1],1,nvar))
    for (j in 1:nvar) coefmean[,1,j] <- apply(coef[,,j],1,mean)
    coefnames <- list(dimnames(coef)[[1]], "Mean", dimnames(coef)[[3]])
  }
  fdnames <- getnames(fd)
  fdnames[[2]] <- "1"
  names(fdnames)[2] <- "Mean"
  names(fdnames)[3] <- paste("Mean",names(fdnames)[3])
  meanfd <- create.fd(coefmean, basis, fdnames)

  return(meanfd)
}
