sumFd <- function(fd)
{
  #  Compute sum function for functional observations

  #  Last modified 6 Feb 2001

  if (!(inherits(fd, "fd"))) stop("Argument FD not a functional data object.")

  coef   <- getcoef(fd)
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  basis  <- getbasis(fd)
  nbasis <- basis$nbasis
  if (ndim == 2) {
    coefsum   <- matrix(apply(coef,1,sum),nbasis,1)
    coefnames <- list(dimnames(coef)[[1]],"Sum")
  } else {
    nvar <- coefd[3]
    coefsum  <- array(0,c(coefd[1],1,nvar))
    for (j in 1:nvar) coefsum[,1,j] <- apply(coef[,,j],1,sum)
    coefnames <- list(dimnames(coef)[[1]], "Sum", dimnames(coef)[[3]])
  }
  fdnames <- getnames(fd)
  fdnames[[2]] <- "1"
  names(fdnames)[2] <- "Sum"
  names(fdnames)[3] <- paste("Sum",names(fdnames)[3])
  sumfd <- create.fd(coefsum, basis, fdnames)

  return(sumfd)
}
