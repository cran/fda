center.fd <- function(fd)
{
  #  remove mean function for functional observations

  #  Last modified 6 Feb 2001

  if (!(inherits(fd, "fd"))) stop('Argument  FD not a functional data object.')

  coef   <- as.array(fd[[1]])
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  basis  <- fd[[2]]
  nbasis <- basis$nbasis
  if (ndim == 2) {
    coefmean <- apply(coef,1,mean)
    coef     <- sweep(coef,1,coefmean)
  } else {
    nvar <- coefd[3]
    for (j in 1:nvar) {
      coefmean <- apply(coef[,,j],1,mean)
      coef[,,j] <- sweep(coef[,,j],1,coefmean)
    }
  }
  fdnames <- fd$fdnames
  names(fdnames)[3] <- paste('Centered',names(fdnames)[3])
  centerfd <- create.fd(coef, basis, fdnames)
  return(centerfd)
}
