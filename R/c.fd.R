'c.fd'<- function(...)
{
#
#   concatenates a number of .fd objects.  It is assumed that all the
#   objects have the same basisfd objects, and that all the coef arrays
#   have the same number of dimensions
#

#  Last modified 6 Feb 2001

  fdlist <- list(...)
  n      <- length(fdlist)
  fd1    <- fdlist[[1]]
  if(n == 1) return(fd1)
  coef    <- getcoef(fd1)
  coefd   <- dim(coefs)
  ndim    <- length(dimcoef)
  basisfd <- getbasis(fd1)
  fdnames <- getnames(fd1)
  #  check that the fd objects are consistent with each other
  if(inherits(fd1, "fd")) stop('Objects must be of class fd')
  for(j in (2:n)) {
    fdj <- fdlist[[j]]
    if(inherits(fdj, "fd")) stop('Objects must be of class fd')
    if(any(unlist(fdj[[2]]) != unlist(basisfd)))
      stop('Objects must all have the same basis')
    if(length(dim(fdj[[1]])) != ndim)
      stop('Objects must all have the same number of multiple functions')
  }
  #  concatenate by concatenate coefficient matrices
  if(ndim == 2) {
    for (j in 2:n) {
      nameslist <- dimnames(coef)
      fdj       <- fdlist[[j]]
      coefj     <- getcoef(fdj)
      coef      <- cbind(coef, coefj)
      nameslist[[2]] <- c(nameslist[[2]], dimnames(coefj)[[2]])
    }
  } else {
    for(j in (2:n)) {
      nameslist <- dimnames(coef)
      fdj       <- fdlist[[j]]
      coefj     <- getcoef(fdj)
      coef      <- c(coef, aperm(coefj, c(1, 3, 2)))
      nameslist[[2]] <- c(nameslist[[2]], dimnames(coefj)[[2]])
    }
    dim(coef) <- c(coefd[1], coefd[3], length(coef)
                  /(coefd[1] * coefd[3]))
    coef <- aperm(coef, c(1, 3, 2))
  }
  dimnames(coef) <- nameslist
  concatfd <- create.fd(coefs, basisfd, fdnames)
  return(concatfd)
}
