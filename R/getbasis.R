getbasis <- function(fd) {
  #  Extracts the basis.fd object from a functional data object FD, or,
  #    if FD is already a basis object, just returns it.

  #  Last modified 4 July 2001
  
  if (inherits(fd, "fd")) {
    basisfd <- fd[[2]]
  }  else {
    if (inherits(fd, "basis.fd")) {
      basisfd <- fd
    } else {
      stop("An object of class fd or basis.fd expected, but not found.")
    }
  }
  return(basisfd)
}
