getnames <- function(fd) {
  #  Extracts the fdnames from a functional data object FD

  #  Last modified 6 Feb 2001

  if (inherits(fd, "fd")) {
    fdnames <- fd[[3]]
  } else {
    stop("An object of class fd or basis.fd expected, but not found.")
  }
  return(fdnames)
}
