is.fd <- function(fd) {
#  check whether FD is a functional data object

#  Last modified 20 Feb 2003

  if (inherits(fd, "fd")) return(TRUE) else return(FALSE)
}
