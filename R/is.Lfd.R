is.Lfd <- function(Lfd) {
#  check whether LFD is a linear differential operator
  if (inherits(Lfd, 'fd') || (is.numeric(Lfd) && Lfd >= 0))  {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
