rangechk <- function(rangeval) {
#  check a range vector argument

#  last modified 16 May 1999

  if (!is.vector(rangeval))       return(FALSE)
  if (length(rangeval) != 2)      return(FALSE)
  if (rangeval[1] >= rangeval[2]) return(FALSE)
  return(TRUE)
}
