"[.fd" <- function(fd, ..., drop=FALSE)
{
  #  select subsets of curves in a functional data object
  
  #  last modified 7 February 2001
  
  my.call  <- match.call()
  nargvals <- nargs() - !missing(drop)
  nind <- nargvals - 1
  if (nind > 2) stop("No more than two subscripts allowed.")
  coef <- getcoef(fd)
  ndim <- length(dim(coef))
  if (ndim == 1) stop("Subscripting not allowed for a single function.")
  fdnames <- fd$fdnames
  ind <- ..1
  if (nind == 1) {
    if (ndim == 2) {
      coefselect <- as.matrix(coef[,ind])
      if (is.null(fdnames[[2]]) != TRUE) fdnames[[2]] <- fdnames[[2]][ind]
    }
    if (ndim == 3) {
      if (nargvals == 3) {
        coefselect <- coef[,ind,,drop=FALSE]
        if (is.null(fdnames[[2]]) != TRUE) fdnames[[2]] <- fdnames[[2]][ind]
        fdnames[[2]] <- fdnames[[2]][ind]
      } else {
        coefselect <- coef[,,ind]
        fdnames[[3]] <- fdnames[[3]][ind]
      }
      if (dim(coefselect)[3] == 1) coefselect <- coefselect[,,1]
    }
  }
  if (nind == 2) {
    if (ndim <  3) stop("Two subscripts only allowed for multiple variables.")
    coefselect <- coef[,..1 ,..2,drop=FALSE]
    if (dim(coefselect)[3] == 1) coefselect <- coefselect[,,1]
    if (is.null(fdnames[[2]]) != TRUE) fdnames[[2]] <- fdnames[[2]][ind]
    fdnames[[3]] <- fdnames[[3]][..2]
  }
  basisfd  <- getbasis(fd)
  fdselect <- create.fd(coefselect, basisfd, fdnames=fdnames)
  return(fdselect)
}
