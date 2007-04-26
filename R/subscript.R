"[.fd" <- function(fdobj, i=TRUE, j=TRUE, drop=FALSE) {
  #  select subsets of curves in a functional data object

  coef <- fdobj$coefs
  ndim <- length(dim(coef))

  if(ndim == 2) {
    coefselect <- coef[, i]
  } else {
    coefselect <- coef[, i, j]
  }
  fd(coefselect, fdobj$basis, names(fdobj$fdnames)) 
}
