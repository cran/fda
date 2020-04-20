symsolve <- function(Asym, Bmat)
{
  #  solves the system ASYM X = BMAT for X where ASYM is symmetric
  #  returns X
  n <- ncol(Asym)
  if (any(is.na(Asym))) stop("Asym has NA values.")
  if (max(abs(Asym-t(Asym)))/max(abs(Asym)) > 1e-10) {
    stop('Argument not symmetric.')
  } else {
    Asym <- (Asym + t(Asym))/2
  }
  Lmat <- chol(Asym)
  temp <- solve(t(Lmat),Bmat)
  Xmat <- backsolve(Lmat,temp)
  return(Xmat)
}
