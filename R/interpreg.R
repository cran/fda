interpreg <- function(x, y, xeval, periodic)
{
  #  This is a wrapper function for the Fortran function interpolate.

  n      <- length(x)
  neval  <- length(xeval)
  ier    <- 0
  yeval  <- rep(0,neval)
  eeval  <- yeval
  ier    <- 0
  if (any(is.na(x))) {
    warning ("X has missing values in INTERPOLATE")
    ier <- 1
  }
  if (any(is.na(y))) {
    warning ("Y has missing values in INTERPOLATE")
    ier <- 1
  }
  if (min(x[2:n]-x[1:(n-1)]) <= 0) {
    warning ("X is not strictly increasing in INTERPOLATE")
    ier <- 1
  }
  if (min(xeval) < x[1] | max(xeval) > x[n]) {
    warning ("XEVAL is out of range in INTERPOLATE")
    IER <- 1
  }
  if (ier != 0) {
    warning("Further processing aborted in INTERPOLATE")
    return(rep(0,neval))
  }

  if (periodic) iper <- 1 else iper <- 0

  result <- .Fortran("interpreg",
         as.integer(n),
         as.double(x),      as.double(y),
         as.integer(neval), as.double(xeval), as.double(yeval),
         as.integer(iper),  as.integer(ier) )

  yeval <- result[[6]]
  ier   <- result[[8]]
  if (ier != 0) warning(paste("IER =",ier," in INTERPOLATE"))

  return (yeval)
}
