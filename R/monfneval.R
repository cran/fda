monfneval <- function(xeval, breakvals, cvec)
{
  #  XEVAL     ... strictly increasing sequence of argument values at
  #                which monotone fn. is to be evaluated.
  #                range of XEVAL must be within range of BREAKVALS
  #  BREAKVALS ... strictly increasing sequence of break values,
  #                first is lower boundary of interval, last is upper boundary
  #  CVEC      ... vector of coefficients of hat functions

  neval <- length(xeval)

  if (neval > 1 & min(xeval[2:neval]-xeval[1:(neval-1)]) <= 0)
      stop("Arguments must be strictly increasing")

  nbreak <- length(breakvals)
  if (nbreak < 2) stop("At least two break values required.")

  if (xeval[1] < breakvals[1])
       stop("Smallest argument out of range.")
  if (xeval[neval] > breakvals[nbreak])
       stop("Largest argument out of range.")

  #  put XEVAL and BREAKVALS into the unit interval

  span      <- breakvals[nbreak] - breakvals[1]
  breaknorm <- (breakvals - breakvals[1])/span
  xnorm     <- (xeval     - breakvals[1])/span

  feval <- rep(0, neval)
  ier   <- 0

  result <- .Fortran("monfneval", as.integer(neval),  as.double(xnorm),
                              as.integer(nbreak), as.double(breaknorm),
                              as.double(cvec),    as.double(feval),
                              as.integer(ier) )
  feval <- result[[6]]
  ier   <- result[[7]]

  if (ier != 0)
    warning(c("Nonzero value of IER returned from Fortran subroutine", ier))

  return (feval)
}
