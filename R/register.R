register <- function(x, y0, y, Dy, D2y,
                     nbreak    = 6,
                     wt        = c(0.5,rep(1, n-2),0.5),
                     breakvals = seq(x[1],x[n],length=nbreak),
                     cvec0     = rep(0,nbreak),
                     lambda    = 0,    iterlim = 20, dbglev = 1,
                     conv      = 1e-3, periodic = FALSE, critno = 2)
{

  #  This is a wrapper function for the Fortran function register.f.

  #  Arguments are:

  #  X         ...  array of N argument values
  #  Y0        ...  N by NVAR array of target function values to fit
  #  Y         ...  N by NVAR array of function values
  #                   for function to be registered
  #  DY        ...  N by NVAR array of first  derivative values
  #  D2Y       ...  N by NVAR array of second derivative values
  #  NBREAK    ...  number of break values defining the monotone warping fn.
  #                 If there is a conflict with the length of BREAKVALS,
  #                 the length of BREAKVALS will be used.
  #  WT        ...  array of length N of weights to be applied to function
  #                 values
  #  BREAKVALS ...  array of break values containing all values of X
  #  CVEC0     ...  NBREAK - 1 starting values for coefficients
  #  LAMBDA    ...  penalty parameter
  #  ITERLIM   ...  number of iteration limit
  #  DBGLEV    ...  level of output of computation history
  #  CONV      ...  convergence criterion
  #  PERIODIC  ...  if T treat the curves as periodic
  #  CRITNO    ...  if 1 least squares, if 2 log eigenvalue ratio

  #  A list is returned with the following elements:

  #  cvecstore  ... the parameters determine the transformations
  #  deltastore ... the estimated shift value for each curve is PERIODIC is T
  #  yreg       ... an array of the same size as Y of registered curve values
  #  hfun       ... warping function values, with domain and range the
  #                 same as the domain of Y

  #  check argument dimensions

  n       <- length(x)

  if (length(dim( y0)) > 2) stop(" Y0 has more than two dimensions")
  if (length(dim(  y)) > 2) stop("  Y has more than two dimensions")
  if (length(dim( Dy)) > 2) stop(" DY has more than two dimensions")
  if (length(dim(D2y)) > 2) stop("D2Y has more than two dimensions")

  y0  <- as.matrix(y0)
  y   <- as.matrix(y)
  Dy  <- as.matrix(Dy)
  D2y <- as.matrix(D2y)

  nvar <- ncol(y0)

  if (ncol(  y) != nvar) stop("  Y must have same no. cols. as Y0")
  if (ncol( Dy) != nvar) stop(" DY must have same no. cols. as Y0")
  if (ncol(D2y) != nvar) stop("D2Y must have same no. cols. as Y0")

  if (nrow(  y) != n) stop("  Y must have same no. row as Y0")
  if (nrow( Dy) != n) stop(" DY must have same no. row as Y0")
  if (nrow(D2y) != n) stop("D2Y must have same no. row as Y0")

  nbreaktemp <- length(breakvals)
  if (nbreak != nbreaktemp) nbreak <- nbreaktemp
  if (length(cvec0) != nbreak) {
       stop("BREAKVALS and CVEC0 lengths are inconsistent") }

  #  put X and BREAKVALS into the unit interval

  span      <- breakvals[nbreak] - breakvals[1]
  breaknorm <- (breakvals - breakvals[1])/span
  xnorm     <- (x         - breakvals[1])/span
  Dynorm    <- Dy*span
  D2ynorm   <- D2y*span*span

  #  check the derivatives to see if they are consistent

  for (ivar in 1:nvar) {
    ratio1 <- derivchk(xnorm,      y[,ivar],  Dynorm[,ivar])
    ratio2 <- derivchk(xnorm, Dynorm[,ivar], D2ynorm[,ivar])
#   print(c(ratio1,ratio2))
    if (ratio1 >= .3 | ratio2 >= .5) warning(
         "Derivatives appear to be inconsistent")
  }

  #  set up arrays for storing transformation values and coefficient values

  ncvec   <- nbreak-1
  ioutlev <- 0
  ier     <- 0
  iter    <- 0
  iconv   <- 0
  history <- matrix(0,2,iterlim+1)
  hfun    <- rep(0,n)
  yreg    <- y

  result <- .Fortran("register",
             as.integer(n),        as.integer(nvar),
             as.double(xnorm),     as.double(y0),
             as.double(y),         as.double(Dynorm),    as.double(D2ynorm),
             as.double(wt),
             as.integer(nbreak),   as.double(breaknorm), as.double(cvec0),
             as.integer(iterlim),  as.integer(iter),
             as.double(hfun),      as.double(yreg),
             as.double(lambda),    as.integer(periodic), as.integer(critno),
             as.double(conv),      as.integer(ioutlev),  as.integer(iconv),
             as.double(history),   as.integer(ier) )

  ier   <- result[[23]]
  if (ier != 0) stop (paste("IER =",ier))
  iconv <- result[[21]]
  if (iconv != 0) warning ("No convergence")
  if (dbglev != 0) {
    iter    <- result[[13]]
    history <- matrix(result[[22]],2,iterlim+1)[,1:(iter+1)]
    cat("Iter.  Crit.     Grad Len.\n")
    for (i in 0:iter)
     cat(format(i),"      ",format(round(history[,i+1],7)),"\n")
  }
  cvec       <- result[[11]]
  deltastore <- cvec[1]
  cvecstore  <- cvec[2:nbreak]
  hfun       <- result[[14]]
  yreg       <- matrix(result[[15]],n,nvar)

  #  compute the registered curve

  return ( list( cvecstore, deltastore, yreg, hfun ) )

}
