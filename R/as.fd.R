as.fd <- function(x, ...) {
  UseMethod('as.fd')
}

as.fd.fdSmooth <- function(x, ...){
  x$fd
}

as.fd.function <- function(x, ...){
# Translate an object of class splinefun to class fd
##
## 1.  check class
##
  objName <- deparse(substitute(x))
  {
    if(length(objName)>1)
      objName <- character(0)
    else
      if(nchar(objName)>33)
        objName <- substring(objName, 1, 33)
  }
  if(!inherits(x, 'function')) 
    stop("'x' (", objName, ") is not of class function")
#
  xenv <- environment(x)
  xz <- get('z', xenv) 
  if(is.null(xz))
    stop("NULL environment of 'x' (", objName,
         ");  therefore, it can NOT have been created by 'splinefun.'")
#  
  if(is.null(xz$method))
    stop("'x' (", objName, ") has a NULL 'method', and therefore",
         " can NOT have been created by 'splinefun.'")
# z$method:  1=periodic, 2=natural, 3=fmm (std B-Splines, I believe)   
#  if(xz$method!=3){
  if(!(xz$method %in% 2:3)){
    msg <- paste("x (", objName, ") ", sep='')
    msg2 <- {
      if(xz$method=="1")
        paste(msg, " uses periodic B-splines, and as.fd ",
              "is programmed\n    to translate only B-splines ",
              "with coincident boundary knots.", sep='')
      else
        paste(msg, "does not use B-splines as required ",
                      "for function 'as.fd'.")
    }
    stop(msg2)
  }
##
## 2.  Create a basis 
##
  Knots <- xz$x
  y.x <- xz$y
  basisobj <- create.bspline.basis(range(Knots), breaks=Knots)
  print(basisobj)
  fdobj    <- fda::fd(matrix(0,basisobj$nbasis,1), basisobj)
  fdParobj <- fdPar(fdobj, lambda=0)
  nKn      <- length(Knots) 
  nobs     <- (2*nKn-1)
  x.       <- seq(Knots[1], Knots[nKn], length=nobs) 
  smooth.basis(x., x(x.), fdParobj)$fd
}

as.fd.smooth.spline <- function(x, ...){
# Translate an object of class smooth.spline to class fd
##
## 1.  check class
##
  objName <- deparse(substitute(x))
  {
    if(length(objName)>1)
      objName <- character(0)
    else
      if(nchar(objName)>33)
        objName <- substring(objName, 1, 33)
  }
  if(!inherits(x, 'smooth.spline')) 
    stop("'x' (", objName, ") is not of class smooth.spline")
##
## 2.  Create a basis 
##
  Kn0 <- x$fit$knot
  x0 <- min(x$x)
  x1 <- max(x$x) 
  Knots <- (x0+(x1-x0)*Kn0[4:(length(Kn0)-3)])
# Don't use 'unique' in case 'x' has coincident interior knots.
#  basis <- create.bspline.basis(breaks=Knots)
  basisobj <- create.bspline.basis(range(Knots), breaks=Knots)
#
  fd(x$fit$coef, basisobj)
}

