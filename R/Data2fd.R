Data2fd <- function(argvals=NULL, y=NULL, basisobj=NULL, nderiv=NULL,
                    lambda=3e-8/diff(as.numeric(range(argvals))),
                    fdnames=NULL, covariates=NULL, method="chol") {
  
  #  Last modified 1 March 2024
  
  #. Tests five situations requiring modification or termination.
  
  argChk <- argvalsySwap(argvals, y, basisobj)
  
  if(!is.numeric(AV <- argChk$argvals)){
    if(is.null(AV))
      stop('is.null(argChk$argvals); should be numeric')
    cat('argChk$argvals is not numeric.\n')
    cat('class(argChk$argvals) = ', class(AV), '\n')
  }
  
  #. Success and smoothing ... S3 object of class fdSmooth is returned
  
  smBasis <- smooth.basisPar(argChk$argvals, argChk$y,
                                 fdobj=basisobj, Lfdobj=nderiv, lambda=lambda,
                                 fdnames=fdnames,
                                 covariates=covariates, method="chol")
  return(smBasis$fd)
  
}

#  -------------------------------------------------------------------------

## 2020-01-16:  Spencer Graves makes argvalsySwap 
## An internal function that tests for 6 situations that require modification
## with a warning, or a terminal error message

argvalsySwap = function(argvals=NULL, y=NULL, basisobj=NULL) {
  
  if (inherits(basisobj, 'basisfd')) rangeval <- basisobj$rangeval
  
  ##. --------------------------------------------------------------------------
  ## Section 1.  if(is.null(y)) use argvals for y
  ##. --------------------------------------------------------------------------
  
  if(is.null(y)){
    if(is.null(argvals)) stop("'y' is missing with no default")
    #   Store argvals as y and alert
    cat("'y' is missing, using 'argvals'\n") 
    y       <- argvals
    argvals <- NULL 
  }
  
  ##. --------------------------------------------------------------------------
  ## Section 2.  test for missing argvals, if so construct a sequence
  ##. --------------------------------------------------------------------------

  dimy <- dim(as.array(y))
  if(is.null(argvals)) {
    # the following code block is run if TRUE
    {  # beginning of code block
      if(is.null(basisobj)){
        basisobj <- create.bspline.basis(basisobj)
      } else {
        if(is.numeric(basisobj)) {
          if(length(basisobj)>1){
            basisobj <- create.bspline.basis(basisobj)
          } else
            basisobj <- create.bspline.basis(norder=basisobj)
        }
        else {
          if(inherits(basisobj, 'fd')){
            basisobj <- basisobj$basis
          } else
            if(inherits(basisobj, 'fdPar'))
              basisobj <- basisobj$fd$basis
        }
      }
    }  #. end of code block
    # This is executed whether or not the previous was
    #. locate the range from basisobj
    a01 <- basisobj$rangeval
    #  if range is null, error message and stop
    if(is.null(a01))
      stop('basisobj does not have a required rangeval component.')
    n <- dimy[1]
    cat(paste("'argvals' is missing;  using seq(", a01[1],
              ", ", a01[2], ", length=", n, ")\n"))
    #  construct the argval sequence
    argvals <- seq(a01[1], a01[2], length=n)
    
    return(list(argvals=argvals, y=y, basisobj=basisobj))
    
  }
  
  ##. --------------------------------------------------------------------------
  ## 3.  swapping y and argvals 
  ##. --------------------------------------------------------------------------

  dima <- dim(as.array(argvals)) 
  { # First line in code block
    if(length(dimy) == length(dima)) {
      if(any(dimy != dima))
        stop("dimensions of 'argvals' and 'y' must be compatible;\n",
             "  dim(argvals) = ", paste(dima, collapse=' x '),
             ";  dim(y) = ", paste(dimy, collapse=' x ') )
      #     Check basisobj
      { # First line in code block
        if(inherits(basisobj, 'fd')) basisobj <- basisobj$basis
        else {
          if(inherits(basisobj, 'fdPar'))
            basisobj <- basisobj$fd$basis
          else {
            if(inherits(basisobj, 'array')) {
              fd.      <- fd(basisobj)
              basisobj <- fd.$basis
            } else { 
              if(inherits(basisobj, 'integer'))
                basisobj <- create.bspline.basis(argvals, norder=basisobj)
              else {
                if(is.null(basisobj))
                  basisobj <- create.bspline.basis(argvals)
                else
                  if(!inherits(basisobj, 'basisfd'))
                    stop("'basisobj' is NOT a functional basis",
                         " object (class 'basisfd');  class = ",
                         class(basisobj)[1])
              }
            }
          }
        }
      }  # Last line in code block
      a01  <- basisobj$rangeval
      arng <- range(argvals)
      if ((rangeval[1]<=arng[1]) && (arng[2]<=rangeval[2])) {
        return(list(argvals=argvals, y=y, basisobj=basisobj))
      }
      #
      yrng <- range(y)
      if ((a01[1]<=yrng[1]) && (yrng[2]<=a01[2])) {
        cat(paste("'argvals' is NOT contained in basisobj$rangeval",
                  ", but 'y' is;  swapping 'argvals' and 'y'.\n"))
        return(list(argvals=y, y=argvals, basisobj=basisobj)) 
      }
      #   Terminal message 
      stop("Neither 'argvals' nor 'y' are contained in ",
           "basisobj$rangeval")
    }
  } # Last line in code block
  
  ##. --------------------------------------------------------------------------
  ## 4.  If(length(dimy) < length(dima)) swap argvals and y
  ##. --------------------------------------------------------------------------
  
  if(length(dimy)<length(dima)) {
    cat(paste("Swapping 'y' and 'argvals', because 'y' is ",
              "simpler,\n  and 'argvals' should be;  now ",
              "dim(argvals) = ", paste(dimy, collapse=" x "),
              ";  dim(y) = ", paste(dima, collapse=" x "),"\n" )) 
    y. <- argvals
    argvals <- y
    y <- y.
    #
    d. <- dima
    dima <- dimy
    dimy <- d.
  }   
  
  if(any(dima != dimy[1:length(dima)]))
    stop("A dimension of 'argvals' does not match 'y':\n",
         "  dim(argvals) = ", paste(dima, collapse=" x "),
         ";  dim(y) = ", paste(dimy, collapse=" x ") ) 
  
  ##. --------------------------------------------------------------------------
  ## 5.  Check compatibility of argvals with basisobj
  ##. --------------------------------------------------------------------------
  
  { # First line in code block
    if(inherits(basisobj, 'fd'))basisobj <- basisobj$basis
    else {
      if(inherits(basisobj, 'fdPar'))
        basisobj <- basisobj$fd$basis
      else {
        if(inherits(basisobj, 'array')){
          fd. <- fd(basisobj)
          basisobj <- fd.$basis
        }
        else { 
          if(inherits(basisobj, 'integer'))
            basisobj <- create.bspline.basis(argvals, norder=basisobj)
          else {
            if(is.null(basisobj))
              basisobj <- create.bspline.basis(argvals)
            else
              if(!inherits(basisobj, 'basisfd'))
                stop("'basisobj' is NOT a functional basis",
                     " object (class 'basisfd');  class = ",
                     class(basisobj)[1])
          }
        }
      }
    }
  } # Last line in code block
  a01 <- basisobj$rangeval
  arng <- range(argvals)
  if((a01[1]<=arng[1]) && (arng[2]<=a01[2])) {
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  }
  #
  stop("'argvals' are not contained in basisobj$rangeval")
}

