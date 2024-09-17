Data2fd <- function(argvals=NULL, y=NULL, basisobj=NULL, nderiv=NULL,
                    lambda=3e-8/diff(as.numeric(range(argvals))),
                    fdnames=NULL, covariates=NULL, method="chol") {
  
  #  Last modified 15 March 2024 by Jim
  
  #. Data2fd is a simplified version of other smoothing functions that allows
  #. neophytes or those in a hurry to smooth data successfully even though 
  #. one or more of its arguments is deficient in commonly occurring ways.
  
  #. First invoke function argvalSwap() to see if any arguments are 
  #. illegitimate, and repair them if possible.
  
  #. Five situations requiring modification or termination:
  #. 1. if argvals and y should be swapped
  #  2. if argvals is now NULL, build argvals and basisobj
  #  3. if the dimensions argvals and y as.array match, build basisobj
  #. 4. if length(dimy) < length(dima) swap argvals and y and their 
  #           dimensions dimy and dima
  #. 5. if basisobj has the wrong class, search other classes for the 
  #     basis object, and change.  Otherwise terminate.
  
  argChk <- argvalsySwap(argvals, y, basisobj)
  
  # check that argvals are numeric
  
  if(!is.numeric(AV <- argChk$argvals)){
    # terminal message
    if(is.null(AV)) stop('is.null(argChk$argvals); should be numeric')
    #. otherwise alert message
    cat('argChk$argvals is not numeric.\n')
    cat('class(argChk$argvals) = ', class(AV), '\n')
    print(AV)
  }  
  
  #. S3 object of class fdSmooth is returned by function smooth.basisPar()
  
  fdSmoothobj <- smooth.basisPar(argChk$argvals, argChk$y,
                                 fdobj=basisobj, Lfdobj=nderiv, lambda=lambda,
                                 fdnames=fdnames,
                                 covariates=covariates, method="chol")
  
  # return only the fd component
  
  return(fdSmoothobj$fd)
  
}

#  -------------------------------------------------------------------------

argvalsySwap = function(argvals=NULL, y=NULL, basisobj=NULL) {
  
#  Last modified 15 March 2024 by Jim
  
  if (inherits(basisobj, 'basisfd')) rangeval <- basisobj$rangeval
  
  ##. --------------------------------------------------------------------------
  ## 1.  If argument argvals is NULL, swap it with y
  ##. --------------------------------------------------------------------------
  
  if(is.null(y)){
    if(is.null(argvals)) stop("'y' is missing with no default")
    #   Store argvals as y and alert
    cat("'y' is missing, using 'argvals'\n") 
    y <- argvals
    argvals <- NULL 
  }
  
  ##. --------------------------------------------------------------------------
  ## 2.  If argument argvals is now NULL, then:
  ##        if basisobj is null, use default basis object
  ##        else if basisobj is numeric, default basis object
  ##                if   length of basis object is > 1 default basis object
  ##                else default basis object with order = numeric value
  ##             else if basisobj is "fd" basisobj = fd$basis
  ##                  else if basisobj is "fdPar" basisobj = fdPar$fd$basis
  ##                       else stop with message
  ##        set rangevalue to basisobj$rangeval but stop if is NULL
  ##        return with warning all three arguments
  ##. --------------------------------------------------------------------------

  dimy <- dim(as.array(y))
  if (is.null(argvals)) {
    {
      if(is.null(basisobj)){
        basisobj <- create.bspline.basis(basisobj)
      } else {
        if(is.numeric(basisobj)) {
          if(length(basisobj) > 1) {
            basisobj <- create.bspline.basis(basisobj)
          } else 
            basisobj <- create.bspline.basis(norder=basisobj)
        }
        else {
          if(inherits(basisobj, 'fd')) {
            basisobj <- basisobj$basis
          } else 
            if(inherits(basisobj, 'fdPar'))
              basisobj <- basisobj$fd$basis
        }
      }
    }
    rangeval <- basisobj$rangeval
    if(is.null(rangeval))
      stop('basisobj does not have a required rangeval component.')
    #    
    n <- dimy[1]
    #. alert message
    cat(paste("'argvals' is missing;  using seq(", rangeval[1],
              ", ", rangeval[2], ", length=", n, ")\n"))       
    argvals <- seq(rangeval[1], rangeval[2], length=n)
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  }
  
  ##. --------------------------------------------------------------------------
  ## 3.  dimy and dima are dimensions of argvals and y as array objects
  ##     If they match, proceed as in step 2 to construct basisobj
  ##     else stop with message
  ##. --------------------------------------------------------------------------

  dima <- dim(as.array(argvals))
  {
    if(length(dimy) == length(dima)){
      if(any(dimy != dima))
        #. terminal message
        stop("dimensions of 'argvals' and 'y' must be compatible;\n",
             "  dim(argvals) = ", paste(dima, collapse=' x '),
             ";  dim(y) = ", paste(dimy, collapse=' x ') )
      #     Check basisobj
      {
        if(inherits(basisobj, 'fd')) basisobj <- basisobj$basis
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
                if(is.null(basisobj)) {
                  basisobj <- create.bspline.basis(argvals)
                }
                else
                  if(!inherits(basisobj, 'basisfd'))
                    #. terminal message
                    stop("'basisobj' is NOT a functional basis",
                         " object (class 'basisfd');  class = ",
                         class(basisobj)[1])
              }
            }
          }
        }
      }
      arng     <- range(argvals)
      rangeval <- basisobj$rangeval
      if ((rangeval[1]<=arng[1]) && (arng[2]<=rangeval[2])) {
        return(list(argvals=argvals, y=y, basisobj=basisobj))
      }
      #
      yrng <- range(y)
      if((rangeval[1]<=yrng[1]) && (yrng[2]<=rangeval[2])) {
        #. alert message
        cat(paste("'argvals' is NOT contained in basisobj$rangeval",
                  ", but 'y' is;  swapping 'argvals' and 'y'.\n"))
        return(list(argvals=y, y=argvals, basisobj=basisobj)) 
      }
      #   Terminal message 
      stop("Neither 'argvals' nor 'y' are contained in ",
           "basisobj$rangeval")
    }
  }
  
  ##. --------------------------------------------------------------------------
  ## 4.  If length(dimy) < length(dima) swap argvals and y and their 
  ##           dimensions dimy and dima
  ##     Then stop if a value in dima is not in dimy 
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
  #  error message if argvals and y are inconsistent
  if(any(dima != dimy[1:length(dima)]))
    #  terminal message
    stop("A dimension of 'argvals' does not match 'y':\n",
         "  dim(argvals) = ", paste(dima, collapse=" x "),
         ";  dim(y) = ", paste(dimy, collapse=" x ") ) 
  
  ##. --------------------------------------------------------------------------
  ## 5.  check basisobj for having the wrong class, and is so 
  #      proceed as above to find an object with the right class
  ##. --------------------------------------------------------------------------
  
  {
    if(inherits(basisobj, 'fd')) basisobj <- basisobj$basis
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
                #. error message if basisobj incorrect class
                stop("'basisobj' is NOT a functional basis",
                     " object (class 'basisfd');  class = ",
                     class(basisobj)[1])
          }
        }
      }
    }
    rangeval <- basisobj$rangeval
  }

  ##. --------------------------------------------------------------------------
  ## 6.  Check compatibility of argvals with basisobj$rangeval
  ##. --------------------------------------------------------------------------
  
  a01  <- basisobj$rangeval
  arng <- range(argvals)
  if ((a01[1] <= arng[1]) && (arng[1] <= a01[2])) {
    return(list(argvals=argvals, y=y, basisobj=basisobj))
  }
  #  error message 
  stop("There are argvals not contained within basisobj$rangeval")
  
}

