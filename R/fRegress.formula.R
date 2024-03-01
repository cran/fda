fRegress.character <- function(y, data=NULL, betalist=NULL,
                               wt=NULL, y2cMap=NULL, SigmaE=NULL,
                               method='fRegress',
                               sep='.', ...) {
  fRegress.formula(y=y, data=data, betalist=betalist,
                   wt=wt, y2cMap=y2cMap, SigmaE=SigmaE,
                   method=method, sep=sep, ...)
}

fRegress.formula <- function(y, data=NULL, betalist=NULL,
                             wt=NULL, y2cMap=NULL, SigmaE=NULL,
                             method='fRegress', sep='.', ...) {
  
  #  FREGRESS.FORMULA  Fits a scalar dependent variable using the concurrent
  #                    functional regression model using inner products
  #                    of functional covariates and functional regression
  #                    functions.
  #
  #  Arguments:
  #  Y        ... A numeric vector or a functional data object that is the 
  #               dependent variable..
  #  DATA     ... a list object of length p with each list
  #               containing an object for an independent variable.
  #               the object may be:
  #                   a functional data object or
  #                   a vector
  #               if XFDLIST is a functional data object or a vector,
  #               it is converted to a list of length 1.
  #  BETALIST ... a list object of length p with each list
  #               containing a functional parameter object for
  #               the corresponding regression function.  If any of
  #               these objects is a functional data object, it is
  #               converted to the default functional parameter object.
  #               if BETALIST is a functional parameter object
  #               it is converted to a list of length 1.
  #  WT       ... a vector of nonnegative weights for observations
  #  Y2CMAP   ... the matrix mapping from the vector of observed values
  #               to the coefficients for the dependent variable.
  #               This is output by function SMOOTH_BASIS.  If this is
  #               supplied, confidence limits are computed, otherwise not.
  #  SIGMAE   ... Estimate of the covariances among the residuals.  This
  #               can only be estimated after a preliminary analysis
  #. METHOD   ... Type of object returned:  
  #                 if "fRegress" results of analysis are returned
  #                 if "model" an argument list fRegressList is returned
  
  #  Last modified 29 January 2024 by Jim Ramsay
  
  ##
  ## 1.  get y = left hand side of the formula
  ##
  
  Formula <- y                    #  character vector
  yName   <- Formula[[2]]         #  name of dependent variable object
  yNm     <- as.character(yName)  #  character strong for name
  #  check name of dependent variable
  if(!inherits(yName, 'name') || (length(yNm) > 1))
    stop('The left hand side of formula must be a simple object; ',
         ' instead, LHS = ', as.character(Formula)[2],
         ', which has class = ', class(yName))
  
  ##
  ## 2.  check the data argument
  ##
  
  dataNames <- names(data)
  #  extract dependent variable 
  yvec <- {
    if (yNm %in% dataNames) data[[yNm]] else get(yNm)  
  }
  #  get the range of the  dependent variable  and 
  #  obtain the dimensions of the coefficient matrix if yvec is of class 'fd'
  #  obtain the dimensions of the yvec if yvec is of class 'numeric'
  trng <- NULL
  {
    if(inherits(yvec, 'fd')){
      ydim <- dim(yvec$coefs)
      if(is.null(ydim) || (length(ydim)<2)) {
        yvec$coefs <- as.matrix(yvec$coefs)
        ydim   <- dim(yvec$coefs)
      }
      ny   <- ydim[2]
      trng <- yvec$basis$rangeval
    } else{
      if(inherits(yvec, 'numeric')){
        ydim <- dim(yvec)
        if(is.null(ydim))
          ny <- length(yvec)
        else
          ny <- ydim[1]
      }
      else
        stop('The left hand side of formula must have class ',
             'numeric or fd;  instead is ', class(yvec))
    }
  }
  
  ##
  ## 3.  check the formula for excessive complexity
  ##
  
  allVars  <- all.vars(Formula)
  xNms     <- allVars[allVars != yNm]
  Terms    <- terms(Formula)
  termLbls <- attr(Terms, 'term.labels')
  oops     <- length(which(!(termLbls %in% xNms))) > 0
  if(oops) stop('formula contains terms that fRegress can not handle; ',
                ' the first one = ', termLbls[[1]])
  #
  k1 <- length(allVars)
  type <- rep(NA,k1)
  names(type) <- allVars
  nbasis      <- type
  if(inherits(yvec, 'fd')){
    type[1] <- yvec$basis$type
    nb      <- yvec$basis$nbasis
    if(!is.null(nb)) nbasis[1] <- nb
  }
  
  ##
  ## 4.  Inventory the right hand side
  ##
  
  k0              <- length(xNms)
  xfdlist0        <- vector('list', k0)
  names(xfdlist0) <- xNms
  xNames          <- xfdlist0
  nVars           <- rep(NA, k0)
  names(nVars)    <- xNms
  oops <- FALSE
  for(i in 1:k0) {
    xNm <- xNms[i]
    xi <- {
      if(xNm %in% dataNames) data[[xNm]] else get(xNm)
    }
    {
      if(class(xi) %in% c('fd', 'fdPar')) {
        xj <- {
          if(inherits(xi, 'fd')) xi
          else xi$fd
        }
        xrng <- xj$basis$rangeval
        {
          if(is.null(trng))
            trng <- xrng
          else
            if(any(xrng != trng)){
              oops <- TRUE
              cat('incompatible rangeval found in ', xNm,
                  '$rangeval = ', paste(xrng, collapse=', '),
                  ' != previous = ', paste(trng, collapse=', '),
                  sep='')
            }
        }
        xdim <- dim(xj$coefs)
        {
          if(is.null(xdim) || (length(xdim)<2)){
            xj$coefs <- as.matrix(xj$coefs)
            xdim <- dim(xj$coefs)
            nxi <- xdim[2]
            nVars[i] <- 1
            xNames[[i]] <- xNm
          }
          else {
            if(length(xdim)<3){
              nxi <- xdim[2]
              nVars[i] <- 1
              xNames[[i]] <- xNm
            }
            else {
              nxi <- xdim[2]
              if(length(xdim)<4){
                nVars[i] <- xdim[3]
                xNmsi <- dimnames(xj$coefs)[[3]]
                {
                  if(is.null(xNmsi))
                    xNames[[i]] <- paste(xNm, 1:xdim[3], sep=sep)
                  else
                    xNames[[i]] <- paste(xNm, xNmsi, sep=sep)
                }
              }
              else {
                oops <- TRUE
                cat(xNm, 'has too many levels:  dim(x$coefs) =',
                    paste(xdim, collapse=', '))
              }
            }
          }
        }
        type[i+1] <- xj$basis$type
        nb <- xj$basis$nbasis
        if(!is.null(nb))nbasis[i+1] <- nb
        xfdlist0[[i]] <- xi
      }
      else {
        if(is.numeric(xi)){
          xdim <- dim(xi)
          {
            if(is.null(xdim) || (length(xdim)<2)){
              nxi <- length(xi)
              nVars[i] <- 1
              xNames[[i]] <- xNm
            }
            else {
              nxi <- xdim[1]
              {
                if(length(xdim)<3){
                  nVars[i] <- xdim[2]
                  xNmsi <- dimnames(xi)[[2]]
                  {
                    if(is.null(xNmsi))
                      xNames[[i]] <- paste(xNm, 1:xdim[2], sep=sep)
                    else
                      xNames[[i]] <- paste(xNm, xNmsi, sep=sep)
                  }
                }
                else{
                  oops <- TRUE
                  cat(xNm, 'has too many levels:  dim(x) =',
                      paste(xdim, collapse=', '))
                }
              }
            }
          }
          xfdlist0[[i]] <- xi
        }
        else {
          if(inherits(xi, 'character'))
            xi <- factor(xi)
          {
            if(inherits(xi, 'factor')) {
              f.i <- formula(paste('~', xNm))
              Xi.df <- data.frame(xi)
              names(Xi.df) <- xNm
              Xi <- (model.matrix(f.i, Xi.df)[, -1, drop=FALSE])
              nxi <- dim(Xi)[1]
              xiNms <- dimnames(Xi)[[2]]
              nVars[i] <- length(xiNms)
              xNmLen <- nchar(xNm)
              xiLvls <- substring(xiNms, xNmLen+1)
              xNames[[i]] <- paste(xNm, xiLvls, sep=sep)
              xfdlist0[[i]] <- Xi
            }
            else{
              oops <- TRUE
              cat('ERROR:  variable', xNm, 'must be of class',
                  'fd, numeric, character or factor;  is', class(xi))
              nxi <- length(xi)
            }
            }
        }
      }
    }
    if(nxi != ny){
      cat('ERROR:  variable', xNm, 'has only',
          nxi, 'observations !=', ny,
          '= the number of observations of yvec.')
      oops <- TRUE
    }
  }
  if (oops) stop('illegal variable on the right hand side.')
  # If no functions found:
  if(is.null(trng)){
    warning("No functions found;  setting rangeval to 0:1")
    trng <- 0:1
  }
  
  ##
  ## 5.  Create xfdlist
  ##
  
  xL.L0   <- rep(1:k0, nVars)
  xNames. <- c('const', unlist(xNames))
  k <- 1+sum(nVars)
  xfdlist <- vector('list', k)
  names(xfdlist) <- xNames.
  #  create constfd for the intercept
  #  xfdlist[[1]] <- create.constant.basis(trng)
  xfdlist[[1]] <- rep(1, ny)
  i1 <- 1
  for(ix in 1:k0) {
    i0  <- i1+1
    xNm <- xNms[ix]
    xi  <- xfdlist0[[ix]]
    {
      if(inherits(xi, 'fd')) {
        if(nVars[ix] < 2) {
          i1            <- i0
          xfdlist[[i0]] <- xi
        } else {
          #          i1 <- (i1+nVars[ix])
          for(i in 1:nVars[ix]){
            i1  <- i1+1
            xii <- xi
            xii$coefs <- xi$coefs[,,i, drop=FALSE]
            xfdlist[[i1]] <- xii
          }
        }
      }
      else {
        if(is.numeric(xi)) {
          if(nVars[ix]<2) {
            i1 <- i0
            xfdlist[[i0]] <- xi
          } else{
            for(i in 1:nVars[ix]) {
              i1 <- i1+1
              xfdlist[[i1]] <- xi[, i]
            }
          }
        }
      }
    }
  }
  
  ##
  ## 6.  check betalist or set up betalist
  ##
  
  {
    if(inherits(betalist, 'list')) {
      #  betalist is an argument
      if(length(betalist) != k)
        stop('length(betalist) = ', length(betalist),
             ';  must be ', k, ' to match length(xfdlist).')
      betaclass <- sapply(betalist, class)
      oops      <- length(which(betaclass != 'fdPar')) > 0
      if(oops)
        stop('If betalist is a list, all components must have class ',
             'fdPar;  component ', oops[1], ' has class ',
             betaclass[oops[1]])
    } else {
      # betalist must be set up
      betalist <- vector('list', k)
      names(betalist) <- xNames.
      for(i in 1:k) {
        if(is.numeric(xfdlist[[i]])) {
          #. ------------------------------------------------------------------
          #  if xfdlist[[i]] is numeric, use a constant basis
          #. ------------------------------------------------------------------
          if(is.numeric(yvec)) {
            #. use constant basis
            bbasis        <- create.constant.basis(c(0,1))
            bfd           <- fd(1, bbasis)
            betalist[[i]] <- fdPar(bfd)
          } else {
            #. basis is set up using that of dependent variable y
            if(inherits(yvec, 'fd')) {
              #  if 'fd' use the basis of dependent variable
              bbasis  <- yvec$basis
              nbbasis <- bbasis$nbasis
              bfd     <- with(yvec, fd(matrix(0,nbbasis,1), 
                                          bbasis, 
                                          fdnames=fdnames))
              betalist[[i]] <- fdPar(bfd)
            } else {
              #  if 'fdPar' use the basis of dependent variable$fd
              bbasis  <- yvec$fd$basis
              nbbasis <- bbasis$nbasis
              bfd           <- with(yvec, fd(matrix(0,nbbasis,1), 
                                                bbasis, 
                                                fdnames=fdnames))
              betalist[[i]] <- with(yvec, fdPar(bfd, Lfd, lambda, 
                                                estimate, penmat))
            }
          }
        } else {
          #. ------------------------------------------------------------------
          #   if xfdlist[[i]] is fd or fdPar, 
          #.  use basis for the independent variable
          #. ------------------------------------------------------------------
          xfdi <- {
            if(i > 1) xfdlist0[[xL.L0[i-1]]] else xfdlist[[1]]
          }
          if(inherits(xfdi, 'fd')) {
            #. xfdi is a functional data object
            bbasis <- xfdi$basis
            nbasis <- bbasis$nbasis
            coef   <- matrix(0,nbasis,1)
            bfd    <- with(xfdi, fd(coef, bbasis, fdnames=fdnames))
            betalist[[i]] <- fdPar(bfd)
          } else {
            if(inherits(xfdi, 'fdPar')) {
              #. xfdi is a fdPar object
              bbasis <- xfdi$fd$basis
              nbasis <- bbasis$nbasis
              coef   <- matrix(0,nbasis,1)
              bfd    <- with(xfdi$fd, fd(coef, bbasis, fdnames=fdnames))
              betalist[[i]] <- with(xfdi, fdPar(bfd, Lfd, lambda,
                                                estimate, penmat))
            } else {
              stop("Object xfdi is neither fd or fdPar")
            }
          }
        }
      }
    }
  }
  
  ##
  ## 7.  extract or set up weight
  ##
  
  {
    if(is.null(wt))
      wt <- rep(1, ny)
    else {
      if(length(wt) != ny)
        stop('length(wt) must match yvec;  length(wt) = ',
             length(wt), ' != number of yvec observations = ', ny)
      if(any(wt<0))
        stop('Negative weights found;  not allowed.')
    }
  }
  xiEnd   <- cumsum(nVars)
  xiStart <- c(1, xiEnd[-1])
  fRegressList <- list(y=yvec, xfdlist=xfdlist, betalist=betalist, wt=wt)
  
  ##
  ## 8.  either output argument list for fRegress() or invoke itcs
  ##
  
  if(method != "fRegress" && method != "model") 
    stop("Argument is neither 'fRegress' nor 'model'")
  if(method == "model") {
    return(fRegressList)
  } else {
    if(inherits(yvec, 'fd')) do.call('fRegress.fd',     fRegressList)
    else                     do.call('fRegress.double', fRegressList)
  }
}
