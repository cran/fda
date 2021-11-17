fRegress.double <- function(y, xfdlist, betalist, wt=NULL,
                             y2cMap=NULL, SigmaE=NULL, returnMatrix=FALSE, ...)
{
  
  #  FREGRESS.DOUBLE  Fits a scalar dependent variable using the concurrent
  #                    functional regression model using inner products
  #                    of functional covariates and functional regression
  #                    functions.
  #
  #  Arguments:
  #  Y        ... an object for the dependent variable,
  #               which  be a numeric vector
  #  XFDLIST  ... a list object of length p with each list
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
  #               with .
  #  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
  #               from a call to function BsplineS.  See this function for
  #               enabling this option.
  #
  #  Returns LIST  ... A list containing seven members with names:
  #    yfdobj      ... first  argument of 
  #    xfdlist     ... second argument of 
  #    betalist    ... third  argument of 
  #    betaestlist ... estimated regression functions
  #    yhatfdobj   ... functional data object containing fitted functions
  #    Cmatinv     ... inverse of the coefficient matrix, needed for
  #                    function .STDERR that computes standard errors
  #    wt          ... weights for observations
  #    df          ... degrees of freedom for fit
  #  This list object is converted to a class with the name ""
  #  function predict. is an example of a method that can be called simply
  #  as predict(List).  In this call List can be any object of the
  #  "".
  
  # Last modified 16 December 2020 by Jim Ramsay
  
  #  check Y and compute sample size N
  
  if (!inherits(y, "numeric")) stop("Y is not a numeric vector.")
    
  #  ----------------------------------------------------------------
  #                   yfdobj is scalar or multivariate
  #  ----------------------------------------------------------------
    
  #  As of 2020, if yfd is an fdPar object, it is converted to an fd object.
  #  The added structure of the fdPar class is not used in any of the fRegress codes.
  #  The older versions of fda package used yfdPar as the name for the first member.

  arglist <- fRegressArgCheck(y, xfdlist, betalist, wt)
  
  yfdobj   <- arglist$yfd
  xfdlist  <- arglist$xfdlist
  betalist <- arglist$betalist
  wt       <- arglist$wt
  
  ymat <- as.matrix(y)
  N    <- dim(ymat)[1]
  p    <- length(xfdlist)
    
  Zmat  <- NULL
  Rmat  <- NULL
  pjvec <- rep(0,p)
  ncoef <- 0
  for (j in 1:p) {
    xfdj       <- xfdlist[[j]]
    if (!inherits(xfdj, "fd")) {
      stop(paste("Independent variable",j,"is not of class fd."))
    }
    xcoef      <- xfdj$coefs
    xbasis     <- xfdj$basis
    betafdParj <- betalist[[j]]
    bbasis     <- betafdParj$fd$basis
    bnbasis    <- bbasis$nbasis
    pjvec[j]   <- bnbasis
    Jpsithetaj <- inprod(xbasis,bbasis)
    Zmat       <- cbind(Zmat,crossprod(xcoef,Jpsithetaj))
    if (betafdParj$estimate) {
      lambdaj    <- betafdParj$lambda
      if (lambdaj > 0) {
        Lfdj  <- betafdParj$Lfd
        Rmatj <- lambdaj*eval.penalty(bbasis,Lfdj)
      } else {
        Rmatj <- matrix(0,bnbasis,bnbasis)
      }
      if (ncoef > 0) {
        zeromat <- matrix(0,ncoef,bnbasis)
        Rmat    <- rbind(cbind(Rmat,       zeromat),
                         cbind(t(zeromat), Rmatj))
      } else {
        Rmat  <- Rmatj
        ncoef <- ncoef + bnbasis
      }
    }
  }
  
  #  -----------------------------------------------------------
  #          set up the linear equations for the solution
  #  -----------------------------------------------------------
  
  #  solve for coefficients defining BETA
  
  if (any(wt != 1)) {
    rtwt   <- sqrt(wt)
    Zmatwt <- Zmat*rtwt
    ymatwt <- ymat*rtwt
    Cmat   <- t(Zmatwt) %*% Zmatwt + Rmat
    Dmat   <- t(Zmatwt) %*% ymatwt
  } else {
    Cmat <- t(Zmat) %*% Zmat + Rmat
    Dmat <- t(Zmat) %*% ymat
  }
  
  eigchk(Cmat)
  
  Cmatinv  <- solve(Cmat)
  
  betacoef <- Cmatinv %*% Dmat
  
#    df <- sum(diag(Zmat %*% Cmatinv %*% t(Zmat)))

    hatvals = diag(Zmat %*% Cmatinv %*% t(Zmat))
    df <- sum(hatvals)
  
  #  set up fdPar object for BETAESTFDPAR
  
  betaestlist <- betalist
  mj2 <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    betafdj    <- betafdParj$fd
    ncoefj     <- betafdj$basis$nbasis
    mj1        <- mj2 + 1
    mj2        <- mj2 + ncoefj
    indexj     <- mj1:mj2
    betacoefj        <- betacoef[indexj]
    betaestfdj       <- betafdj
    betaestfdj$coefs <- as.matrix(betacoefj)
    betaestfdParj    <- betafdParj
    betaestfdParj$fd <- betaestfdj
    betaestlist[[j]] <- betaestfdParj
  }
  
  #  set up fd object for predicted values
  
  yhatmat <- matrix(0,N,1)
  for (j in 1:p) {
    xfdj <- xfdlist[[j]]
    if (inherits(xfdj, "fd")) {
      xbasis     <- xfdj$basis
      xnbasis    <- xbasis$nbasis
      xrng       <- xbasis$rangeval
      nfine      <- max(501,10*xnbasis+1)
      tfine      <- seq(xrng[1], xrng[2], len=nfine)
      deltat     <- tfine[2]-tfine[1]
      xmat       <- eval.fd(tfine, xfdj, 0, returnMatrix)
      betafdParj <- betaestlist[[j]]
      betafdj    <- betafdParj$fd
      betamat    <- eval.fd(tfine, betafdj, 0, returnMatrix)
      fitj       <- deltat*(crossprod(xmat,betamat) -
                              0.5*(outer(xmat[1,    ],betamat[1,    ]) +
                                     outer(xmat[nfine,],betamat[nfine,])))
      yhatmat    <- yhatmat + fitj
    } else{
      betaestfdParj <- betaestlist[[j]]
      betavecj      <- betaestfdParj$fd$coefs
      yhatmat       <- yhatmat + xfdj %*% t(betavecj)
    }
  }
  yhatfdobj <- yhatmat
  
    # Calculate OCV and GCV scores

    OCV = sum( (ymat-yhatmat)^2/(1-hatvals)^2 )
    GCV = sum( (ymat-yhatmat)^2 )/( (sum(1-hatvals))^2 )
  
  #  -----------------------------------------------------------------------
  #        Compute pointwise standard errors of regression coefficients
  #               if both y2cMap and SigmaE are supplied.
  #  -----------------------------------------------------------------------
  
  
  if (!(is.null(y2cMap) || is.null(SigmaE))) {
    
    #  check dimensions of y2cMap and SigmaE
    
    y2cdim <- dim(y2cMap)
    if (y2cdim[2] != dim(SigmaE)[1])  stop(
      "Dimensions of Y2CMAP not correct.")
    
    
    #  compute linear mapping c2bMap takinging coefficients for
    #  response into coefficients for regression functions
    
    c2bMap <- Cmatinv %*% t(Zmat)
    y2bmap <- c2bMap
    bvar   <- y2bmap %*% as.matrix(SigmaE) %*% t(y2bmap)
    betastderrlist <- vector("list",p)
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj     <- betabasisj$nbasis
      mj1        <- mj2 + 1
      mj2        <- mj2 + ncoefj
      indexj     <- mj1:mj2
      bvarj      <- bvar[indexj,indexj]
      betarng    <- betabasisj$rangeval
      nfine      <- max(c(501,10*ncoefj+1))
      tfine      <- seq(betarng[1], betarng[2], len=nfine)
      bbasismat  <- eval.basis(tfine, betabasisj, 0, returnMatrix)
      bstderrj   <- sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
      bstderrfdj <- smooth.basis(tfine, bstderrj, betabasisj)$fd
      betastderrlist[[j]] <- bstderrfdj
    }
  } else {
    betastderrlist = NULL
    bvar           = NULL
    c2bMap         = NULL
  }
  
  #  -----------------------------------------------------------------------
  #                  Set up output list object
  #  -----------------------------------------------------------------------
  
  fRegressList <-
    list(yfdobj         = y,
         xfdlist        = xfdlist,
         betalist       = betalist,
         betaestlist    = betaestlist,
         yhatfdobj      = yhatfdobj,
         Cmat           = Cmat,
         Dmat           = Dmat,
         Cmatinv        = Cmatinv,
         wt             = wt,
         df             = df,
		 GCV			= GCV,
		 OCV			= OCV,
         y2cMap         = y2cMap,
         SigmaE         = SigmaE,
         betastderrlist = betastderrlist,
         bvar           = bvar,
         c2bMap         = c2bMap)
  
  class(fRegressList) <- "fRegress"
  
  return(fRegressList)
  
}