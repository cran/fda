fRegress.stderr <- function (y, y2cMap, SigmaE, returnMatrix = FALSE, ...) 
{
  
  #  FREGRESS.STDERR  computes standard error estimates for regression
  #       coefficient functions estimated by function FREGRESS.
  #
  #  Arguments:
  #
  #  Y            ... a list object produced by function FREGRESS with class name
  #                   'fRegress".  
  #                   This is indicated by Y in the arguments since R syntax
  #                   requires all of tghe fRegress family of functions to
  #                   use this notation.
  #  Y2CMAP       ... the matrix mapping from the vector of observed values
  #                   to the coefficients for the dependent variable.
  #                   This is output by function SMOOTH_BASIS.  If this is
  #                   supplied, confidence limits are computed, otherwise not.
  #  SIGMAE       ... Estimate of the covariances among the residuals.  This
  #                   can only be estimated after a preliminary analysis
  #                   with FREGRESS.
  #  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
  #               from a call to function BsplineS.  See this function for
  #               enabling this option.
  #
  #  Returns:
  #
  #  BETASTDERRLIST ... a list object, each list containing a fdPar object
  #                     for the standard error of a regression function.
  #  BVAR           ... the symmetric matrix of sampling variances and
  #                     covariances for the matrix of regression coefficients
  #                     for the regression functions.  These are stored
  #                     column-wise in defining BVARIANCE.
  #  C2BMAP         ... the matrix mapping from response variable coefficients
  #                     to coefficients for regression coefficients
  
  #  Last modified 10 August 2020
  
  #  retrieve objects from y
  
  yfdPar         <- y$yfdPar
  xfdlist        <- y$xfdlist
  betalist       <- y$betalist
  betaestlist    <- y$betaestlist
  yhatfdobj      <- y$yhatfdobj
  Cmat           <- y$Cmat
  Dmat           <- y$Dmat
  Cmatinv        <- y$Cmatinv
  wt             <- y$wt
  df             <- y$df
  betastderrlist <- y$betastderrlist
  YhatStderr     <- y$YhatStderr
  Bvar           <- y$Bvar
  c2bMap         <- y$c2bMap
  
  #  get number of independent variables
  
  p <- length(xfdlist)

  #  compute number of coefficients
  
  ncoef <- 0
  for (j in 1:p) {
    betaParfdj <- betalist[[j]]
    ncoefj <- betaParfdj$fd$basis$nbasis
    ncoef <- ncoef + ncoefj
  }
  
  if (inherits(yfdPar, "fdPar") || inherits(yfdPar, "fd")) {
    
    #  ----------------------------------------------------------------
    #           YFDPAR is functional for a functional parameter
    #  ----------------------------------------------------------------
    
    if (inherits(yfdPar, "fd")) 
      yfdPar <- fdPar(yfdPar)
    
    #  get number of replications and basis information for YFDPAR
    
    yfd <- yfdPar$fd
    N <- dim(yfd$coefs)[2]
    ybasisobj <- yfdPar$fd$basis
    rangeval <- ybasisobj$rangeval
    ynbasis <- ybasisobj$nbasis
    ninteg <- max(501, 10 * ynbasis + 1)
    tinteg <- seq(rangeval[1], rangeval[2], len = ninteg)
    deltat <- tinteg[2] - tinteg[1]
    
    ybasismat <- eval.basis(tinteg, ybasisobj, 0, returnMatrix)
    
    #  compute BASISPRODMAT
    
    basisprodmat <- matrix(0, ncoef, ynbasis * N)
    
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj <- betabasisj$nbasis
      bbasismatj <- eval.basis(tinteg, betabasisj, 0, returnMatrix)
      xfdj <- xfdlist[[j]]
      tempj <- eval.fd(tinteg, xfdj, 0, returnMatrix)
      mj1 <- mj2 + 1
      mj2 <- mj2 + ncoefj
      indexj <- mj1:mj2
      mk2 <- 0
      for (k in 1:ynbasis) {
        mk1 <- mk2 + 1
        mk2 <- mk2 + N
        indexk <- mk1:mk2
        tempk <- bbasismatj * ybasismat[, k]
        basisprodmat[indexj, indexk] <- deltat * crossprod(tempk, 
                                                           tempj)
      }
    }
    
    #  check dimensions of Y2CMAP
    
    y2cdim <- dim(y2cMap)
    if (y2cdim[1] != ynbasis || y2cdim[2] != dim(SigmaE)[1]) 
      stop("Dimensions of Y2CMAP not correct.")
    
    #  compute variances of regression coefficient function values
    
    Varc <- y2cMap %*% SigmaE %*% t(y2cMap)
    CVar <- kronecker(Varc, diag(rep(1, N)))
    C2BMap <- Cmatinv %*% basisprodmat
    Bvar <- C2BMap %*% CVar %*% t(C2BMap)
    nplot <- max(51, 10 * ynbasis + 1)
    tplot <- seq(rangeval[1], rangeval[2], len = nplot)
    betastderrlist <- vector("list", p)
    PsiMatlist <- vector("list", p)
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj <- betabasisj$nbasis
      mj1 <- mj2 + 1
      mj2 <- mj2 + ncoefj
      indexj <- mj1:mj2
      bbasismat <- eval.basis(tplot, betabasisj, 0, returnMatrix)
      PsiMatlist <- bbasismat
      bvarj <- Bvar[indexj, indexj]
      bstderrj <- sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
      bstderrfdj <- smooth.basis(tplot, bstderrj, betabasisj)$fd
      betastderrlist[[j]] <- bstderrfdj
    }
    
    #  compute estimated variance-covariance matrix over plotting grid
    #  of fitted values
    
    YhatStderr <- matrix(0, nplot, N)
    B2YhatList <- vector("list", p)
    for (iplot in 1:nplot) {
      YhatVari <- matrix(0, N, N)
      tval <- tplot[iplot]
      for (j in 1:p) {
        Zmat <- eval.fd(tval, xfdlist[[j]])
        betabasisj <- betalist[[j]]$fd$basis
        PsiMatj <- eval.basis(tval, betabasisj)
        B2YhatMapij <- t(Zmat) %*% PsiMatj
        B2YhatList[[j]] <- B2YhatMapij
      }
      m2j <- 0
      for (j in 1:p) {
        m1j <- m2j + 1
        m2j <- m2j + betalist[[j]]$fd$basis$nbasis
        B2YhatMapij <- B2YhatList[[j]]
        m2k <- 0
        for (k in 1:p) {
          m1k <- m2k + 1
          m2k <- m2k + betalist[[k]]$fd$basis$nbasis
          B2YhatMapik <- B2YhatList[[k]]
          YhatVari <- YhatVari + 
            B2YhatMapij %*% Bvar[m1j:m2j,m1k:m2k] %*% t(B2YhatMapik)
        }
      }
      YhatStderr[iplot, ] <- matrix(sqrt(diag(YhatVari)), 1, N)
    }
  }
  
  else {
    
    #  ----------------------------------------------------------------
    #                   YFDPAR is scalar or multivariate
    #  ----------------------------------------------------------------
    
    ymat <- as.matrix(yfdPar)
    N <- dim(ymat)[1]
    Zmat <- NULL
    for (j in 1:p) {
      xfdj <- xfdlist[[j]]
      if (inherits(xfdj, "fd")) {
        xcoef <- xfdj$coefs
        xbasis <- xfdj$basis
        betafdParj <- betalist[[j]]
        bbasis <- betafdParj$fd$basis
        Jpsithetaj <- inprod(xbasis, bbasis)
        Zmat <- cbind(Zmat, t(xcoef) %*% Jpsithetaj)
      }
      else if (inherits(xfdj, "numeric")) {
        Zmatj <- xfdj
        Zmat <- cbind(Zmat, Zmatj)
      }
    }
    
    #  compute linear mapping c2bMap takinging coefficients for
    #  response into coefficients for regression functions
    
    c2bMap <- Cmatinv %*% t(Zmat)
    y2bmap <- c2bMap
    Bvar <- y2bmap %*% as.matrix(SigmaE) %*% t(y2bmap)
    betastderrlist <- vector("list", p)
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj <- betabasisj$nbasis
      mj1 <- mj2 + 1
      mj2 <- mj2 + ncoefj
      indexj <- mj1:mj2
      bvarj <- Bvar[indexj, indexj]
      xfdj <- xfdlist[[j]]
      if (inherits(xfdj, "fd")) {
        betarng <- betabasisj$rangeval
        ninteg <- max(c(501, 10 * ncoefj + 1))
        tinteg <- seq(betarng[1], betarng[2], len = ninteg)
        bbasismat <- eval.basis(tinteg, betabasisj, 0, 
                                returnMatrix)
        bstderrj <- sqrt(diag(bbasismat %*% bvarj %*% 
                                t(bbasismat)))
        bstderrfdj <- smooth.basis(tinteg, bstderrj, 
                                   betabasisj)$fd
      }
      else {
        bsterrj <- sqrt(diag(bvarj))
        onebasis <- create.constant.basis(betabasisj$rangeval)
        bstderrfdj <- fd(t(bstderrj), onebasis)
      }
      betastderrlist[[j]] <- bstderrfdj
    }
    
    #  compute estimated variance-covariance matrix for fitted values
    
    B2YhatList <- vector("list", p)
    YhatVari <- matrix(0, N, N)
    for (j in 1:p) {
      betabasisj <- betalist[[j]]$fd$basis
      Xfdj <- xfdlist[[j]]
      B2YhatMapij <- inprod(Xfdj, betabasisj)
      B2YhatList[[j]] <- B2YhatMapij
    }
    m2j <- 0
    for (j in 1:p) {
      m1j <- m2j + 1
      m2j <- m2j + betalist[[j]]$fd$basis$nbasis
      B2YhatMapij <- B2YhatList[[j]]
      m2k <- 0
      for (k in 1:p) {
        m1k <- m2k + 1
        m2k <- m2k + betalist[[k]]$fd$basis$nbasis
        B2YhatMapik <- B2YhatList[[k]]
        YhatVari <- YhatVari + B2YhatMapij %*% Bvar[m1j:m2j, 
                                                    m1k:m2k] %*% t(B2YhatMapik)
      }
    }
    YhatStderr <- matrix(sqrt(diag(YhatVari)), N, 1)
  }
  
  #  return output object of class fRegress
  
  fRegressList <- list(yfdPar = y$yfdPar, xfdlist = y$xfdlist, 
                       betalist = y$betalist, betaestlist = y$betaestlist, yhatfdobj = y$yhatfdobj, 
                       Cmat = y$Cmat, Dmat = y$Dmat, Cmatinv = y$Cmatinv, wt = y$wt, 
                       df = y$df, y2cMap = y2cMap, SigmaE = SigmaE, betastderrlist = betastderrlist, 
                       YhatStderr = YhatStderr, Bvar = Bvar, c2bMap = c2bMap)
  
  class(fRegressList) = "fRegress"
  return(fRegressList)
}