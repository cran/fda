predict.fRegress <- function (object, newdata = NULL, se.fit = FALSE, 
          interval = c("none", "confidence", "prediction"), level = 0.95, ...) 
{
  
  #  Last modfied by Jim Ramsay 10 August 2020
  
  #  compute predicted values 
  
  yhatfd <- object$yhatfdobj
  if (is.null(newdata)) {
    pred <- yhatfd
  } else {
    betaestlist <- object$betaestlist
    p <- length(betaestlist)
    for (j in 1:p) {
      if (inherits(betaestlist[[j]], "fdPar")) 
        betaestlist[[j]] <- betaestlist[[j]]$fd
    }
    Nnew <- dim(newdata[[1]]$coefs)[2]
    if (inherits(yhatfd, "fd") || inherits(yhatfd, "fdpar")) {
      for (j in 1:p) {
        xi <- newdata[[j]]
        bi <- betaestlist[[j]]
        if (j == 1) {
          pred <- bi * xi
        }
        else {
          pred <- pred + bi * xi
        }
      }
    }
    else {
      for (j in 1:p) {
        xi <- newdata[[j]]
        bi <- betaestlist[[j]]
        if (j == 1) {
          pred <- inprod(xi, bi)
        }
        else {
          pred <- pred + inprod(xi, bi)
        }
      }
    }
  }
  
  #  check that standard errors of predicted values are required
  
  int <- match.arg(interval)
  need.se <- (se.fit || (int != "none"))
  if (!need.se) {
    return(pred)
  }
  else {
    
    #  compute variance-covariance matrix over plotting grid
    
    Bvar = object$Bvar
    if (is.null(Bvar)) 
      stop(paste("Standard error for predict object cannot be computed", 
                 " without preliminary use of function fRegress.stderr()."))
    ncoef <- 0
    for (j in 1:p) {
      betafdj <- betaestlist[[j]]
      ncoefj <- betafdj$basis$nbasis
      ncoef <- ncoef + ncoefj
    }
    
    if (inherits(yhatfd, "fdPar") || inherits(yhatfd, "fd")) {

      #  functional dependent variable case
      
      nplot <- 101
      rangeval <- yhatfd$basis$rangeval
      tplot <- seq(rangeval[1], rangeval[2], len = nplot)
      YhatStderr <- matrix(0, nplot, N)
      B2YhatList <- vector("list", p)
      for (iplot in 1:nplot) {
        YhatVari <- matrix(0, N, N)
        tval <- tplot[iplot]
        for (j in 1:p) {
          Zmat <- eval.fd(tval, newdata[[j]])
          betabasisj <- betaestlist[[j]]$basis
          PsiMatj <- eval.basis(tval, betabasisj)
          B2YhatMapij <- t(Zmat) %*% PsiMatj
          B2YhatList[[j]] <- B2YhatMapij
        }
        m2j <- 0
        for (j in 1:p) {
          m1j <- m2j + 1
          m2j <- m2j + betaestlist[[j]]$basis$nbasis
          B2YhatMapij <- B2YhatList[[j]]
          m2k <- 0
          for (k in 1:p) {
            m1k <- m2k + 1
            m2k <- m2k + betaestlist[[k]]$basis$nbasis
            B2YhatMapik <- B2YhatList[[k]]
            YhatVari <- YhatVari + B2YhatMapij %*% Bvar[m1j:m2j, 
                                                        m1k:m2k] %*% t(B2YhatMapik)
          }
        }
        YhatStderr[iplot, ] <- matrix(sqrt(diag(YhatVari)), 
                                      1, N)
      }
    }
    else {
      
      #  scale dependent variable case
      
      ymat <- as.matrix(yhatfd)
      N <- dim(ymat)[1]
      B2YhatList <- vector("list", p)
      YhatVari <- matrix(0, N, N)
      for (j in 1:p) {
        betabasisj <- betaestlist[[j]]$basis
        Xfdj <- newdata[[j]]
        B2YhatMapij <- inprod(Xfdj, betabasisj)
        B2YhatList[[j]] <- B2YhatMapij
      }
      m2j <- 0
      for (j in 1:p) {
        m1j <- m2j + 1
        m2j <- m2j + betaestlist[[j]]$basis$nbasis
        B2YhatMapij <- B2YhatList[[j]]
        m2k <- 0
        for (k in 1:p) {
          m1k <- m2k + 1
          m2k <- m2k + betaestlist[[k]]$basis$nbasis
          B2YhatMapik <- B2YhatList[[k]]
          YhatVari <- YhatVari + B2YhatMapij %*% Bvar[m1j:m2j, 
                                                      m1k:m2k] %*% t(B2YhatMapik)
        }
      }
      YhatStderr <- matrix(sqrt(diag(YhatVari)), N, 1)
    }
    return(list(pred = pred, YhatStderr = YhatStderr))
  }
}
