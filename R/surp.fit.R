surp.fit <- function(x, dataList) {
  
  argvals <- dataList$argvals 
  Wbin    <- dataList$Wbin 
  wtvec   <- dataList$wtvec 
  Kmat    <- dataList$Kmat
  Zmat    <- dataList$Zmat 
  Phimat  <- dataList$Phimat 
  n       <- length(argvals)
  M       <- dataList$M
  K       <- dim(Phimat)[2]
  Bmat    <- matrix(x, K, M-1)
  
  logM     <- log(M)
  onewrd   <- all(wtvec == 1)
  Xmat     <- Phimat %*% Bmat %*% t(Zmat)
  expXmat  <- M^Xmat
  sumexpXmat <- as.matrix(apply(expXmat,1,sum))
  Pmat     <- expXmat/(sumexpXmat %*% matrix(1,1,M))
  Smat     <- -Xmat + (log(sumexpXmat) %*% matrix(1,1,M))/logM
  Rmat     <- Wbin - Smat
  vecBmat  <- matrix(Bmat,K*(M-1),1,byrow=TRUE)
  vecRmat  <- matrix(Rmat,n*M,    1,byrow=TRUE)
  vecKmat  <- kronecker(diag(rep(1,M-1)),Kmat)
  fitscale <- 1
  if (!onewrd) {
    vecwtmat <- diag(rep(wtvec,M))
    PENSSE   <- t(vecRmat) %*% diag(as.numeric(wtvec)) %*% vecRmat/fitscale + 
      t(vecBmat) %*% vecKmat %*% vecBmat
  } else {
    PENSSE   <- t(vecRmat) %*% vecRmat/fitscale + t(vecBmat) %*% vecKmat %*% vecBmat
  }
  DvecXmatDvecB <- kronecker(Zmat,Phimat)
  DvecSmatDvecX <- matrix(0,n*M,n*M)
  m2 <- 0
  for (m in 1:M) {
    m1 <- m2 + 1
    m2 <- m2 + n
    m4 <- 0
    for (l in 1:M) {
      m3 <- m4 + 1
      m4 <- m4 + n
      diagPl <- diag(Pmat[,l])
      DvecSmatDvecX[m1:m2,m3:m4] <- diagPl
    }
  }
  DvecSmatDvecX <- DvecSmatDvecX - diag(rep(1,n*M))
  DvecSmatDvecB <- DvecSmatDvecX %*% DvecXmatDvecB
  if (!onewrd) {
    DPENSSE  <- -2*t(DvecSmatDvecB) %*% vecwtmat %*% vecRmat/fitscale + 
      2*vecKmat %*% vecBmat
  } else {
    DPENSSE  <- -2*t(DvecSmatDvecB) %*% vecRmat/fitscale + 2*vecKmat %*% vecBmat
  }
  if (!onewrd) {
    D2PENSSE <-  
      2*((DvecSmatDvecB) %*% vecwtmat %*% DvecSmatDvecB)/fitscale + 2*vecKmat
  } else {
    D2PENSSE <-  2*(t(DvecSmatDvecB) %*% DvecSmatDvecB)/fitscale  + 2*vecKmat
  }
  
  return(list(
    PENSSE   = PENSSE, 
    DPENSSE  = DPENSSE, 
    D2PENSSE = D2PENSSE)
  )
  
}
