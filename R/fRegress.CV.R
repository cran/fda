fRegress.CV <- function(y, xfdlist, betalist, wt=NULL, CVobs=1:N, ...)
{

# FREGRESS.CV computes cross-validated error sum of squares
# for scalar or functional responses. NOTE: ordinary and
# generalized cross validation scores are now returned by fRegress
# when scalar responses are used.

# last modified 29 October 2009 by Jim Ramsay

argList  <- fRegressArgCheck(y, xfdlist, betalist, wt)
yfdPar   <- argList$yfdPar
xfdlist  <- argList$xfdlist
betalist <- argList$betalist
wt       <- argList$wt 

p <- length(xfdlist)
N <- dim(xfdlist[[1]]$coef)[2]

M <- length(CVobs)

if (inherits(yfdPar, "numeric"))  {

    yvec   <- yfdPar
    SSE.CV <- 0
    errfd  <- c()
    for (m in 1:M) {
      i        <- CVobs[m]  
      xfdlisti <- vector("list",p)
      for (j in 1:p) {
        xfdj          <- xfdlist[[j]]
        if (inherits(xfdj, "numeric")) {
          betafdParj <- betalist[[j]]
          betafdj    <- betafdParj$fd
          basisj     <- betafdj$basis
          betarangej <- basisj$rangeval
          conbasisj  <- create.constant.basis(betarangej)
          xfdj       <- fd(matrix(xfdj,1,N), conbasisj)
        }
        basisj <- xfdj$basis
        coefj  <- xfdj$coefs
        if (dim(coefj)[1] == 1) coefj <- matrix(coefj[-i],1,N-1)
        else                    coefj <- as.matrix(coefj[,-i])
        xfdlisti[[j]] <- fd(coefj,basisj)
      }
      yveci         <- yvec[-i]
      wti           <- wt[-i]
      fRegressListi <- fRegress(yveci, xfdlisti, betalist, wti)
      betaestlisti  <- fRegressListi$betaestlist
      yhati <- 0
      for (j in 1:p) {
        betafdParj <- betaestlisti[[j]]
        betafdj    <- betafdParj$fd
        xfdj       <- xfdlist[[j]]
        bbasisj    <- betafdj$basis
        rangej     <- bbasisj$rangeval
        nfine      <- max(101, bbasisj$nbasis*10+1)
        tfine      <- seq(rangej[1], rangej[2], len=nfine)
        delta      <- tfine[2]-tfine[1]
        betavec    <- eval.fd(tfine, betafdj)
        xveci      <- eval.fd(tfine, xfdj[i])
        yhati      <- yhati + delta*(sum(xveci*betavec) -
                                    0.5*( xveci[1]    *betavec[1] +
                                          xveci[nfine]*betavec[nfine] ))
      }
      errfd[i] = yvec[i] - yhati;
      SSE.CV <- SSE.CV + errfd[i]^2
    }
 } else { 
    yfd      <- yfdPar$fd
    SSE.CV   <- 0
    errcoefs <- c()
    for(m in 1:N){
      i <-  CVobs[m]
      txfdlist <- xfdlist           
      for(k in 1:p){
        txfdlist[[k]] <- xfdlist[[k]][-i]
      }
      wti = wt[-i]
      tres <- fRegress(yfd[-i],txfdlist,betalist,wti)
      yhat <- 0                       
      for(k in 1:p){
        yhat <- yhat + xfdlist[[k]][i]*tres$betaestlist[[k]]$fd
      }
      errfdi   <- yfd[i] - yhat
      SSE.CV   <- SSE.CV + inprod(errfdi,errfdi)
      errcoefs <- cbind(errcoefs,errfdi$coefs)
    }
    errfd <- fd(errcoefs,errfdi$basis)
    names(errfd$fdnames)[[3]] <- "Xval Errors"
}
return(list(SSE.CV=SSE.CV,errfd.cv=errfd))
}


