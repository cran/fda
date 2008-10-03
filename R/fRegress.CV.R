fRegress.CV <- function(y, xfdlist, betalist, ...){

# FREGRESS.CV computes cross-validated error sum of squares
# only for scalar dependent variable

# last modified 2008.12.19 by Spencer
# Previously modified 15 December 2008 by Jim
  yvec <- y
  if (!inherits(yvec, "numeric"))
    stop("Dependent variable is not scalar.")

  N <- length(yvec)
  p <- length(xfdlist)
  betafdPar <- betalist[[2]]
  betarange <- betafdPar$fd$basis$rangeval
  SSE.CV    <- 0
  for (i in 1:N) {
    xfdlisti <- vector("list",p)
    for (j in 1:p) {
        xfdj   <- xfdlist[[j]]
        if (inherits(xfdj, "numeric")) {
            xfdj <- fd(matrix(xfdj,1,N), create.constant.basis(betarange))
        }
        basisj <- xfdj$basis
        coefj  <- xfdj$coefs
        if (dim(coefj)[1] == 1) coefj <- matrix(coefj[-i],1,N-1)
        else                    coefj <- as.matrix(coefj[,-i])
        xfdlisti[[j]] <- fd(coefj,basisj)
    }
    yveci         <- yvec[-i]
    fRegressListi <- fRegress(yveci, xfdlisti, betalist)
    betaestlisti  <- fRegressListi$betaestlist
    yhati <- 0
    for (j in 1:p) {
        betafdParj <- betaestlisti[[j]]
        betafdj    <- betafdParj$fd
        xfdj       <- xfdlist[[j]]
        if (inherits(xfdj, "numeric")) {
            xfdj <- fd(matrix(xfdj,1,N), create.constant.basis(betarange))
        }
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
    SSE.CV <- SSE.CV + (yvec[i] - yhati)^2
  }

  return(SSE.CV)

}

