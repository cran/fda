AmpPhaseDecomp <- function(xfd, yfd, wfd)
{
#  Computes the amplitude-phase decomposition for a registration.

#  Arguments:
#  XFD  ...  FD object for unregistered functions
#  YFD  ...  FD object for registered functions
#  Wfd  ...  FD object for W functions determining warping functions

#  Returns:
#  MS.amp ... mean square for amplitude variation 
#  MS.pha ... mean square for amplitude variation 
#  RSQR   ... squared correlation measure of prop. phase variation 
#  C      ... constant C

#  Last modified 24 November 2008

xbasis  <- xfd$basis
nxbasis <- xbasis$nbasis
nfine   <- 10*nxbasis
xrng    <- xbasis$rangeval
tfine   <- seq(xrng[1],xrng[2],len=nfine)
delta   <- tfine[2]-tfine[1]

wfine   <- eval.fd(tfine, wfd)
xfine   <- eval.fd(tfine, xfd)
yfine   <- eval.fd(tfine, yfd)
efine   <- exp(wfine)
mufine  <- apply(xfine, 1, mean)
etafine <- apply(yfine, 1, mean)
N       <- dim(xfine)[2]
rfine   <- yfine - outer(etafine,rep(1,N))

intetasqr <- delta*(sum(etafine^2)-0.5*(etafine[1]^2 + etafine[nfine]^2))
intmusqr  <- delta*(sum(mufine^2) -0.5*(mufine[1]^2  + mufine[nfine]^2))

Cnum <- matrix(0,nfine,1)
for (i in 1:nfine) {
    Dhi     <- efine[i,]
    Syi     <- yfine[i,]^2
    Cnum[i] <- cov(Dhi, Syi)
}
intCnum <- delta*(sum(Cnum)-0.5*(Cnum[1]+Cnum[nfine]))
intysqr <- rep(0,N)
intrsqr <- rep(0,N)
for (i in 1:N) {
    intysqr[i] <- delta*(sum(yfine[,i]^2)-0.5*(yfine[1,i]^2 + yfine[nfine,i]^2))
    intrsqr[i] <- delta*(sum(rfine[,i]^2)-0.5*(rfine[1,i]^2 + rfine[nfine,i]^2))
}
C      <- 1 + intCnum/mean(intysqr)
MS.amp <- C*mean(intrsqr)
MS.pha <- C*intetasqr - intmusqr
RSQR   <- MS.pha/(MS.amp+MS.pha)

return(list("MS.amp" = MS.amp, "MS.pha" = MS.pha, "RSQR" = RSQR, "C" = C)) 

}

