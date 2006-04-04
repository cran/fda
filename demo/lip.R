#  -----------------------------------------------------------------------
#                       Lip Movement Data
#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  These are rather simple data, involving the movement of the lower lip
#  while saying "bob".  There are 20 replications and 51 sampling points.
#  The data are used to illustrate two techniques:  landmark registration
#  and principal differental analysis.  
#  Principal differential analysis estimates a linear differential equation
#  that can be used to describe not only the observed curves, but also a 
#  certain number of their derivatives.  
#  For a rather more elaborate example of principal differential analysis, 
#  see the handwriting data.
#  -----------------------------------------------------------------------

#  Last modified 21 March 2006

#  ----------------  input the data  ------------------------

liptime  <- seq(0,1,.02)
liprange <- c(0,1)

#  -------------  create the fd object -----------------
#       use 31 order 6 splines so we can look at acceleration

nbasis <- 51
norder <- 6
lipbasis <- create.bspline.basis(liprange, nbasis, norder)

#  ------------  apply some light smoothing to this object  -------

Lfdobj   <- int2Lfd(4)
lambda   <- 1e-12
lipfdPar <- fdPar(lipbasis, Lfdobj, lambda)

lipfd <- smooth.basis(liptime, lip, lipfdPar)$fd
names(lipfd$fdnames) = c("Normalized time", "Replications", "mm")

#  set up plotting arrangements for one and two panel displays allowing
#  for larger fonts

#  ---------  plot the functions and their accelerations  -----

par(mfrow=c(2,1), mar=c(5,5,4,2), pty="m", ask=FALSE)
plot(lipfd,        main="Lip Position", cex=1.2)
plot(lipfd, Lfd=2, ylab="mm/sec/sec", main="Lip Acceleration", cex=1.2)

#  -----------------------------------------------------------------------
#       Register the data using the two landmarks defined by the minimum
#        and the right elbow.
#       Manually identify these points in each curve
#  -----------------------------------------------------------------------

nmarks <- 2

lipmat   <- eval.fd(liptime,lipfd)

lipmeanfd <- mean.fd(lipfd)

par(mfrow=c(1,1),pty="m")
lipmarks <- matrix(0,20,nmarks)
index <- 1:20
for (i in index) {
  plot(liptime, lipmat[,i], xlab="", ylab="", main=paste("Curve",i))
  indexi <- identify(liptime, lipmat[,i], n=nmarks)
  lipmarks[i,] <- liptime[indexi]
}

lipmeanmarks <- apply(lipmarks,2,mean)

#  -------------   register the curves  --------------------

#  First create a basis object for the warping function
#  it has order 4 (piecewise cubic) and two interior knots
#  positioned at the mean landmark values since
#  NBASIS = NORDER + # interior knots

wnbasis <- 6
wnorder <- 4
wbreaks <- c(0,lipmeanmarks,1)
warpbasis <- create.bspline.basis(liprange, wnbasis, wnorder, wbreaks);
WfdPar    <- fdPar(fd(matrix(0,wnbasis,1), warpbasis), 2, 1e-4)

lipreglist <- landmarkreg(lipfd, lipmarks, lipmeanmarks, WfdPar)

lipregfd   <- lipreglist$regfd
lipwarpfd  <- lipreglist$warpfd

#  plot unregistered and registered curves

par(mfrow=c(1,2), pty="s")

plot(lipfd, main="Unregistered")
lines.fd(lipmeanfd, lty=2)
abline(v=lipmeanmarks,lty=2)

plot(lipregfd, main="Registered")
lines.fd(lipmeanfd, lty=2)
abline(v=lipmeanmarks,lty=2)

#  plot warping functions and deformations

par(mfrow=c(1,2), pty="s")
plot(lipwarpfd, href=FALSE, main="Warping Functions")
abline(0,1,lty=2)
hmat <- eval.fd(liptime, lipwarpfd)
defmat <- hmat - outer(liptime,rep(1,20))
matplot(liptime,defmat,type="l",lty=1, 
        xlab="Normalized time", ylab="Warped Normalized time",
        main="Deformation Functions")
abline(h=0,lty=2)

#  ------------  carry out a pca and plot results  -------------------

lambda    <- 1e-6
pcafdPar  <- fdPar(lipbasis, 2, lambda)
lippca.fd <- pca.fd(lipfd, nharm=3, pcafdPar)

par(mfrow=c(1,1),pty="m")
plot.pca.fd(lippca.fd)

lipeigvals <- lippca.fd[[2]]
plot(1:19, log10(lipeigvals[1:19]), type="b",
     xlab="Eigenvalue Number", ylab="", main="Log10 Eigenvalues")

#  ---------------------------------------------------------------------
#                    Principal differential analysis  
#  ---------------------------------------------------------------------

#  set up a second order linear differnetial equation solution

pdabasisfd <- create.bspline.basis(liprange, nbasis=21)
betafdPar  <- fdPar(pdabasisfd)

#  set up list of functional parameter objects for weight fns.

bwtlist = vector("list", 2)
bwtlist[[1]] <- betafdPar
bwtlist[[2]] <- betafdPar

xfdlist <- list(lipfd)

pdaList <- pda.fd(xfdlist, bwtlist)

#  plot weight functions

bwtestlist <- pdaList$bwtlist

par(mfrow=c(2,1),pty="m")
for (j in 1:2) {
	bfdParj <- bwtestlist[[j]]
	plot(bfdParj$fd)
}

#  compute forcing functions

Lfdest <- Lfd(2, bwtestlist)

force        <- eval.fd(liptime, lipfd, Lfdest)
lipaccel     <- eval.fd(liptime, lipfd, 2)
lipmeanaccel <- apply(lipaccel, 1, mean)

par(mfrow=c(1,1),ask=FALSE)
yrange <- c(min(min(lipmeanaccel),min(force)),
            max(max(lipmeanaccel),max(force)))
matplot(liptime, force, type="l", lty=1, ylim=yrange)
lines(liptime, lipmeanaccel, lty=4, lwd=2)

#  plot the mean forcing function along with second deriv.

forcemean <- apply(force, 1, mean)

plot(liptime, forcemean, type="l", lty=1, ylim=yrange)
lines(liptime, lipmeanaccel, lty=4)

#  solve equation

result <- odesolv(bwtestlist)
xp <- result[[1]]
yp <- result[[2]]

#  plot the two solutions

par(mfrow=c(2,1),pty="m")
pltrng <- c(min(yp[1,,]), max(yp[1,,]))
matplot(xp,t(yp[1,,]), type="l", lty=1, ylim=pltrng, main="Function")
abline(h=0, lty=2)
pltrng <- c(min(yp[2,,]), max(yp[2,,]))
matplot(xp,t(yp[2,,]), type="l", lty=1, ylim=pltrng, main="Derivative")
abline(h=0, lty=2)

#  plot fit to each curve

lipmat   <- eval.fd(liptime, lipfd)
D2lipmat <- eval.fd(liptime, lipfd, 2)

umat <- matrix(0,length(liptime),2)
umat[,1] <- approx(xp, t(yp[1,1,]), liptime)$y
umat[,2] <- approx(xp, t(yp[1,2,]), liptime)$y

par(mfrow=c(1,2),pty="s",ask=TRUE)
index <- 1:20
for (i in index) {
    plot(liptime, force[,i], type="l",
          ylim=c(-1000,1000), xlab="Normalized Time", ylab="", 
          main=paste("Record",i,"Forcing Fn."))
    lines(liptime, D2lipmat[,i],lty=4)
    abline(h=0,lty=2)
    xhat <- lipmat[,i] - lsfit(umat, lipmat[,i], int=FALSE)$residual
    matplot(liptime, cbind(xhat, lipmat[,i]), type="l", lty=c(1,2),
          xlab="Normalized Time", ylab="", main="Function")
}

