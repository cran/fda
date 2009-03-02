###
###
### Ramsey, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 11  Functional Models and Dynamics
###
library(fda)

##
## Section 11.1  Introduction to Dynamics
##

# section 11.1.2 Interpreting Second Order Linear Dynamics
#  Figure 11.1

x = seq(-1, 1, length=201)
x2 <- x^2
op = par(xpd=TRUE)
plot(x, x2, xlim=c(-1,1), axes=FALSE, lwd=2, type='l', xlab='', ylab='')
axis(1, 0)
axis(2, 0, las=1)
lines(c(-1, 1), c(0, 0), lty='dashed')
lines(c(0, 0), c(0, 1), lty='dashed')
text(-1.2, .5, expression(beta[0]), cex=2 )
text(0, -0.2, expression(beta[1]), cex=2)
text(-0.4, 0.8, "Increasing\nOscillations", cex=2)
text(0.4, 0.8, "Decreasing\nOscillations", cex=2)
text(-0.7, 0.005, 'Exponential\nGrowth', cex=2)
text(0.7, 0.005, 'Exponential\nDecay', cex=2)
arrows(-0.3, 0.2, -0.4, 0.16, length=0.1, lwd=2)
text(-0.29, 0.2, "d=0", adj=0, cex=2)
par(op)

##
## Section 11.2 Principal Differential Analysis for Linear Dynamics
##
#  (no computations in this section)

##
## Section 11.3 Principal Differential Analysis of the Lip Data
##
# Figure 11.2
matplot(liptime, lip, type = 'l',
        xlab='Normalized Time', ylab='lip position (mm)')

# Prep for Figure 11.3

lipfd = smooth.basisPar(liptime, lip, 6, Lfdobj=int2Lfd(4),
                         lambda=1e-12)$fd
names(lipfd$fdnames) = c("time(seconds)", "replications", "mm")

lipbasis  = lipfd$basis
bwtlist   = list(fdPar(lipbasis,2,0),fdPar(lipbasis,2,0))

xfdlist1   = list(lipfd)

pdaList   = pda.fd(xfdlist1, bwtlist)
bwtestlist= pdaList$bwtlist

bwtestlist[[1]]$fd$fdnames = list('time','rep','beta0')
bwtestlist[[2]]$fd$fdnames = list('time','rep','beta1')

dfd        = 0.25*bwtestlist[[2]]$fd^2 - bwtestlist[[1]]$fd
dfd$fdnames= list('time','rep','discriminant')

# The top two panels of Figure 11.3

plot.pda.fd(pdaList,whichdim=3,cex.axis=1.5,cex.lab=1.5,lwd=2)
plot(pdaList,whichdim=3,cex.axis=1.5,cex.lab=1.5,lwd=2)

# Figure 11.3

op = par(mfrow=c(3,1))
plot(bwtestlist[[1]]$fd,cex.lab=1.5,cex.axis=1.5,lwd=2,main="beta 0")
plot(bwtestlist[[2]]$fd,cex.lab=1.5,cex.axis=1.5,lty=2,lwd=2,main="beta 1")
plot(dfd,cex.lab=1.5,cex.axis=1.5,lwd=2,main="discriminant")
par(op)

# Figure 11.4

pda.overlay(pdaList)

##
## Section 11.4 PDA of the Handwriting Data
##
fdabasis= create.bspline.basis(norder=7, breaks=handwritTime)

fdafd0 = fd(array(0, c(fdabasis$nbasis, dim(handwrit)[-1])), fdabasis)
lambda = 1e8
fdaPar = fdPar(fdafd0, 5, lambda)

smoothList= smooth.basis(handwritTime, handwrit, fdaPar)
fdafd     = smoothList$fd
df        = smoothList$df
gcv       = smoothList$gcv

#  Add suitable names for the dimensions of the data.

fdafd$fdnames[[1]] <- "Milliseconds"
fdafd$fdnames[[2]] <- "Replications"
fdafd$fdnames[[3]] <- "Metres"

xfdlist = list(fdafd[,1],fdafd[,2])

pdaPar    = fdPar(fdabasis,2,1)
pdaParlist= list(pdaPar, pdaPar)
bwtlist   = list( list(pdaParlist,pdaParlist),
                  list(pdaParlist,pdaParlist) )

# This can take some time to compute, so time it:
pdaTime   = system.time({
pdaList   = pda.fd(xfdlist, bwtlist)} )
pdaTime/60
# 20 minutes on a dual core;  2009.03.21

# Figure 11.5

eigenres = eigen.pda(pdaList)

##
## Section 11.5 Registration and PDA
##

WfdPar = fdPar(lipfd$basis,2,1e-16)
lipmeanmarks= mean(lipmarks)

lipreglist  = landmarkreg(lipfd, as.matrix(lipmarks),
                         lipmeanmarks, WfdPar)
Dlipregfd   = register.newfd(deriv.fd(lipfd,1),
                         lipreglist$warpfd, type='direct')
D2lipregfd  = register.newfd(deriv.fd(lipfd,2),
                         lipreglist$warpfd, type='direct')
xfdlist2     = list(-Dlipregfd,-lipreglist$regfd)
lipregpda   = fRegress(D2lipregfd, xfdlist2, bwtlist)

pdalist3 = pda.fd(lipreglist$regfd,bwtlist)

bwtestlist2 = lipregpda$betaestlist

# Figure 11.6

op = par(mfrow=c(2,1))
plot(bwtestlist[[1]]$fd,lwd=2,cex.lab=1.5,cex.axis=1.5,
     ylim=c(-200, 1300))
lines(pdalist3$bwtlist[[1]]$fd,lwd=2,lty=3)

plot(bwtestlist[[2]]$fd,lwd=2,cex.lab=1.5,cex.axis=1.5)
lines(pdalist3$bwtlist[[2]]$fd,lwd=2,lty=3)
par(op)

##
## Section 11.6 Details for pda.fd, eigen.fd, pda.overlay
##              and register.newfd
##



##
## Section 11.7 Some Things to Try
##
# (exercises for the reader)

##
## Section 11.8  More to Read
##
