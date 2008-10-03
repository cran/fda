###
###
### Ramsey, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 8.  Registration: Aligning Features
###         for Samples of Curves
###
library(fda)

##
## Section 8.1 Amplitude and Phase Variation
##
age  = growth$age
nage <- length(age)
ageRng <- range(age)
agefine <- seq(ageRng[1], ageRng[2], length=101)

hgtf = growth$hgtf
ncasef <- dim(hgtf)[2]

norder <- 6
nbasis <- nage + norder - 2
wbasis <- create.bspline.basis(ageRng, nbasis, norder, age)

Lfdobj    <- 3          #  penalize curvature of acceleration
lambda    <- 10^(-0.5)  #  smoothing parameter
cvecf     <- matrix(0, nbasis, ncasef)
Wfd0      <- fd(cvecf, wbasis)
growfdPar <- fdPar(Wfd0, Lfdobj, lambda)

growthMon <- smooth.monotone(age, hgtf, growfdPar)
# (wait for an interative fit to each of 54 girls)

Wfd        <- growthMon$Wfd
betaf      <- growthMon$beta
hgtfhatfd  <- growthMon$yhatfd

children <- 10

# Figure 8.1
op <- par(mfrow=c(2,1))


# Need plot.monfd?
plot(growthMon[1:children])
# Error in xy.coords(x, y, xlabel, ylabel, log) :

# *****




accelFit <- predict(growthMon, agefine, 2)
str(accelFit)

# *** lambda = 10^(-.5) was too big



matplot(agefine, accelFit[, 1:children], type='l',
        ylim=c(-4, 2))
lines(agefine, apply(accelFit[, 1:children], 1, mean),
      lty='dashed', lwd=2)














par(op)

# Figure 8.2

# transfer from fda or afda?  ?????









##
## Section 8.2 Time-Warping Functions and Registration
##

# Figure 8.3

par(mfrow=c(2,2),pty="m",ask=TRUE)
for (i in children) {
    hgtffit <- eval.fd(age,     hgtfhatfd[i])
    hgtfvec <- eval.fd(agefine, hgtfhatfd[i])
    velfvec <- eval.fd(agefine, hgtfhatfd[i], 1)
    accfvec <- eval.fd(agefine, hgtfhatfd[i], 2)
    plot(age, hgtf[,i], xlim=c(1,18), ylim=c(60,200),
         xlab="", ylab=paste("Height for female",i))
    lines(agefine, hgtfvec)
    resi <- hgtf[,i] - hgtffit
    ind  <- resi >= -.7 & resi <= .7
    plot(age[ind], resi[ind], type="b", xlim=c(1,18), ylim=c(-.7,.7),
         xlab="", ylab="Residuals")
    abline(h=0, lty=2)
    plot(agefine, velfvec, type="l", xlim=c(1,18), ylim=c(0,15),
         xlab="Years", ylab="Velocity")
    plot(agefine, accfvec, type="l", xlim=c(1,18), ylim=c(-4,2),
         xlab="Years", ylab="Acceleration")
    abline(h=0, lty=2)
}

# ???????????














##
## Section 8.3 Time-Warping Functions and Registration
##

# Use locator(1) to select the age
# of the center of the pubertal growth spurt
# for the first 10 girls
PGSctr = rep(0,10)
agefine = seq(1,18,len=101)
par(mfrow=c(1,1), ask=TRUE)
for (icase in 1:10) {
accveci = predict(accfd[icase], agefine)
plot(agefine,accveci,"l", ylim=c(-6,4),
xlab="Year", ylab="Height Accel.",
main=paste("Case",icase))
lines(c(1,18),c(0,0),lty=2)
PGSctr[icase] = locator(1)$x
}

PGSctrmean = mean(PGSctr)
wbasisLM = create.monomial.basis(c(1,18), 3)
WfdParLM = fdPar(wbasisLM)

landmarkList = landmarkreg(accfd, PGSctr, PGSctrmean,
WfdParLM, TRUE)
accregfdLM = landmarkList$regfd
warpfdLM = landmarkList$warpfd

# -> bottom pannel of Figure 8.1


##
## Section 8.4 Continuous Registration with Function register.fd
##

wbasisCR = create.bspline.basis(c(1,18), 5, 4,
c(1,PGSmeanctr,18))
Wfd0CR = fd(matrix(0,5,10),wbasisCR)
regList = register.fd(mean(accregfdLM),
accregfdLM, Wfd0CR)
accregfdCR = regList$regfd
WfdCR = regList$Wfd

# Figure 8.4




# Figure 8.5






##
## Section 8.5 A Decomposition into Amplitude and Phase Sums of Squares
##

AmpPhasList = AmpPhaseDecomp(accffd, accregfdLM, Wfd)
RSQR = AmpPhasList$RSQR




##
## 8.6 Registering the Chinese Handwriting Data
##

# Figure 8.6





##
## 8.7 Details for Functions landmarkreg and register.fd
##
help(landmarkreg)
help(register.fd)

##
## Section 8.8 Some Things to Try
##
# (exercises for the reader)

##
## Section 8.8  More to Read
##
