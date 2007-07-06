###
###
### Ramsey & Silverman (2002) Applied Functional Data Analysis (Springer)
###
### ch. 6.  Human growth 
###
library(fda)

##
## sec. 6.1.  Introduction
##

##
## sec. 6.2.  Height measurements at three scales
##
str(growth)

# pp.  84-85.  Figure 6.1.  the first 10 females
# of the Berkeley growth study
op <- par(cex=)
with(growth, matplot(age, hgtf[, 1:10], type="b", pch="o", 
                     ylab="Height (cm.)") )
par(op)
# Figure 6.2.  Heights of one boy during one school year

# A monotone smooth requires some effort.  
#  (1) set up the basis

nbasis.onechild   <- 33
# Establish a B-spline basis
# with nbasis.onechild basis function
hgtbasis <- with(onechild, 
  create.bspline.basis(range(day), nbasis.onechild))

tst <- create.bspline.basis(onechild$day)
#  set up the functional data object for W <- log Dh

# Start by creating a functional 0 from hgtbasis
cvec0 <- rep(0,nbasis.onechild)
Wfd0   <- fd(cvec0, hgtbasis)

#  set parameters for the monotone smooth
# with smoothing lambda = 1e-1 or 1e-12 

WfdPar.1  <- fdPar(Wfd0, 2, lambda=.1)
WfdPar.12  <- fdPar(Wfd0, 2, lambda=1e-12)

#  --------------   carry out the monotone smooth  ---------------

# The monotone smooth is
# beta[1]+beta[2]*integral(exp(Wfdobj)),
#     where Wfdobj is a functional data object
smoothList.1   <- with(onechild, 
  smooth.monotone(x=day, y=height, WfdParobj=WfdPar.1) )
smoothList.12   <- with(onechild, 
  smooth.monotone(x=day, y=height, WfdParobj=WfdPar.12) )

#str(smoothList)
#attach(smoothList)

# Create a fine grid at which to evaluate the smooth
dayfine  <- with(onechild, seq(day[1],day[length(day)],len=151))

# eval.monfd = integral(exp(Wfdobj))
# This is monotonically increasing, since exp(Wfdobj)>0.  
#yhat     <- with(smoothList,
#             beta[1] + beta[2]*eval.monfd(onechild$day, Wfdobj))
yhatfine.1 <- with(smoothList.1, 
       beta[1] + beta[2]*eval.monfd(dayfine, Wfdobj))
yhatfine.12 <- with(smoothList.12, 
       beta[1] + beta[2]*eval.monfd(dayfine, Wfdobj))

plot(onechild, ylab="Height (cm.)") # raw data
lines(dayfine, yhatfine.1, lwd=2) # lambda=0.1:  reasonable
lines(dayfine, yhatfine.12, lty=2)
# lambda=1e-12:too close to a straight line

# p.  86, Figure 6.3.  Growth of the length of the tibia of a newborn
# data not available





##
## sec. 6.3.  Velocity and acceleration
##

# p. 87, Figure 6.4.  Estimated growth rate of the first girl 

nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
# 35 
rng.growth <- range(growth$age)
# 1 18 
wbasis.growth <- create.bspline.basis(rng.growth,
                   nbasis.growth, norder.growth,
                   growth$age)

#  starting values for coefficient

cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth  <- fd(cvec0.growth, wbasis.growth)

Lfdobj.growth    <- 3          #  penalize curvature of acceleration
lambda.growth    <- 10^(-0.5)  #  smoothing parameter
growfdPar <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.growth)

#  ---------------------  Now smooth the data  --------------------

smoothGirl1 <- with(growth, smooth.monotone(x=age,
        y=hgtf[, 1], WfdParobj=growfdPar, conv=0.001, active=TRUE, 
        dbglev=0) )
#*** The default active = c(FALSE, rep(TRUE, ncvec-1))
#    This tells smooth.monotone to estimate
#    force the first coeffeicient of W to 0,
#    and estimate only the others.
#*** That starts velocity unrealistically low.
#*** Fix this with active=TRUE 
#    to estimate all elements of W.  

# Create a fine grid at which to evaluate the smooth
agefine  <- with(growth, seq(age[1], age[nage], len=151))

# Lfdobj = 1 for first derivative, growth rate or velocity 
smoothG1.1 <- with(smoothGirl1,
       beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=1))

plot(agefine, smoothG1.1, type="l",
     xlab="Year", ylab="Growth velocity (cm/year)")
axis(3, labels=FALSE)
axis(4, labels=FALSE)


# p. 88, Figure 6.5, Estimated growth velocity of a 10-year old boy
str(onechild)

nDays <- dim(onechild)[1]
norder.oneCh <- 6
nbasis.oneCh <- nDays+norder.oneCh-2
rng.days <- range(onechild$day)

# B-spline basis 
wbasis.oneCh <- create.bspline.basis(rng.days,
       nbasis.oneCh, norder.oneCh, onechild$day)

# starting values for coefficients
cvec0.oneCh <- matrix(0, nbasis.oneCh, 1)
Wfd0.oneCh <- fd(cvec0.oneCh, wbasis.oneCh)

Lfdobj.oneCh    <- 3
#  penalize curvature of acceleration
growfdPar.oneCh1000 <- fdPar(Wfd0.oneCh, Lfdobj.oneCh,
                 lambda=1000 )
growfdPar.oneCh100 <- fdPar(Wfd0.oneCh, Lfdobj.oneCh,
                 lambda=100 )
growfdPar.oneCh10 <- fdPar(Wfd0.oneCh, Lfdobj.oneCh,
                 lambda=10 )
growfdPar.oneCh1 <- fdPar(Wfd0.oneCh, Lfdobj.oneCh,
                 lambda=1 )
growfdPar.oneCh.3 <- fdPar(Wfd0.oneCh, Lfdobj.oneCh,
                 lambda=sqrt(0.1) )
growfdPar.oneCh.001 <- fdPar(Wfd0.oneCh, Lfdobj.oneCh,
                 lambda=.001)

# now smooth the data

zmat.oneCh <- matrix(1, nDays, 1)
smoothOneCh1000 <- with(onechild, smooth.monotone(x=day,
       y=height, WfdParobj=growfdPar.oneCh1000, zmat=zmat.oneCh,
       conv=0.001, active=TRUE, dbglev=0) )
smoothOneCh100 <- with(onechild, smooth.monotone(x=day,
       y=height, WfdParobj=growfdPar.oneCh100, 
       conv=0.001, active=TRUE, dbglev=0) )
smoothOneCh10 <- with(onechild, smooth.monotone(x=day,
       y=height, WfdParobj=growfdPar.oneCh10, zmat=zmat.oneCh,
       conv=0.001, active=TRUE, dbglev=0) )
smoothOneCh1 <- with(onechild, smooth.monotone(x=day,
       y=height, WfdParobj=growfdPar.oneCh1, zmat=zmat.oneCh,
       conv=0.001, active=TRUE, dbglev=0) )
smoothOneCh.3 <- with(onechild, smooth.monotone(x=day,
       y=height, WfdParobj=growfdPar.oneCh.3, zmat=zmat.oneCh,
       conv=0.001, active=TRUE, dbglev=0) )
smoothOneCh.001 <- with(onechild, smooth.monotone(x=day,
       y=height, WfdParobj=growfdPar.oneCh.001, zmat=zmat.oneCh,
       conv=0.001, active=TRUE, dbglev=0) )

dayFine <- with(onechild, seq(day[1], day[nDays], len=151))

sm.OneCh1000 <- with(smoothOneCh1000,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )
sm.OneCh100 <- with(smoothOneCh100,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )
sm.OneCh10 <- with(smoothOneCh10,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )
smoothOneCh1 <- with(smoothOneCh1,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )
smoothOneCh.3 <- with(smoothOneCh.3,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )
smoothOneCh.001 <- with(smoothOneCh.001,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )

plot(dayFine, sm.OneCh1000, type="l", xlab="day",
     ylab="Growth velocity (cm/day)" )
axis(3, labels=FALSE)
axis(4, labels=FALSE)
# Too smooth 

plot(dayFine, sm.OneCh100, type="l", xlab="day",
     ylab="Growth velocity (cm/day)" )
axis(3, labels=FALSE)
axis(4, labels=FALSE)
# Close to Figure 6.5

plot(dayFine, sm.OneCh10, type="l", xlab="day",
     ylab="Growth velocity (cm/day)" )
axis(3, labels=FALSE)
axis(4, labels=FALSE)
# Shows more features than in Figure 6.5.  

plot(dayFine, smoothOneCh1, type="l", xlab="day",
     ylab="Growth velocity (cm/day)" )
axis(3, labels=FALSE)
axis(4, labels=FALSE)

plot(dayFine, smoothOneCh.3, type="l", xlab="day",
     ylab="Growth velocity (cm/day)" )
axis(3, labels=FALSE)
axis(4, labels=FALSE)

plot(dayFine, smoothOneCh.001, type="l", xlab="day",
     ylab="Growth velocity (cm/day)" )
axis(3, labels=FALSE)
axis(4, labels=FALSE)

# p. 88, Figure 6.6.  Estimated growth velocity of a baby

# Data not available.






# p. 89, Figure 6.7.
#Estimated growth acceleration for 10 girls in the Berkeley Growth Study
# Similar to Figure 6.4, but acceleration not velocity
# and for 10 girls not one

nage <- length(growth$age)
norder.growth <- 6
nbasis.growth <- nage + norder.growth - 2
# 35 
rng.growth <- range(growth$age)
# 1 18 
wbasis.growth <- create.bspline.basis(rng.growth,
                   nbasis.growth, norder.growth,
                   growth$age)
str(wbasis.growth)
#  starting values for coefficient

cvec0.growth <- matrix(0,nbasis.growth,1)
Wfd0.growth  <- fd(cvec0.growth, wbasis.growth)

Lfdobj.growth    <- 3          #  penalize curvature of acceleration
lambda.growth    <- 10^(-0.5)  #  smoothing parameter
growfdPar <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.growth)

#  ---------------------  Now smooth the data  --------------------

ncasef <- 10
smoothGirls <- vector("list", ncasef)
for(icase in 1:ncasef){
  smoothGirls[[icase]] <- with(growth, smooth.monotone(x=age,
      y=hgtf[, icase], WfdParobj=growfdPar, conv=0.001,
      active=TRUE, dbglev=0) )
  cat(icase, "")
}

# Create a fine grid at which to evaluate the smooth
nptsFine <- 151

agefine  <- with(growth, seq(age[1], age[nage], len=nptsFine))
smoothGirlsHt <- array(NA, dim=c(nptsFine, ncasef))
smoothGirlsAcc <- smoothGirlsVel <- smoothGirlsHt
for(icase in 1:ncasef)
  smoothGirlsHt[, icase] <- with(smoothGirls[[icase]],
       beta[1]+beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=0))

with(growth, plot(age, hgtf[, 1]))
lines(agefine, smoothGirlsHt[, 1])

for(icase in 1:ncasef)
  smoothGirlsVel[, icase] <- with(smoothGirls[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=1) )

tstVel <- diff(smoothGirlsHt[, 1])/diff(agefine)
quantile(smoothGirlsHt[, 1])
quantile(smoothGirlsVel[, 1])
quantile(tstVel)

plot(smoothGirlsVel[, 1])
lines(tstVel)
# GOOD:  eval.monfd computes velocity as expected 

for(icase in 1:ncasef)
  smoothGirlsAcc[, icase] <- with(smoothGirls[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=2) )

tstAcc <- diff(smoothGirlsVel[, 1])/diff(agefine)
quantile(smoothGirlsAcc[, 1])
quantile(tstAcc)

plot(smoothGirlsAcc[, 1])
lines(tstAcc)
# GOOD:  eval.monfd computes acceleration as expected 

op <- par(cex=1.5)
matplot(agefine, smoothGirlsAcc, type="l", ylim=c(-4, 2),
        xlab="Year", ylab="Growth acceleration (cm/year^2)")
lines(agefine, apply(smoothGirlsAcc, 1, mean), lwd=3)
par(op)
# OK, but the image doesn't quite match Figure 6.7.
# It's either smoother or it represents different girls (?) 










#############################################
lambda.gr.1 <- .1
growfdPar.1 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr.1)

ncasef <- 10
smoothGirls.1 <- vector("list", ncasef)
for(icase in 1:ncasef){
  smoothGirls.1[[icase]] <- with(growth, smooth.monotone(x=age,
      y=hgtf[, icase], WfdParobj=growfdPar.1, conv=0.001,
      active=TRUE, dbglev=0) )
  cat(icase, "")
}
    
agefine  <- with(growth, seq(age[1], age[nage], len=nptsFine))
smoothGirlsAcc.1 <- array(NA, dim=c(nptsFine, ncasef))
for(icase in 1:ncasef)
  smoothGirlsAcc.1[, icase] <- with(smoothGirls.1[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=2) )

matplot(agefine, smoothGirlsAcc.1, type="l", ylim=c(-4, 2),
        xlab="Year", ylab="Growth acceleration (cm/year^2)")
lines(agefine, apply(smoothGirlsAcc.1, 1, mean), lwd=3)
# Still possibly oversmoothed,
# but perhaps not as bad as with lambda = 10^(-0.5)???

#############################################
lambda.gr2.3 <- .03
growfdPar2.3 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr2.3)

ncasef <- 10
smoothGirls2.3 <- vector("list", ncasef)
for(icase in 1:ncasef){
  smoothGirls2.3[[icase]] <- with(growth, smooth.monotone(x=age,
      y=hgtf[, icase], WfdParobj=growfdPar2.3, conv=0.001,
      active=TRUE, dbglev=0) )
  cat(icase, "")
}
    
agefine  <- with(growth, seq(age[1], age[nage], len=nptsFine))
smoothGirlsAcc2.3 <- array(NA, dim=c(nptsFine, ncasef))
for(icase in 1:ncasef)
  smoothGirlsAcc2.3[, icase] <- with(smoothGirls2.3[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=2) )

matplot(agefine, smoothGirlsAcc2.3, type="l", ylim=c(-4, 2),
        xlab="Year", ylab="Growth acceleration (cm/year^2)")
lines(agefine, apply(smoothGirlsAcc2.3, 1, mean), lwd=3)
abline(h=0, lty="dashed")
# Good match


op <- par(cex=1.5)
matplot(agefine, smoothGirlsAcc2.3, type="l", ylim=c(-4, 2),
        xlab="age", ylab="Growth acceleration (cm/year^2)",
        bty="n")
lines(agefine, apply(smoothGirlsAcc2.3, 1, mean), lwd=3)
abline(h=0, lty="dashed")
par(op)

#############################################
lambda.gr2 <- .01
growfdPar2 <- fdPar(Wfd0.growth, Lfdobj.growth, lambda.gr2)

ncasef <- 10
smoothGirls2 <- vector("list", ncasef)
for(icase in 1:ncasef){
  smoothGirls2[[icase]] <- with(growth, smooth.monotone(x=age,
      y=hgtf[, icase], WfdParobj=growfdPar2, conv=0.001,
      active=TRUE, dbglev=0) )
  cat(icase, "")
}
    
agefine  <- with(growth, seq(age[1], age[nage], len=nptsFine))
smoothGirlsAcc2 <- array(NA, dim=c(nptsFine, ncasef))
for(icase in 1:ncasef)
  smoothGirlsAcc2[, icase] <- with(smoothGirls2[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=2) )

matplot(agefine, smoothGirlsAcc2, type="l", ylim=c(-4, 2),
        xlab="Year", ylab="Growth acceleration (cm/year^2)")
lines(agefine, apply(smoothGirlsAcc2, 1, mean), lwd=3)
# Clearly undersmoothed, but not as much as
# with any smaller lambda 

##
## sec. 6.4.  An equation for growth 
##

# p. 91, Figure 6.8.  Relative acceleration and its integral

#Estimated growth acceleration for 10 girls in the Berkeley Growth Study
# Similar to Figure 6.7, but acceleration / velocity

smoothGirlsVel2.3 <- array(NA, dim=c(nptsFine, ncasef))
for(icase in 1:ncasef)
  smoothGirlsVel2.3[, icase] <- with(smoothGirls2.3[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=1) )

smoothGirlsRelAcc2.3 <- (smoothGirlsAcc2.3 /
                         smoothGirlsVel2.3) 

matplot(agefine, smoothGirlsRelAcc2.3, type="l", ylim=c(-2, 0.5),
        xlab="Year", ylab="w(t)")
abline(h=0, lty="dashed")

diff(range(quantile(diff(agefine))))
d.agefine <- mean(diff(agefine))
smoothGirlsIntRelAcc2.3 <- (apply(smoothGirlsRelAcc2.3, 2, cumsum)
                            * d.agefine)
str(smoothGirlsIntRelAcc2.3)

matplot(agefine, smoothGirlsIntRelAcc2.3, type="l", ylim=c(-4, 0.5),
        xlab="Year", ylab="W(t)")
abline(h=0, lty="dashed")

##
## sec. 6.5.  Timing or phase variation in growth 
##

# p. 92, Figure 6.9.  Acceleration curves
#  differing in amplitude only and phase only

# read numbers off the figure:  
Fig6.9 <- cbind(Age=c(5, 6, 9, 11, 13, 15, 17, 20),
                accel=c(-.5, -.4, -.5, 0, 1.4, -2.8, -.6, 0) )
(n.Fig6.9 <- dim(Fig6.9)[1])

Fig6.9basis <- create.bspline.basis(range(Fig6.9[, 1]), 
                           breaks=Fig6.9[, 1])
str(Fig6.9basis)

#**** I can't make this work 

Fig6.9s <- smooth.basisPar(argvals=Fig6.9[, 1], y=Fig6.9[, 2],
               fdobj=Fig6.9basis)
#Error in df < n : comparison (3) is possible only for atomic and list types

Fig6.9basis. <- create.bspline.basis(range(Fig6.9[, 1]), n.Fig6.9,
                           breaks=Fig6.9[, 1])
Fig6.9basis.10 <- create.bspline.basis(range(Fig6.9[, 1]), 10,
                           breaks=Fig6.9[, 1])
Fig6.9basis.12 <- create.bspline.basis(range(Fig6.9[, 1]), 12,
                           breaks=Fig6.9[, 1])
Fig6.9basis.6 <- create.bspline.basis(range(Fig6.9[, 1]), 6,
                           breaks=Fig6.9[, 1])
str(Fig6.9basis)
str(Fig6.9basis.)
str(Fig6.9basis.6)
str(Fig6.9basis.10)
str(Fig6.9basis.12)

Fig6.9sP <- smooth.basisPar(Fig6.9[, 1], Fig6.9[, 2], Fig6.9basis)
#Error in df < n : comparison (3) is possible only for atomic and list types
Fig6.9sP. <- smooth.basisPar(Fig6.9[, 1], Fig6.9[, 2], Fig6.9basis.)
# OK
Fig6.9sP.6 <- smooth.basisPar(Fig6.9[, 1], Fig6.9[, 2], Fig6.9basis.6)
#Error in bsplineS(evalarg, breaks, norder, nderiv) : 
#	NORDER less than 1.
Fig6.9sP.10 <- smooth.basisPar(Fig6.9[, 1], Fig6.9[, 2], Fig6.9basis.10)
#Error in df < n : comparison (3) is possible only for atomic and list types
Fig6.9sP.12 <- smooth.basisPar(Fig6.9[, 1], Fig6.9[, 2], Fig6.9basis.12)
#Error in df < n : comparison (3) is possible only for atomic and list types

Fig6.9sP1 <- smooth.basisPar(Fig6.9[, 1], Fig6.9[, 2], Fig6.9basis.,
                             lambda=1)

str(Fig6.9sP)
plot(Fig6.9sP$fd)
points(Fig6.9)
Age <- seq(5, 20, len=61)
lines(Age, eval.fd(Age, Fig6.9sP$fd)) 
# straight line segments, not smooth curves ????? 
plot(Fig6.9sP$fd)
Age <- seq(5, 20, len=61)
lines(Age, eval.fd(Age, Fig6.9sP$fd)) 
# straight line segments, not smooth curves ????? 

# Start by creating a functional 0 from hgtbasis
Fig6.9vec0 <- rep(0,n.Fig6.9)
Fig6.9fd0   <- fd(Fig6.9vec0, Fig6.9basis)
#  set parameters for the monotone smooth
# with smoothing lambda = 1e-1 or 1e-12 

Fig6.9Par.1  <- fdPar(Fig6.9fd0, 2, lambda=.1)
Fig6.9Par.12  <- fdPar(Fig6.9fd0, 2, lambda=1e-12)

Fig6.9s <- smooth.basis(Fig6.9[, 1], Fig6.9[, 2], Fig6.9basis)

Fig6.9fd. <- data2fd(Fig6.9[, 2], Fig6.9[, 1], Fig6.9basis)

str(Fig6.9fd.)
plot(Fig6.9fd.)
#  crazy between 17 and 20

Fig6.9s.1 <- smooth.basisPar(Fig6.9[, 1], Fig6.9[, 2], Fig6.9basis,
                           lambda=.1)
plot(Fig6.9s.1$fd)
plot(Fig6.9s$fd)

nptsFine6.9 <- 151
ageFine6.9 <- seq(5, 20, length=nptsFine6.9)
plot(Fig6.9fd., Fig6.9[,1])
plot(Fig6.9fd., ageFine6.9)
#  crazy between 17 and 20

Fig6.9fine <- eval.fd(ageFine6.9, Fig6.9fd.)
Fig6.9s. <- eval.fd(ageFine6.9, Fig6.9s$fd)
str(Fig6.9s.)
range(Fig6.9s.)

str(Fig6.9fine)
range(Fig6.9fine)

plot(Fig6.9)
lines(ageFine6.9, Fig6.9s.)
# over smoothing or ...??? 

plot(Fig6.9, ylim=range(Fig6.9fine))
lines(ageFine6.9, Fig6.9s.)
# over smoothing or ...??? 

str(ageFine6.9)
str(Fig6.9fine)
plot(Fig6.9fine[, 1], 3*Fig6.9fine[, 2], type="l")
lines(Fig6.9fine)
points(Fig6.9)
plot(Fig6.9, type="l")

#****** give up on this for now ... :(

# p. 93, Figure 6.10.  Time warping functions
#  for 10 Berkeley girls

str(smoothGirlsAcc)

smGirlsAccAvg <- apply(smoothGirlsAcc, 1, mean)
str(smGirlsAccAvg)

warp10GirlsMon <- vector('list', 10)
for(i in 1:10)
  warp10GirlsMon[[i]] <- register.fd(smoothGirlsMean, smoothGirls[[i]]$Wfdobj,
                                     growfdPar2.3)   
agefine  <- with(growth, seq(age[1], age[nage], len=nptsFine))
GirlsAccReg2.3 <- array(NA, dim=c(nptsFine, ncasef))
for(i in 1:10){
  bi2 <- smoothGirls[[i]]$beta[2]
  GirlsAccReg2.3[, i] <- bi2*eval.monfd(agefine, 
                                        warp10GirlsMon[[i]]$regfd, Lfdobj=2)
}
op <- par(cex=1.5)
matplot(agefine, GirlsAccReg2.3, type="l", ylim=c(-4, 2),
        xlab="warped age", ylab="Growth acceleration (cm/yr^2)",
        bty="n")
lines(agefine, apply(GirlsAccReg2.3, 1, mean), lwd=3)
abline(h=0, lty="dashed")
par(op)

# Try crit = 1, because the default crit = 2 didn't seem great 
warp10GirlsMonl <- vector('list', 10)
for(i in 1:10)
  warp10GirlsMonl[[i]] <- register.fd(smoothGirlsMean, smoothGirls[[i]]$Wfdobj,
                                     growfdPar2.3, crit=1)   
agefine  <- with(growth, seq(age[1], age[nage], len=nptsFine))
GirlsAccReg2.3l <- array(NA, dim=c(nptsFine, ncasef))
for(i in 1:10){
  bi2 <- smoothGirls[[i]]$beta[2]
  GirlsAccReg2.3l[, i] <- bi2*eval.monfd(agefine, 
                                        warp10GirlsMonl[[i]]$regfd, Lfdobj=2)
}
op <- par(cex=1.5)
matplot(agefine, GirlsAccReg2.3, type="l", ylim=c(-4, 2),
        xlab="warped age", ylab="Growth acceleration (cm/yr^2)",
        bty="n")
lines(agefine, apply(GirlsAccReg2.3, 1, mean), lwd=3)
abline(h=0, lty="dashed")
par(op)
                                

  
#warpmat = monfn(agefine, Wfd);
#warpmat = 1 + 17.*warpmat./(ones(101,1)*warpmat(101,:));

#for i = 1:length(index)
#   subplot(1,2,1)
#  plot(agefine, yvec(:,i), '-', agefine, y0vec, '--', agefine, yregmat(:,i), '-');
#  axis('square')
#  title(['Case ',num2str(i)])
#  subplot(1,2,2)
#  plot(agefine, warpmat(:,i), '-', agefine, agefine, '--')
#   axis('square')
#   pause
#end

  

  
for(icase in 1:ncasef)
  smoothGirlsAcc2.3[, icase] <- with(smoothGirls2.3[[icase]],
     beta[2]*eval.monfd(agefine, Wfdobj, Lfdobj=2) )






                              


op <- par(cex=1.5)
plot(range(agefine), c(-4, 2), type="n", 
matplot(agefine, smoothGirlsAcc2.3, type="l",
        xlab="age (registered to mean)",
        ylab="Growth acceleration (cm/year^2)", bty="n")

     


     
lines(agefine, apply(smoothGirlsAcc2.3, 1, mean), lwd=3)
abline(h=0, lty="dashed")
par(op)

















nptsCoarse <- 18

ageCoarse  <- with(growth, seq(age[1], age[nage], len=nptsCoarse))
smoothGirlsAcc2.3c <- array(NA, dim=c(nptsCoarse, ncasf))
for(icase in 1:ncasef)
  smoothGirlsAcc2.3c[, icase] <- with(smoothGirls2.3[[icase]],
     beta[2]*eval.monfd(ageCoarse, Wfdobj, Lfdobj=2) )

matplot(agefine, smoothGirlsAcc2.3, type="l", ylim=c(-4, 2),
        xlab="Year", ylab="Growth acceleration (cm/year^2)")
for(icase in 1:ncasef)
  text(ageCoarse, smoothGirlsAcc2.3c[, icase], icase)

goodGirls <- c(1:2, 5:6, 8:9)
matplot(agefine, smoothGirlsAcc2.3[, goodGirls], type="l",
        ylim=c(-4, 2),
        xlab="Year", ylab="Growth acceleration (cm/year^2)")
for(icase in goodGirls)
  text(ageCoarse, smoothGirlsAcc2.3c[, icase], icase)

refGirl <- 9


#smoothGirls2.3


smoothGirls2.3
tst9.1 <- with(growth, register.fd(smoothGirls2.3[[9]]$Wfdobj,
                                   smoothGirls2.3[[1]]$Wfdobj) )








matplot(ageCoarse, smoothGirlsAcc2.3, type="l", ylim=c(-4, 2),
        xlab="Year", ylab="Growth acceleration (cm/year^2)")
lines(agefine, apply(smoothGirlsAcc2.3, 1, mean), lwd=3)
abline(h=0, lty="dashed")


hgtRng <- range(growth$hgtf[, 1:10])
with(growth, matplot(age, hgtf[, 1:10], type="b", pch="o", 
                     ylab="Height (cm.)") )
