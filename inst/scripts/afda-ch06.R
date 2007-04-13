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

# pp.  84-85.  Figure 6.1.  the first 10 females of the Berkeley growth study
with(growth, matplot(age, hgtf[, 1:10], type="b", pch="o", 
                     ylab="Height (cm.)") )

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
# beta[1]+beta[2]*integral(exp(
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
lines(dayfine, yhatfine.12, lty=2)# lambda=1e-12:too close to a straight line

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
tst <- with(onechild, smooth.monotone(x=day,
       y=height, WfdParobj=growfdPar.oneCh1000, 
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

sm.OneCh10 <- with(smoothOneCh10,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )
smoothOneCh1 <- with(smoothOneCh1,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )
smoothOneCh.3 <- with(smoothOneCh.3,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )
smoothOneCh.001 <- with(smoothOneCh.001,
   beta[2]*eval.monfd(dayFine, Wfdobj, Lfdobj=1) )

plot(dayFine, sm.OneCh10, type="l", xlab="day",
     ylab="Growth velocity (cm/day)" )
axis(3, labels=FALSE)
axis(4, labels=FALSE)

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

# p. 88, Figure 6.6.  Estiamted growth velocity of a baby

# Data not available.






# p. 89, Figure 6.7.
#Estimated growth acceleration for 10 girls in the Berkeley Growth Study
# Similar to Figure 6.4, but acceleration not velocity
# for 10 girls not one

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

zmat.growth <- matrix(1,nage,1)

cvecf <- matrix(0, nbasis, ncasef)
betaf <- matrix(0, 2, ncasef)
RMSEf <- matrix(0, 1, ncasef)

for(icase in 1:10){
  hgtf.i     <- growth$hgtf[,icase]
  smoothGirl.i <- with(growth, smooth.monotone(x=age,
          y=hgtf.i, WfdParobj=growfdPar, zmat=zmat.growth,
	  conv=0.001, active=TRUE, dbglev=0)
   Wfd     <- smoothList$Wfdobj
   beta    <- smoothList$beta
   Flist   <- smoothList$Flist
   iternum <- smoothList$iternum
   cvecm[,icase] <- Wfd$coefs
   betam[,icase] <- beta
   hgthat <- beta[1] + beta[2]*monfn(age, Wfd)
   RMSE   <- sqrt(mean((hgt - hgthat)^2*wgt)/mean(wgt))
   RMSEm[icase] <- RMSE
   cat(c(icase, iternum),paste("  ",round(Flist$f,4),
       "  ",round(RMSE, 4)))




}

