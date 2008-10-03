###
###
### Ramsey, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 4  How to Build Functional Data Objects
###
library(fda)
##
## Section 4.1 Adding Coefficients to Bases to Define Functions
##
#  4.1.1 Coefficient Vectors, Matrices and Arrays
daybasis65 = create.fourier.basis(c(0,365), 65)
# dummy coefmat
coefmat = matrix(0, 65, 35, dimnames=list(
     daybasis65$names, CanadianWeather$place) )
tempfd. = fd(coefmat, daybasis65)

# 4.1.2 Labels for Functional Data Objects

fdnames      = list("Age (years)", "Child", "Height (cm)")

# or

fdnames      = vector('list', 3)
fdnames[[1]] = "Age (years)"
fdnames[[2]] = "Child"
fdnames[[3]] = "Height (cm)"

station      = vector('list', 35)
station[[ 1]]= "St. Johns"
#.
#.
#.
station[[35]] = "Resolute"

# Or:
station = as.list(CanadianWeather$place)

fdnames = list("Day", "Weather Station" = station,
    "Mean temperature (deg C)")

##
## 4.2 Methods for Functional Data Objects
##
bspl2 = create.bspline.basis(norder=2)
plot(bspl2)

tstFn0 = fd(c(-1, 2), bspl2)
tstFn1 = fd(c(1, 3), bspl2)

plot(tstFn0)
plot(tstFn1)

fdsumobj = tstFn0+tstFn1
plot(fdsumobj)

fddifobj = tstFn1-tstFn0
plot(fddifobj)

fdprdobj = tstFn0 * tstFn1
plot(fdprdobj)
# NOT good:  must approximate a parabola over [0, 1]
# with a straight line (i.e., same basis set)

fdsqrobj = tstFn0^2
plot(fdsqrobj)
# Again:  NOT good:  approx. a parabola with a straight line

a   = 0.5
fd.a= tstFn0^a
#Error in backsolve(Lmat, temp) :
#  NA/NaN/Inf in foreign function call (arg 4)
# ---> square root of negative numbers

a   = (-1)
fd.a= tstFn0^a
plot(fd.a)
# Approximate a hyperbola, 1/x, with a straight line:
# nonsense, as expected

Tempbasis = create.fourier.basis(c(0, 365), 65)
Tempfd    = Data2fd(day.5, CanadianWeather$dailyAv[,,'Temperature.C'],
                    Tempbasis)
meanTempfd= mean(Tempfd)
sumTempfd = sum(Tempfd)
plot((meanTempfd-sumTempfd*(1/35))^2)
# round off error, as it should be.

# Figure 4.1.
plot(Tempfd)
lines(meanTempfd, lwd=3)
# Looks right

meanTempVec = eval.fd(day.5, meanTempfd)
meanTempVec.= predict(meanTempfd, day.5)
all.equal(meanTempVec, meanTempVec.)
lines(day.5, meanTempVec, col='red', lty='dashed', lwd=4)
# Looks right

DmeanTempVec = eval.fd(day.5, meanTempfd, 1)
DmeanTempVec.= predict(meanTempfd, day.5, 1)
all.equal(DmeanTempVec, DmeanTempVec.)
plot(day.5, DmeanTempVec, type='l')

accelLfd     = int2Lfd(2)

harmaccelLfd = vec2Lfd(c(0,c(2*pi/365)^2, 0), c(0, 365))
LmeanTempVec = eval.fd(day.5, meanTempfd, harmaccelLfd)
LmeanTempVec.= predict(meanTempfd, day.5, harmaccelLfd)
all.equal(LmeanTempVec, LmeanTempVec.)

# daytime = day.5
# JJindex = dayOfYearShifted

tempmat  = daily$tempav[dayOfYearShifted, ]
tempbasis= create.fourier.basis(c(0,365),65)

temp.fdSmooth = smooth.basis(day.5, tempmat, tempbasis)
tempfd        = temp.fdSmooth$fd
tempfd$fdnames= list("Day (July 2 to June 30)",
    "Weather Station",
    "Mean temperature (deg. C)")
plot(tempfd, col=1, lty=1, xlab='Day (July 1 to June 30)',
     ylab='Mean temperature (deg. C)')

# Section 4.2.1 Illustration: Sinusoidal Coefficients
# Figure 4.2

basis13  = create.bspline.basis(c(0,10), 13)
tvec     = seq(0,1,len=13)
sinecoef = sin(2*pi*tvec)
sinefd   = fd(sinecoef, basis13, list("t","","f(t)"))
op       = par(cex=1.2)
plot(sinefd, lwd=2)
points(tvec*10, sinecoef, lwd=2)
par(op)

##
## Section 4.3 Smoothing using Regression Analysis
##

# Section 4.3.1 Plotting the January Thaw

# Figure 4.3

# This assumes the data are in "MtlDaily.txt"
# in the working directory getwd()
# MtlDaily = matrix(scan("MtlDaily.txt",0),34,365)
# thawdata = t(MtlDaily[,16:47])
thawdata = t(MontrealTemp[, 16:47])
daytime  = ((16:47)+0.5)
plot(daytime, apply(thawdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Temperature (deg C)", cex=1.2)

thawbasis   = create.bspline.basis(c(16,48),7)
thawbasismat= eval.basis(thawbasis, daytime)

# Figure 4.4

thawcoef = solve(crossprod(thawbasismat),
    crossprod(thawbasismat,thawdata))
thawfd   = fd(thawcoef, thawbasis,
    list("Day", "Year", "Temperature (deg C)"))
plot(thawfd, lty=1, lwd=2, col=1)

# Figure 4.5

plotfit.fd(thawdata[,1], daytime, thawfd[1],
           lty=1, lwd=2, main='')

##
## Section 4.4 The Linear Differential Operator or Lfd Class
##
omega          = (2*pi/365)
thawconst.basis= create.constant.basis(thawbasis$rangeval)

betalist       = vector("list", 3)
betalist[[1]]  = fd(0, thawconst.basis)
betalist[[2]]  = fd(omega^2, thawconst.basis)
betalist[[3]]  = fd(0, thawconst.basis)
harmaccelLfd.  = Lfd(3, betalist)

accelLfd = int2Lfd(2)

harmaccelLfd.thaw = vec2Lfd(c(0,omega^2,0), thawbasis$rangeval)
all.equal(harmaccelLfd.[-1], harmaccelLfd.thaw[-1])

class(accelLfd)
class(harmaccelLfd)

Ltempmat = eval.fd(day.5, tempfd, harmaccelLfd)
Ltempmat.= predict(tempfd, day.5, harmaccelLfd)
all.equal(Ltempmat, Ltempmat.)

D2tempfd = deriv.fd(tempfd, 2)
Ltempfd  = deriv.fd(tempfd, harmaccelLfd)

##
## Section 4.5 Bivariate Functional Data Objects:
##             Functions of Two Arguments
##
Bspl2 = create.bspline.basis(nbasis=2, norder=1)
Bspl3 = create.bspline.basis(nbasis=3, norder=2)

corrmat = array(1:6/6, dim=2:3)
bBspl2.3= bifd(corrmat, Bspl2, Bspl3)

##
## Section4.6 The structure of the fd and Lfd Classes
##
help(fd)
help(Lfd)

##
## Section 4.7 Some Things to Try
##
# (exercises for the reader)
