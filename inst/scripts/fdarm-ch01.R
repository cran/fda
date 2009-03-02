###
###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 1.  Introduction
###
library(fda)
##
## Section 1.1.  What are functional data?
##
# Prep for Figure 1.1
age     = growth$age
age.rng = range(age)
agefine = seq(age.rng[1], age.rng[2], length=501)

# Monotone smooth (see section 5.4.2)
# B-spline of order 6 = quintic polynomials
# so the acceleration will be cubic
gr.basis = create.bspline.basis(norder=6, breaks=growth$age)

# Consider only the first 10 girls
children= 1:10
ncasef  = length(children)

# starting values for coeficients
cvecf           = matrix(0, gr.basis$nbasis, ncasef)
dimnames(cvecf) = list(gr.basis$names,
              dimnames(growth$hgtf)[[2]][children])
# Create an initial functional data object
gr.fd0  = fd(cvecf, gr.basis)
# Create an initial functional parameter object
gr.fdPar= fdPar(gr.fd0, Lfdobj=3, lambda=1/sqrt(10))
# Lfdobj = 3 to penalize the rate of change of acceleration

#  -----------------  Now smooth the female data  --------------------

# Figure 1.15 was created with lambda = 10^(-1.5);
# we use that also for Figure 1.1

# Figure 1.1

gr.fdPar1.5 = fdPar(gr.fd0, Lfdobj=3, lambda=10^(-1.5))
hgtfmonfd   = with(growth, smooth.monotone(age, hgtf[,children],
                                           gr.fdPar1.5) )

with(growth, matplot(age, hgtf[, children], pch='o',
                     xlab='Age (years)', ylab='Height (cm)',
                     ylim=c(60, 183)) )
hgtf.vec1.5 = predict(hgtfmonfd$yhatfd, agefine)
matlines(agefine, hgtf.vec1.5, lty=1)

# Figure 1.2

accfvec1.5 = predict(hgtfmonfd$yhatfd, agefine, Lfdobj=2)
matplot(agefine, accfvec1.5, type='l', lty=1, ylim=c(-4, 2),
        xlab='Age (years)', ylab=expression(Acceleration (cm/yr^2)),
        xlim=c(1, 18), las=1)
abline(h=0, lty='dotted')
lines(agefine, rowMeans(accfvec1.5), lty='dashed', lwd=2)

# Figure 1.3

plot(nondurables, ylim=c(0, 120), 
     xlab = '', ylab='Nondurable Goods Index', las=1)

# Figure 1.4

plot(log10(nondurables), xlab = '', 
     ylab=expression(Log[10]~~Nondurable~~Goods~~Index), las=1 )
abline(lm(log10(nondurables) ~ index(nondurables)), lty='dashed')

# Figure 1.5

op = par(mfrow=c(2,1), mar=c(4, 4, 1, 2)+.1)
plot(Tray47~Time, refinery, pch='*', xlab='', ylab='Tray 47 level')
plot(Reflux~Time, refinery, pch='*', xlab='Time (min)',
     ylab='Reflux flow')
par(op)

##
## Section 1.2.  Multivariate functional data
##

# Figure 1.6

Time = c(0, as.numeric(dimnames(gait)[[1]]), 1)
# Interpolate a common value for 0 & 1
gait01        = 0.5*(gait[20,,]+gait[1,,])
Gait          = array(NA, dim=c(22, 39, 2))
dimnames(Gait)= c(list(Time), dimnames(gait)[2:3])
Gait[1,,]     = gait01
Gait[2:21,,]  = gait
Gait[22,,]    = gait01

before = 17:22
Before = c(seq(-0.225, -0.025, by=0.05), 0)

after = 1:6
After = c(1, seq(1.025, 1.225, by=0.05))

op = par(mfrow=c(2,1), mar=c(4, 4, 1, 2)+.1)
matplot(Time, Gait[,,'Hip Angle'], ylab='Hip Angle (degrees)',
        type='l', xlim=c(-0.25, 1.25), lty='solid', xlab='')
matlines(Before, Gait[before,,'Hip Angle'], lty='dotted')
matlines(After, Gait[after,,'Hip Angle'], lty='dotted')
abline(v=0:1, lty='dotted')

matplot(Time, Gait[,,'Knee Angle'], ylab='Knee Angle (degrees)',
        type='l', xlim=c(-0.25, 1.25), lty='solid',
        xlab='Time (portion of gait cycle)')
matlines(Before, Gait[before,,'Knee Angle'], lty='dotted')
matlines(After, gait[after,,'Knee Angle'], lty='dotted')
abline(v=0:1, lty='dotted')
par(op)

# Figure 1.7

gaitMean = apply(gait, c(1, 3), mean)

xlim = range(c(gait[,,1], gaitMean[, 1]))
ylim = range(c(gait[,,2], gaitMean[, 2]))

plot(gait[, 1, ], type='b', xlim=xlim, ylim=ylim,
     xlab='Hip angle (degrees)', ylab='Knee angle (degrees)',
     pch='.', cex=3, lwd=1.5)
points(gaitMean, type='b', lty='dotted', pch='.',
       cex=3, lwd=1.5)
i4 = seq(4, 20, 4)
text(gait[i4, 1, ], labels=LETTERS[c(2:5, 1)])
text(gaitMean[i4, ], labels=LETTERS[c(2:5, 1)])

# Figure 1.8

matplot(100*handwrit[, , 1], 100*handwrit[, , 2], type="l",
        lty='solid', las=1, xlab='', ylab='')

# Figure 1.9

mark = seq(1, 601, 12)

i    = 1
StatSci1 = StatSciChinese[, i, ]
# Where does the pen leave the paper?
thresh = quantile(StatSci1[, 3], .8)

sel1 = (StatSci1[, 3] < thresh)
StatSci1[!sel1, 1:2] = NA
plot(StatSci1[, 1:2], type='l', lwd=2)
points(StatSci1[mark, 1], StatSci1[mark, 2], lwd=2)

##
## Section 1.3.  Functional models for nonfunctional data
##

# Figure 1.10:

# Based on a functional logistic regression.
# Unfortunately, the current 'fda' package does NOT include
# code to estimate a functional logistic regression model.

##
## Section 1.4.  Some functional data analyses
##

fig1.11Stns = c('Montreal', 'Edmonton', 'Pr. Rupert', 'Resolute')
fig1.11Temp = CanadianWeather$dailyAv[, fig1.11Stns, 'Temperature.C']

# Comment from the end of section 1.5.1:
# "the temperature data in Figure 1.11 were fit using
#  smoothing splines".

Temp.fourier   = create.fourier.basis(c(0, 365), 13)
fig1.11Temp.fd = Data2fd(day.5, fig1.11Temp, Temp.fourier)

# Figure 1.11

plot(fig1.11Temp.fd, day.5, axes=FALSE, col=1, lwd=2,
     xlab='', ylab='Mean Temperature (deg C)')
axis(2, las=1)
axisIntervals(labels=monthLetters)

monthIndex = rep(1:12, daysPerMonth)
monthAvTemp= matrix(NA, 12, 4, dimnames=list(
                                   month.abb, fig1.11Stns))
for(i in 1:4)
  monthAvTemp[, i] = tapply(fig1.11Temp[, i], monthIndex, mean)

StnLtrs = substring(fig1.11Stns, 1, 1)
matpoints(monthMid, monthAvTemp, pch=StnLtrs, lwd=2, col=1)

legend('bottom', paste(fig1.11Stns, ' (', StnLtrs, ')', sep=''),
       lty=1:4, col=1, lwd=2)

# Figure 1.12

dayrange     = c(0,365)
Lbasis       = create.constant.basis(dayrange,
                                  axes=list("axesIntervals"))
Lcoef        = matrix(c(0,(2*pi/365)^2,0),1,3)
bfdobj       = fd(Lcoef,Lbasis)
bwtlist      = fd2list(bfdobj)
harmaccelLfd = Lfd(3, bwtlist)

Lfd.Temp = deriv(fig1.11Temp.fd, harmaccelLfd)

plot(Lfd.Temp, day.5, axes=FALSE, xlab='', ylab='L-Temperature',
     col=1, lwd=2)
axis(2, las=1)
axisIntervals(labels=monthLetters)
legend('bottom', paste(fig1.11Stns, ' (', StnLtrs, ')', sep=''),
       lty=1:4, lwd=2)

##
## Section 1.5.  The first steps in a functional data analysis
##

# Figure 1.13

P.RupertPrecip = CanadianWeather$dailyAv[ ,
                        'Pr. Rupert', 'Precipitation.mm']

plot(day.5, P.RupertPrecip, axes=FALSE, pch='.', cex=2,
     xlab="", ylab="Preciipitation (mm)")
axis(2, las=1)
axisIntervals(labels=monthLetters)

TempBasis = create.fourier.basis(c(0, 365), 13)
P.Rupert.Prec.fd = smooth.basisPar(day.5, P.RupertPrecip,
                 TempBasis, harmaccelLfd, lambda=10^7)
lines(P.Rupert.Prec.fd, lwd=2)

# Figure 1.14

matplot(pinchtime, pinchraw, type='l', xlab='seconds',
        ylab='Force (N)')

# Figure 1.15

(i11.7 = which(abs(agefine-11.7) == min(abs(agefine-11.7)))[1])

# Use the fit from Figure 1.1

hgtf.vel1.5 = predict(hgtfmonfd$yhatfd, agefine, 1)
hgtf.acc1.5 = predict(hgtfmonfd$yhatfd, agefine, 2)

plot(hgtf.vel1.5, hgtf.acc1.5, type='n', xlim=c(0, 12), ylim=c(-5, 2),
     xlab='Velocity (cm/yr)', ylab=expression(Acceleration (cm/yr^2)),
     las=1)
for(i in 1:10){
  lines(hgtf.vel1.5[, i], hgtf.acc1.5[, i])
  points(hgtf.vel1.5[i11.7, i], hgtf.acc1.5[i11.7, i])
}
abline(h=0, lty='dotted')

##
## Section 1.6.  Exploring variability in functional data
##
# (no data analysis in this section)

##
## Section 1.7.  Functional linear models
##
# (no data analysis in this section)

##
## Section 1.8.  Using derivatives in functional data analysis
##
# (no data analysis in this section)

##
## Section 1.9.  Concluding remarks
##
# (no data analysis in this section)

##
## Section 1.10.  Some Things to Try
##
