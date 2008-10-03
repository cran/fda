###
###
### Ramsey, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 5.  Smoothing: Computing Curves from Noisy Data
###
library(fda)

heightmat = growth$hgtf
##
## Section 5.1.  Regression Splines: Smoothing by Regression Analysis
##
heightbasis12 = create.bspline.basis(c(1,18), 12, 6)

basismat = eval.basis(growth$age, heightbasis12)

heightcoef = lsfit(basismat, heightmat,
                   intercept=FALSE)$coef

heightList = smooth.basis(growth$age, heightmat,
                          heightbasis12)
heightfd   = heightList$fd
height.df  = heightList$df
height.gcv = heightList$gcv

age           = growth$age
heightbasismat= eval.basis(age, heightbasis12)
y2cMap        = solve(crossprod(heightbasismat),
                       t(heightbasismat))

##
## Section 5.2.  Data Smoothing with Roughness Penalties
##

# section 5.2.2 The Roughness Penalty Matrix R
tempbasis = create.fourier.basis(c(0,365),65)
harmaccelLfd = vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
Rmat = eval.penalty(tempbasis, harmaccelLfd)

# section 5.2.4 Defining Smoothing by Functional Parameter Objects

norder     = 6
nbasis     = length(age) + norder - 2
heightbasis= create.bspline.basis(c(1,18), nbasis, norder, age)

heightfdPar   = fdPar(heightbasis, 4, 0.01)

heightfdSmooth= smooth.basis(age, heightmat, heightfdPar)
heightfd      = heightfdSmooth$fd

# section 5.2.5 Choosing Smoothing Parameter lambda

loglam        = seq(-6, 0, 0.25)
Gcvsave       = rep(NA, length(loglam))
names(Gcvsave)= loglam
Dfsave  = Gcvsave
for(i in 1:length(loglam)){
  hgtfdPari  = fdPar(heightbasis, Lfdobj=4, 10^loglam[i])
#          penalize curvature of acceleration (Lfdojb=4)
  hgtSm.i    = smooth.basis(age, heightmat, hgtfdPari)
  Gcvsave[i] = sum(hgtSm.i$gcv)
  Dfsave[i]  = hgtSm.i$df
}

# Figure 5.1.
plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )

##
## 5.3.  Case Study: The Log Precipitation Data
##
logprecav = CanadianWeather$dailyAv[
         dayOfYearShifted, , 'log10precip']

dayrange  = c(0,365)
daybasis  = create.fourier.basis(dayrange, 365)

Lcoef        = c(0,(2*pi/diff(dayrange))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, dayrange)

loglam       = seq(4,9,0.25)
nlam         = length(loglam)
dfsave       = rep(NA,nlam)
names(dfsave)= loglam
gcvsave      = dfsave

for (ilam in 1:nlam) {
  cat(paste('log10 lambda =',loglam[ilam],'\n'))
  lambda       = 10^loglam[ilam]
  fdParobj     = fdPar(daybasis, harmaccelLfd, lambda)
  smoothlist   = smooth.basis(day.5, logprecav,
                            fdParobj)
  dfsave[ilam] = smoothlist$df
  gcvsave[ilam]= sum(smoothlist$gcv)
}

# Figure 5.2.
plot(loglam, gcvsave, type='b', ylab='GCV Criterion',
     xlab=expression(log[10](lambda)) )

lambda   = 1e6
fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd = logprec.fit$fd
fdnames = list("Day (July 1 to June 30)",
               "Weather Station" = CanadianWeather$place,
               "Log 10 Precipitation (mm)")
logprec.fd$fdnames = fdnames
plot(logprec.fd)

# plotfit.fd:  Pauses between plots
# *** --->>> input required (e.g., click on the plot)
#            to advance to the next plot
plotfit.fd(logprecav, day.5, logprec.fd)

##
## Section 5.4 Positive, Monotone, Density
##             and Other Constrained Functions
##

# sec. 5.4.1.  Positive Smoothing
lambda     = 1e3
WfdParobj  = fdPar(daybasis, harmaccelLfd, lambda)
VanPrec    = CanadianWeather$dailyAv[
  dayOfYearShifted, 'Vancouver', 'Precipitation.mm']
VanPrecPos = smooth.pos(day.5, VanPrec, WfdParobj)
Wfd        = VanPrecPos$Wfdobj
Wfd$fdnames= list("Day (July 1 to June 30)",
      "Weather Station" = CanadianWeather$place,
                   "Log 10 Precipitation (mm)")

precfit = exp(eval.fd(day.5, Wfd))
plot(day.5, VanPrec, type="p", cex=1.2,
     xlab="Day (July 1 to June 30)",
     ylab="Millimeters",
     main="Vancouver's Precipitation")
lines(day.5, precfit,lwd=2)

# 5.4.2.  Monotone smoothing

day    = infantGrowth[, 'day']
tib    = infantGrowth[, 'tibiaLength']
n      = length(tib)
nbasis = 42

Wbasis = create.bspline.basis(c(1,n), nbasis)
Wfd0   = fd(matrix(0,nbasis,1), Wbasis)
WfdPar = fdPar(Wfd0, 2, 1e-4)
result = smooth.monotone(day, tib, WfdPar)
Wfd    = result$Wfd
beta   = result$beta
dayfine = seq(1,n,len=151)
tibhat  = beta[1]+beta[2]*eval.monfd(dayfine ,Wfd)
Dtibhat =        beta[2]*eval.monfd(dayfine, Wfd, 1)
D2tibhat=        beta[2]*eval.monfd(dayfine, Wfd, 2)

op <- par(mfrow=c(1,2), mar=c(5,5,3,2), lwd=2)

#  plot height

plot(day, tib, type = "p", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Length (mm)")
lines(dayfine, tibhat, lwd=2)

#  plot velocity

plot(dayfine, Dtibhat, type = "l", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Velocity (mm/day)")

#  plot acceleration

plot(dayfine, D2tibhat, type = "l", cex=1.2, las=1,
     xlab="Day", ylab="Tibia Acceleration (mm/day/day)")
lines(c(1,n),c(0,0),lty=2)

par(op)

# 5.4.3.  Probability Density Functions
NR <- length(ReginaPrecip)
plot(1:NR, sort(ReginaPrecip), xlab='Rank of rainfall',
     ylab='Ordered daily rainfall (mm)' )

sel2.45 <- ((2 <= ReginaPrecip) & (ReginaPrecip <= 45))
RegPrec <- sort(ReginaPrecip[sel2.45])
N <- length(RegPrec)

Wknots  = RegPrec[round(N*seq(1/N,1,len=11),0)]
Wnbasis = length(Wknots) + 2
Wbasis  = create.bspline.basis(range(RegPrec),13,4,Wknots)

Wlambda     = 1e-1
WfdPar      = fdPar(Wbasis, 2, Wlambda)
densityList = density.fd(RegPrec, WfdPar)
Wfd         = densityList$Wfdobj
C.           = densityList$C

Zfine = seq(RegPrec[1],RegPrec[N],len=201)
Wfine = eval.fd(Zfine, Wfd)
Pfine = exp(Wfine)/C.

plot(Zfine, Pfine, type='l', xlab='Precipitation (mm)',
     ylab='Probability Density')
abline(v=Wknots, lty='dashed')

##
## Section 5.5 Assessing the Fit to the Data
##

logprecmat = eval.fd(day.5, logprec.fd)
logprecres = logprecav - logprecmat
#  across stations
logprecvar1 = apply(logprecres^2, 1, sum)/35
#  across time
logprecvar2 = apply(logprecres^2, 2, sum)/(365-12)

# Figure 5.7
plot(sqrt(logprecvar2), xlab='Station Number',
     ylab='Standard Deviation across Day')
rt <- which(CanadianWeather$place %in%
     c("Winnipeg", 'Regina', 'Churchill', 'Montreal', 'St. Johns'))
lft <- which(CanadianWeather$place %in%
             c('Yellowknife', 'Resolute', 'Vancouver', 'Iqaluit',
               'Pr. George', 'Pr. Rupert') )
below <- which(CanadianWeather$place %in% 'Edmonton')
top <- which(CanadianWeather$place %in% 'Halifax')

text(rt, sqrt(logprecvar2[rt]), labels=CanadianWeather$place[rt],
     pos=4)
text(lft, sqrt(logprecvar2[lft]), labels=CanadianWeather$place[lft],
     pos=2)
text(below, sqrt(logprecvar2[below]), labels=CanadianWeather$place[below],
     pos=1)
text(top, sqrt(logprecvar2[top]), labels=CanadianWeather$place[top],
     pos=3)

# Figure 5.8
logstddev.fit = smooth.basis(day.5, log(logprecvar1)/2, fdParobj)
logstddev.fd = logstddev.fit$fd
logprecvar1fit = exp(eval.fd(day.5, logstddev.fd))

plot(day.5, sqrt(logprecvar1), xlab='Day',
     ylab='Standard seviation across stations')
lines(day.5, logprecvar1fit, lwd=2)

##
## Section 5.6 Details for the fdPar Class and smooth.basis Function
##
help(fdPar)
help(smooth.basis)

##
## Section 5.8 Some Things to Try
##
# (exercises for the reader)

##
## Section 5.7 More to Read
##
