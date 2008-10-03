###
###
### Ramsey, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. Chapter 7  Exploring Variation: Functional Principal
###                and Canonical Components analysis
###
library(fda)

##
## Section 7.1 An Overview of Functional PCA
##
#  (no computations in this section)

##
## Section 7.2 PCA with Function pca.fd
##
#  Create logprec.fd, copy from chapter 6

# Section 7.2.1 PCA of the Log Precipitation Data
logprecav = CanadianWeather$dailyAv[
         dayOfYearShifted, , 'log10precip']
dayrange  = c(0,365)
daybasis  = create.fourier.basis(dayrange, 365)

Lcoef        = c(0,(2*pi/diff(dayrange))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, dayrange)

lambda   = 1e6
fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd = logprec.fit$fd

logprec.pcalist = pca.fd(logprec.fd, 2)
print(logprec.pcalist$values)

# Figure 7.1
plot.pca.fd(logprec.pcalist)
plot(logprec.pcalist, expand=.5)

# Figure 7.2
logprec.rotpcalist = varmx.pca.fd(logprec.pcalist)
plot.pca.fd(logprec.rotpcalist, expand=.5)

# Figure 7.3

#  plot.pca.fd(..., type='scores')???





# Section 7.2.2 PCA of Log Precipitation Residuals
# logprecres = residuals from
# the smooths of the log precipitation curves in Chapter 5.

logprecav = CanadianWeather$dailyAv[
         dayOfYearShifted, , 'log10precip']

dayrange  = c(0,365)
daybasis  = create.fourier.basis(dayrange, 365)
Lcoef        = c(0,(2*pi/diff(dayrange))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, dayrange)
lambda   = 1e6
fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd = logprec.fit$fd
fdnames = list("Day (July 1 to June 30)",
               "Weather Station" = CanadianWeather$place,
               "Log 10 Precipitation (mm)")
logprec.fd$fdnames = fdnames

logprecmat = eval.fd(day.5, logprec.fd)
logprecres = logprecav - logprecmat

# Figure 7.4
logprecres.fd = smooth.basis(day.5, logprecres,
    fdParobj)$fd
plot(logprecres.fd, lwd=2, col=4, lty=1, cex=1.2,
     xlim=c(0,365), ylim=c(-0.07, 0.07),
     xlab="Day", ylab="Residual (log 10 mm)")

# Figure 7.5
logprec.pca1 = pca.fd(logprecres.fd, 1)
plot(logprec.pca1, expand=0.01)
# ???????????????


##
## Section 7.3 More Functional PCA Features
##
#  (no computations in this section)

##
## Section 7.4 PCA of joint X-Y Variation in Handwriting
##
nharm = 3

# Need 'fdafd' created in Section 1.2
fdapcaList = pca.fd(fdafd, nharm)
plot.pca.fd(fdapcaList)
fdarotpcaList = varmx.pca.fd(fdapcaList)
plot.pca.fd(fdarotpcaList)

fdaeig = fdapcaList$values
neig = 12
x = matrix(1,neig-nharm,2)
x[,2] = nharm:neig
y = log10(fdaeig[nharm:neig])
c = lsfit(x,y,int=FALSE)$coef

# Figure 7.6
op <- par(mfrow=c(1,1),cex=1.2)
plot(1:neig, log10(fdaeig[1:neig]), "b",
     xlab="Eigenvalue Number",
     ylab="Log10 Eigenvalue")
lines(1:neig, c[1]+ c[2]*(1:neig), lty=2)
par(op)

# Figure 7.7 varimax rotation ...







##
## Section 7.5 Exploring Functional Covariation
##             with Canonical Correlation Analysis
##

ccabasis = create.fourier.basis(dayrange, 3)

#  need tempfd ???

ccalist = cca.fd(tempfd, logprecfd, 3, ccabasis, ccabasis)
#  Figure 7.8 & 7.9







##
## Section 7.6 Details for the pca.fd and cca.fd Functions
##
help(pca.fd)
help(cca.fd)

##
## Section 7.7 Some Things to Try
##
# (exercises for the reader)

##
## Section 7.8 More to Read
##
