###
###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch. 9.  Functional Linear Models for Scalar Response
###

library(fda)

##
## Section 9.1 Functional Linear regression with a Scalar response
##

#  (no computations in this section)

##
## Section 9.2 A Scalar Response Model for Log Annual Precipitation
##

annualprec   = log10(apply(daily$precav,2,sum))

tempbasis65  = create.fourier.basis(c(0,365),65)
tempSmooth65 = smooth.basis(day.5, daily$tempav, tempbasis65)
tempfd65     = tempSmooth65$fd

##
## Section 9.3 Setting Up the Functional Linear Model
##
#  (no computations in this section)

##
## Section 9.4 Three Estimates of the Regression Coefficient
##             Predicting Annual Precipitation
##

templist      = vector("list",2)
templist[[1]] = rep(1,35)
templist[[2]] = tempfd65

# 9.4.1 Low Dimensional Regression Coefficient Function beta

conbasis   = create.constant.basis(c(0,365))
betabasis5 = create.fourier.basis(c(0,365),5)
betalist1  = vector("list",2)
betalist1[[1]] = fdPar(conbasis)
betalist1[[2]] = fdPar(betabasis5)

fRegressList1 = fRegress(annualprec,templist,betalist1)

betaestlist1  = fRegressList1$betaestlist
tempbetafd1   = betaestlist1[[2]]$fd

# Figure 9.1

plot(tempbetafd1, xlab="Day", ylab="Beta for temperature")

coef(betaestlist1[[1]])
# 0.0095 as in the book

annualprechat1 = fRegressList1$yhatfdobj
annualprecres1 = annualprec - annualprechat1
SSE1.1  = sum(annualprecres1^2)
SSE0    = sum((annualprec - mean(annualprec))^2)
(RSQ1   = (SSE0-SSE1.1)/SSE0)
# 0.80 as in the book
(Fratio1 = ((SSE0-SSE1.1)/5)/(SSE1.1/29))
# 22.6 as in the book

# 9.4.2 Coefficient beta Estimate Using a Roughness Penalty

Lcoef = c(0,(2*pi/365)^2,0)

harmaccelLfd = vec2Lfd(Lcoef, c(0,365))

# refit with 35 terms rather than 5 in the fourier basis
betabasis35 = create.fourier.basis(c(0, 365), 35)
lambda      = 10^12.5
betafdPar.  = fdPar(betabasis35, harmaccelLfd, lambda)

betalist2      = betalist1
betalist2[[2]] = betafdPar.

annPrecTemp    = fRegress(annualprec, templist, betalist2)
betaestlist2   = annPrecTemp$betaestlist
annualprechat2 = annPrecTemp$yhatfdobj

print(annPrecTemp$df)

SSE1.2 = sum((annualprec-annualprechat2)^2)
(RSQ2 = (SSE0 - SSE1.2)/SSE0)
# 0.75 as in the book

(Fratio2 = ((SSE0-SSE1.2)/3.7)/(SSE1.2/30.3))
# 25.1 as in the book

# Figure 9.2

plot(annualprechat2, annualprec)
abline(lm(annualprec~annualprechat2), lty='dashed')

# Figure 9.3
# ... see section 9.4.4 below ...

plot(betaestlist2[[2]]$fd)

# Compare with the constant fit:

betalist      = betalist1
betalist[[2]] = fdPar(conbasis)
fRegressList  = fRegress(annualprec, templist, betalist)
betaestlist   = fRegressList$betaestlist

annualprechat = fRegressList$yhatfdobj
SSE1 = sum((annualprec-annualprechat)^2)

(RSQ = (SSE0 - SSE1)/SSE0)
# 0.49 as in the book

(Fratio = ((SSE0-SSE1)/1)/(SSE1/33))
# 31.3 as in the book


# 9.4.3 Choosing Smoothing Parameters

loglam = seq(5,15,0.5)
nlam   = length(loglam)
SSE.CV = rep(NA,nlam)
for (ilam in 1:nlam) {
  lambda     = 10^(loglam[ilam])
  betalisti  = betalist2
  betafdPar2 = betalisti[[2]]
  betafdPar2$lambda = lambda
  betalisti[[2]] = betafdPar2
  fRegi          = fRegress.CV(annualprec, templist, betalisti)
  SSE.CV[ilam]   = fRegi$SSE.CV
}

plot(loglam, SSE.CV, type="b", lwd=2,
     xlab="log smoothing parameter lambda",
     ylab="Cross-validation score", cex=1.2)

# 9.4.4 Confidence Intervals

resid   = annualprec - annualprechat2
SigmaE. = sum(resid^2)/(35-annPrecTemp$df)
SigmaE  = SigmaE.*diag(rep(1,35))
y2cMap  = tempSmooth65$y2cMap

stderrList = fRegress.stderr(annPrecTemp, y2cMap, SigmaE)

betafdPar      = betaestlist2[[2]]
betafd         = betafdPar$fd
betastderrList = stderrList$betastderrlist
betastderrfd   = betastderrList[[2]]

# Figure 9.3

plot(betafd, xlab="Day", ylab="Temperature Reg. Coeff.",
     ylim=c(-6e-4,1.2e-03), lwd=2)

lines(betafd+2*betastderrfd, lty=2, lwd=1)
lines(betafd-2*betastderrfd, lty=2, lwd=1)

# Section 9.4.5 Scalar Response Models by Functional Principal Components

daybasis365   = create.fourier.basis(c(0, 365), 365)
lambda        = 1e6
tempfdPar365  = fdPar(daybasis365, harmaccelLfd, lambda)
tempSmooth365 = smooth.basis(day.5, daily$tempav,
                              tempfdPar365)
tempfd = tempSmooth365$fd

lambda    = 1e0
tempfdPar = fdPar(daybasis365, harmaccelLfd, lambda)
temppca   = pca.fd(tempfd, 4, tempfdPar)
harmonics = temppca$harmonics

pcamodel = lm(annualprec~temppca$scores)
pcacoefs = summary(pcamodel)$coef
betafd   = pcacoefs[2,1]*harmonics[1] + pcacoefs[3,1]*harmonics[2] +
           pcacoefs[4,1]*harmonics[3]
coefvar  = pcacoefs[,2]^2
betavar  = coefvar[2]*harmonics[1]^2 + coefvar[3]*harmonics[2]^2 +
           coefvar[4]*harmonics[3]^2

# Figure 9.5

plot(betafd, xlab="Day", ylab="Regression Coef.",
     ylim=c(-6e-4,1.2e-03), lwd=2)
lines(betafd+2*sqrt(betavar), lty=2, lwd=1)
lines(betafd-2*sqrt(betavar), lty=2, lwd=1)

##
## Section 9.5 Statistical Tests
##

F.res = Fperm.fd(annualprec, templist, betalist)
F.res$Fobs

F.res$qval

##
## Section 9.6 Some Things to Try
##
# (exercises for the reader)

##
## Section 9.7  More to Read
##
