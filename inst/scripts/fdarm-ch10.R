###
###
### Ramsey, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch.  10.  Linear Models for Functional Responses
###
library(fda)

##
## Section 10.1  Functional Responses and an Analysis of Variance Model
##
#  Section 10.1.1 Climate Region Effects on Temperature
regions = unique(CanadianWeather$region)
p = length(regions) + 1
regionList = vector("list", p)
regionList[[1]] = c(rep(1,35),0)
for (j in 2:p) {
  xj = CanadianWeather$region == regions[j-1]
  regionList[[j]] = c(xj,1)
}

# tempfd from chapter 9


Lcoef = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))
daybasis365=create.fourier.basis(c(0, 365), 365)
lambda =1e6
tempfdPar365 =fdPar(daybasis365, harmaccelLfd, lambda)
tempSmooth365 <- smooth.basis(day.5, daily$tempav,
                              tempfdPar365)
tempfd = tempSmooth365$fd



# None of the code in the rest of section 10.1.1 has been tested

coef = tempfd$coef
coef36 = cbind(coef,matrix(0,65,1))
temp36fd = fd(coef36,tempbasis,tempfd$fdnames)

betabasis = create.fourier.basis(c(0, 365), 11)
betafdPar = fdPar(betabasis)
betaList = vector("list",p)
for (j in 1:p) betaList[[j]] = betafdPar

fRegressList = fRegress(temp36fd, regionList,
betaList)
betaestList = fRegressList$betaestlist
regionFit = fRegressList$yhatfd
regions = c("Canada", regions)
par(mfrow=c(2,3),cex=1)
for (j in 1:p) plot(betaestList[[j]]$fd, lwd=2,
                    xlab="Day (July 1 to June 30)",
                    ylab="", main=regions[j])
plot(regionFit, lwd=2, col=1, lty=1,
     xlab="Day", ylab="", main="Prediction")




# 10.1.2 Trends in Sea Bird Populations on Kodiak Island

# Need the data

# Figure 10.2

fooddummy = matrix(0,13,1)
foodindex = c(1,2,5,6,12,13)
fooddummy[foodindex] = 1
fooddummy = rbind(fooddummy, fooddummy)
birddummy = diag(rep(1,13))

Zmat = matrix(0,28,15)
Zmat[1:26,1] = rep(1,26)
Zmat[1:26,2] = fooddummy
Zmat[ 1:13,3:15] = birddummy
Zmat[14:26,3:15] = birddummy
Zmat[27, foodindex+2 ] = 1
Zmat[28,-(foodindex+2)] = 1

logbirdcoef = logbirdfd$coefs
logbirdcoef0 = cbind(logbirdcoef,
matrix(0,nbirdbasis,2))
logbirdfd0 = fd(logbirdcoef0,birdbasis)

p = 15
xfdlist = vector("list",p)
for (j in 1:p) xfdlist[[j]] = Zmat[,j]
betalist = vector("list",p)
foodbasis = create.bspline.basis(c(0,19),5)
betalist[[1]] = fdPar(foodbasis)
betalist[[2]] = fdPar(foodbasis)
birdbasis = create.constant.basis(c(0,19))
for (j in 3:p) betalist[[j]] = fdPar(birdbasis)

fRegressList = fRegress(logbirdfd0,xfdlist,betalist)
betaestlist = fRegressList$betaestlist
yhatfdobj = fRegressList$yhatfdobj

# Figure 10.3 ???




# Section 10.1.3 Choosing Smoothing Parameters

loglam = seq(-2,0,0.25)
SSE.CV = rep(0,length(loglam))
betafdPari = betafdPar
for(i in 1:length(loglam)){
  betafdPari$lambda = 10^loglam[i]
  betalisti = betalist
  for (j in 1:2) betalisti[[j]] <- betafdPari

# *** Need logbirdfd0

  SSE.CV[i] <- fRegress.CV(logbirdfd0, xfdlist,
                           betalisti,CVobs=1:26)$SSE.CV
}

# Figure 10.4



##
## Section 10.2 Functional Responses with Functional Predictors:
##              The Concurrent Model
##
#  Section 10.2.1 Estimation for the Concurrent Model

#  Section 10.2.2 Confidence Intervals for Regression Functions

yhatmat = eval.fd(plotyear, yhatfdobj)
ymat = eval.fd(plotyear, logbirdfd0)

rmat = ymat - yhatmat
SigmaE = var(t(rmat))
stddevE = sqrt(diag(SigmaE))
SigmaE = diag(stddevEˆ2)

birdbasismat = eval.basis(plotyear, birdbasis)
y2cMap = solve(crossprod(birdbasismat)), t(birdbasismat))

stderrList = fRegress.stderr(fRegressList, y2cMap,
     SigmaE)
betastderrlist = stderrList$betastderrlist

par(mfrow=c(2,1),ask=FALSE)
titlelist = vector("list", p)
titlelist[[1]] = "Intercept"
titlelist[[2]] = "Feed effect"
plotbeta(betaestlist, betastderrlist,
        titlelist=titlelist, index=1:2)

# Section 10.2.3 Knee Angle Predicted from Hip Angle

xfdlist = list(rep(1,39), hipfd)

# *** Where's hipfd ???

betafdPar = fdPar(gaitbasis, harmaccelLfd)
betalist = list(betafdPar,betafdPar)
fRegressList = fRegress(kneefd, xfdlist, betalist)
kneehatfd = fRegressList$yhatfd
betaestlist = fRegressList$betaestlist

kneemat = eval.fd(gaittime, kneefd)

kneehatmat = eval.fd(gaittime, kneehatfd)
resmat = kneemat - kneehatmat
SigmaE = cov(t(resmat))

kneefinemat = eval.fd(gaitfine, kneefd)
kneemeanvec = eval.fd(gaitfine, kneemeanfd)
kneehatfinemat = eval.fd(gaitfine, kneehatfd)
resmat = kneefinemat - kneehatfinemat
resmat0 = kneefinemat -
kneemeanvec %*% matrix(1,1,ncurve)
SSE0 = apply((resmat0)ˆ2, 1, sum)
SSE1 = apply(resmatˆ2, 1, sum)
Rsqr = (SSE0-SSE1)/SSE0

fRegressList1 = fRegress(kneefd, xfdlist, betalist,
           y2cMap, SigmaE)

fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
betastderrlist = fRegressList2$betastderrlist
titlelist = list("Intercept", "Hip coefficient")
plotbeta(betaestlist, betastderrlist, gaitfine, titlelist)

##
## Section 10.3 Beyond the Concurrent Model
##
#  (no computations in this section)

##
## Section 10.4 A Functional Linear Model for Swedish Mortality
##
betabasis = create.bspline.basis(c(0,80),23)
beta0Par = fdPar(betabasis, 2, 1e-5)
beta1sPar = fdPar(betabasis, 2, 1e3)
beta1tPar = fdPar(betabasis, 2, 1e3)
betaList = list(beta0Par, beta1sPar, beta1tPar)

linmodSmooth = linmod(NextYear, LastYear, betaList)

# Where's LastYear?  ???

# Figure 10.11






##
## Section 10.5 Permutation Tests of Functional Hypotheses
##
#  Section 10.5.1 Functional t-Tests

tperm.fd(hgtmfd,hgtffd)

# Figure 10.12

# Section 10.5.2 Functional F-Tests

F.res = Fperm.fd(temp36fd, regionlist, betaList)









##
## 10.6 Details for R Functions fRegress, fRegress.CV and fRegress.stderr
##
help(fRegress)
help(fRegress.CV)
help(fRegress.stderr)

##
## 10.7 Details for Function plotbeta
##
help(plotbeta)

##
## 10.8 Details for Function linmod
##
help(linmod)

##
## 10.9 Details for Functions Fperm.fd and tperm.fd
##
help(Fperm.fd)
help(tperm.fd)

##
## Section 10.10 Some Things to Try
##
# (exercises for the reader)

##
## Section 10.11  More to Read
##
