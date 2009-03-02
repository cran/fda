###
###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###
### ch.  10.  Linear Models for Functional Responses
###

library(fda)

##
## Section 10.1  Functional Responses and an Analysis of Variance Model
##

####### TEMPERATURE ~ REGION #######

#  Section 10.1.1 Climate Region Effects on Temperature

regions.         = unique(CanadianWeather$region)
p                = length(regions.) + 1
regionList       = vector("list", p)
names(regionList)= c('Canada', regions.)
regionList[[1]]  = c(rep(1,35),0)
for (j in 2:p) {
  xj             = (CanadianWeather$region == regions.[j-1])
  regionList[[j]]= c(xj,1)
}

# tempfd from chapter 9

Lcoef       = c(0,(2*pi/365)^2,0)
harmaccelLfd= vec2Lfd(Lcoef, c(0,365))
tempbasis   = create.fourier.basis(c(0, 365), 65)
lambda      = 1e6
tempfdPar65 = fdPar(tempbasis, harmaccelLfd, lambda)
tempShifted = daily$tempav[dayOfYearShifted, ]
tempSmooth65= smooth.basis(day.5, tempShifted, tempfdPar65)
tempfd      = tempSmooth65$fd

coef    = tempfd$coef
coef36  = cbind(coef,matrix(0,65,1))
temp36fd= fd(coef36,tempbasis,tempfd$fdnames)

betabasis      = create.fourier.basis(c(0, 365), 11)
betafdPar      = fdPar(betabasis)
betaList       = vector("list",p)
names(betaList)= regions.
for (j in 1:p) betaList[[j]] = betafdPar

fRegressList= fRegress(temp36fd, regionList, betaList)
betaestList = fRegressList$betaestlist
regionFit   = fRegressList$yhatfd
regions     = c("Canada", regions.)

# Figure 10.1

op          = par(mfrow=c(2,3),cex=1)
for (j in 1:p) plot(betaestList[[j]]$fd, lwd=2,
                    xlab="Day (July 1 to June 30)",
                    ylab="", main=regions[j])
plot(regionFit, lwd=2, col=1, lty=1,
     xlab="Day (July 1 to June 30)", ylab="", main="Prediction")
par(op)

####### SEA BIRDS ON KODIAK ISLAND #######

# 10.1.2 Trends in Sea Bird Populations on Kodiak Island

#  select only the data for sites Uyak and Uganik, which have data
#  from 1986 to 2005, except for 1998

sites = c('Uganik', 'Uyak')
sel   = seabird$Bay %in% sites
UU    = seabird[sel,]

# Drop 2 species with many NAs

NAs  = sapply(UU, function(x)sum(is.na(x)))
NAs. = which(NAs > 2)
birdindex= (1:15)[-NAs.]
birds = names(UU)[birdindex]

#  Compute mean counts taken over both sites and transects

meanCounts = matrix(NA, 20, 13)
dimnames(meanCounts) = list(1986:2005, birds)

for(i in 1:20){
  sel = (UU$Year == rownames(meanCounts)[i])
  meanCounts[i, ] = sapply(UU[sel, birds], mean, na.rm=TRUE)
}

selYear = !is.na(meanCounts[,1])
logCounts = log10(meanCounts[selYear,])

#  time vectors in years and in indices in 1:20

yearObs   = as.numeric(rownames(logCounts))
yearCode  = (1:20)[selYear]

# Figure 10.2

shellfishindex = c(1,2,5,6,12,13)
fishindex      = (1:13)[-shellfishindex]
ylim = range(logCounts)

op = par(mfrow=c(2,1), mar=c(2, 4, 4, 1)+.1)

matplot(yearObs, logCounts[, shellfishindex], xlab='', ylab='',
        ylim=ylim, main='Shellfish Diet', type='b', col=1)
meanShellfish = apply(meanCounts[, shellfishindex], 1, mean)
lines(yearObs, log10(meanShellfish[!is.na(meanShellfish)]), lwd=3)
abline(h=0, lty='dotted')

matplot(yearObs, logCounts[, fishindex], xlab='', ylab='',
        ylim=ylim, main='Fish Diet', type='b', col=1)
meanFish = apply(meanCounts[, shellfishindex], 1, mean)
lines(yearObs, log10(meanFish[!is.na(meanFish)]), lwd=3)
abline(h=0, lty='dotted')
par(op)

#  Compute mean counts taken over transects only within sites
#  so we have 2 observations for each bird species each year.
#  Two of these counts are zero, and are replaced by 1/(2*n)

meanCounts2 = matrix(NA, 20, 26)

for(i in 1:20) for (j in 1:2) {
  sel = (UU$Year == rownames(meanCounts)[i] & as.character(UU$Bay) == sites[j])
  meanCountsij = sapply(UU[sel, birds], mean, na.rm=TRUE)
  n = sum(sel)
  if (n > 0) {
    meanCountsij[meanCountsij == 0] = 1/(2*n)
  }
  meanCounts2[i,(j-1)*13+(1:13)] = meanCountsij
}

selYear2   = !is.na(meanCounts2[, 1])
yearCode  = (1:20)[selYear2]
all.equal(yearCode, c(1:12, 14:20))

logCounts2 = log10(meanCounts2[selYear2,])

#  Represent log mean counts exactly with a polygonal basis

birdbasis = create.polygonal.basis(yearCode)
#birdlist = smooth.basis(yearCode, logCounts, birdbasis)
birdlist2 = smooth.basis(yearCode, logCounts2, birdbasis)

#birdfd  = birdlist$fd
birdfd2 = birdlist2$fd

#yearfine = seq(1, 20, len=191)
#birdmatS = eval.fd(yearfine, birdfd[ shellfishindex])
#birdmatF = eval.fd(yearfine, birdfd[-shellfishindex])
#birdvecS = apply(birdmatS, 1, mean)
#birdvecF = apply(birdmatF, 1, mean)

#  -----------------------------------------------------------------
#  After some preliminary analyses we determined that there was no
#  contribution from either site or food*site interaction.
#  Now we use a reduced model with only a feed effect,
#  but we add bird effects, which were seen in the plot to be
#  strong.  Birds are nested within feed groups, and either their
#  effects must sum to zero within each group, or we must designate
#  a bird in each group as a baseline, and provide dummy variables
#  for the remainder.  We opt for the latter strategy.
#  -----------------------------------------------------------------

#  The design matrix contains an intercept dummy variable, a
#  feed dummy variable, and dummy variables for birds, excluding
#  the second bird in each group, which turns out to be the each
#  group's most abundant species, and which is designated as the
#  baseline bird for that group.

#meanlogCounts = apply(logCounts,2,mean)
#meanlogCounts[shellfishindex]
#meanlogCounts[-shellfishindex]

# 15 columns for the intercept + diet + 13 bird species
# 26 rows for the 26 (species - bay) combinations
Zmat0 = matrix(0,26,15)

#  Intercept or baseline effect

Intercept = rep(1,26)

#  Crustacean/Mollusc feeding effect:  a contrast between the two groups

foodindex = c(1,2,5,6,12,13)
fooddummy = c(2*rep(1:13 %in% foodindex, 2)-1)

#  Bird effect, one for each species

birddummy = diag(rep(1,13))
birdvarbl = rbind(birddummy,birddummy)

#  fill the columns of the design matrix

Zmat0[,1]    = Intercept
Zmat0[,2]    = fooddummy
Zmat0[,3:15] = birdvarbl

#  Two extra dummy observations are added to the functional data
#  object for log counts, and two additional rows are added to
#  the design matrix to force the bird effects within each diet
#  group to equal 0.

birdfd3 = birdfd2
birdfd3$coefs = cbind(birdfd3$coefs, matrix(0,19,2))

Zmat = rbind(Zmat0, matrix(0,2,15))
Zmat[27,shellfishindex+2]  = 1
Zmat[28,fishindex+2] = 1

p = 15
xfdlist = vector("list",p)
names(xfdlist) = c("const", "diet", birds)
betalist = xfdlist
for (j in 1:p) xfdlist[[j]] = Zmat[,j]

#  set up the functional parameter object for (the regression fns.
#  use cubic b-spline basis for intercept and food coefficients

betabasis1 = create.bspline.basis(c(1,20),21,4,yearCode)
lambda = 10
betafdPar1 = fdPar(betabasis1,2,lambda)
betalist[[1]] = betafdPar1
betalist[[2]] = betafdPar1
betabasis2 = create.constant.basis(c(1,20))
betafdPar2 = fdPar(betabasis2)
for (j in 3:15) betalist[[j]] = betafdPar2

birdRegress = fRegress(birdfd3, xfdlist, betalist)
betaestlist = birdRegress$betaestlist

# Figure 10.3 is produced in Section 10.2.2 below
# after estimating the smoothing parameter in Section 10.1.3
#
# Here we plot the regression parameters
# without the confidence intervals.

op = par(mfrow=c(2,1))
plot(betaestlist$const$fd)
plot(betaestlist$diet$fd)
par(op)

##
## Section 10.1.3 Choosing Smoothing Parameters
##

#  Choose the level of smoothing by minimizing cross-validated
#  error sums of squares.

loglam = seq(-2,4,0.25)
SSE.CV = rep(0,length(loglam))
betafdPari = betafdPar1
for(i in 1:length(loglam)){
    print(loglam[i])
    betafdPari$lambda = 10^loglam[i]
    betalisti = betalist
    for (j in 1:2) betalisti[[j]] = betafdPari
    CVi = fRegress.CV(birdfd3, xfdlist, betalisti,CVobs=1:26)
    SSE.CV[i] = CVi$SSE.CV
}

#  Figure 10.4

plot(loglam,SSE.CV,type='b',cex.lab=1.5,cex.axis=1.5,lwd=2,
  xlab='log smoothing parameter',ylab='cross validated sum of squares')

#  Cross-validation is minimized at something like lambda = sqrt(10),
#  although the discontinous nature of the CV function is disquieting.

betafdPar1$lambda = 10^0.5
for (j in 1:2) betalist[[j]] = betafdPar1

#  carry out the functional regression analysis

#fRegressList = fRegress(birdfd3, xfdlist, betalist)
fitShellfish.5 = fRegress(birdfd3, xfdlist, betalist)

#  plot regression functions

betanames = list("Intercept", "Food Effect")

birdBetaestlist = fitShellfish.5$betaestlist
#betaestlist = fRegressList$betaestlist

op = par(mfrow=c(2,1), cex=1.2)
for (j in 1:2) {
    betaestParfdj = birdBetaestlist[[j]]
    betaestfdj    = betaestParfdj$fd
    betaestvecj   = eval.fd(yearCode, betaestfdj)
	  plot(yearObs, betaestvecj, type="l", lwd=4, col=4,
           xlab="Year", ylab="Temp.",
           main=betanames[[j]])
}
par(op)

#  plot predicted functions

#yhatfdobj = fRegressList$yhatfdobj
birdYhatfdobj = fitShellfish.5$yhatfdobj

plotfit.fd(logCounts2, yearCode, birdYhatfdobj$fd[1:26])
# *** Click on the plot to advance to the next ...

##
## Section 10.2 Functional Responses with Functional Predictors:
##              The Concurrent Model
##

#  Section 10.2.2 Confidence Intervals for Regression Functions

#yhatmat = eval.fd(yearCode, fitShellfish.5$yhatfdobj$fd)
birdYhatmat = eval.fd(yearCode, birdYhatfdobj$fd[1:26])
rmatb   = logCounts2 - birdYhatmat
SigmaEb = var(t(rmatb))

#y2cMap = birdSmoothPar$y2cMap
y2cMap.bird = birdlist2$y2cMap

birdStderrList = fRegress.stderr(fitShellfish.5, y2cMap.bird,
                             SigmaEb)
birdBeta.sdList = birdStderrList$betastderrlist

op = par(mfrow=c(2,1))
plotbeta(birdBetaestlist[1:2], birdBeta.sdList[1:2])
par(op)

####### KNEE ~ HIP #######

# Section 10.2.3 Knee Angle Predicted from Hip Angle

gaittime = seq(0.5,19.5,1)
gaitrange = c(0,20)
gaitfine = seq(0,20,len=101)

harmaccelLfd20 = vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
gaitbasis = create.fourier.basis(gaitrange, nbasis=21)

gaitLoglam = seq(-4,0,0.25)
nglam   = length(gaitLoglam)

# First select smoothing for the raw data

gaitSmoothStats = array(NA, dim=c(nglam, 3),
      dimnames=list(gaitLoglam, c("log10.lambda", "df", "gcv") ) )
gaitSmoothStats[, 1] = gaitLoglam

#  loop through smoothing parameters

for (ilam in 1:nglam) {
  gaitSmooth = smooth.basisPar(gaittime, gait, gaitbasis,
                   Lfdobj=harmaccelLfd20, lambda=10^gaitLoglam[ilam])
  gaitSmoothStats[ilam, "df"]  = gaitSmooth$df
  gaitSmoothStats[ilam, "gcv"] = sum(gaitSmooth$gcv)
  # note: gcv is a matrix in this case
}

#  display and plot GCV criterion and degrees of freedom

gaitSmoothStats
plot(gaitSmoothStats[, c(1, 3)], type='b')

#  set up plotting arrangements for one and two panel displays
#  allowing for larger fonts

op = par(mfrow=c(2,1))
plot(gaitSmoothStats[, c(1, 3)], type="b", log="y")
plot(gaitSmoothStats[, 1:2], type="b", log="y")
par(op)

#    GCV is minimized with lambda = 10^(-1.5).

gaitSmooth = smooth.basisPar(gaittime, gait,
       gaitbasis, Lfdobj=harmaccelLfd20, lambda=10^(-1.5))
gaitfd = gaitSmooth$fd

names(gaitfd$fdnames) = c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] = c("Hip", "Knee")

hipfd  = gaitfd[,1]
kneefd = gaitfd[,2]

# Figure 10.5

kneefdMean = mean(kneefd)

op = par(mfrow=c(3,1))
plot(kneefdMean, xlab='', ylab='', ylim=c(0, 80),
     main='Mean Knee Angle', lwd=2)
abline(v=c(7.5, 14.7), lty='dashed')
plot(deriv(kneefdMean), xlab='', ylab='',
     main='Knee Angle Velocity', lwd=2)
abline(v=c(7.5, 14.7), h=0, lty='dashed')
plot(deriv(kneefdMean, 2), xlab='', ylab='',
     main='Knee Angle Acceleration', lwd=2)
abline(v=c(7.5, 14.7), h=0, lty='dashed')
par(op)

# Figure 10.6

phaseplanePlot(gaitfine, kneefdMean,
               labels=list(evalarg=gaittime, labels=1:20),
               xlab='Knee Velocity', ylab='Knee Acceleration')

# Set up a  functional linear regression

xfdlist   = list(const=rep(1,39), hip=hipfd)
betafdPar = fdPar(gaitbasis, harmaccelLfd20)
betalist  = list(const=betafdPar, hip=betafdPar)

gaitRegress= fRegress(kneefd, xfdlist, betalist)

# Figure 10.7
op = par(mfrow=c(2,1))

# Intercept
betaestlist = gaitRegress$betaestlist
kneeIntercept = predict(betaestlist$const$fd, gaitfine)

# mean knee angle
kneeMean = predict(kneefdMean, gaitfine)

# Plot intercept & mean knee angle
ylim1 = range(kneeIntercept, kneeMean)
plot(gaitfine, kneeIntercept, ylim=ylim1, lwd=2,
     main="Intercept and Mean Knee Angle", type='l',
     xlab='', ylab='')
lines(gaitfine, kneeMean, lty='dashed')
abline(h=0, v=c(7.5, 14.7), lty='dashed')

# Hip coefficient
hipCoef = predict(betaestlist$hip$fd, gaitfine)

# Squared multiple correlation
kneehatfd = gaitRegress$yhatfd$fd
#kneemat = predict(kneefd, gaittime)
kneehatmat = eval.fd(gaittime, kneehatfd)
#resmat. = kneemat - kneehatmat
resmat. = gait[,,'Knee Angle'] - kneehatmat
SigmaE = cov(t(resmat.))

kneefinemat   = eval.fd(gaitfine, kneefd)
kneemeanvec   = eval.fd(gaitfine, mean(kneefd))
kneehatfinemat= eval.fd(gaitfine, kneehatfd)
resmat        = kneefinemat - kneehatfinemat
ncurve        = dim(gait)[2]
resmat0 = kneefinemat - kneemeanvec %*% matrix(1,1,ncurve)
SSE0 = apply((resmat0)^2, 1, sum)
SSE1 = apply(resmat^2, 1, sum)
knee.R2 = (SSE0-SSE1)/SSE0

# Plot Hip Coefficient & Squared Multiple Correlation

ylim2=c(0, max(hipCoef, knee.R2))
plot(gaitfine, hipCoef, lwd=2, xlab='', ylab='', ylim=ylim2, type='l',
     main='Hip Coefficient and Squared Multiple Correlation')
abline(v=c(7.5, 14.7), lty='dashed')
lines(gaitfine, knee.R2, lty='dashed')

# done
par(op)

# Figure 10.8

gaitbasismat = eval.basis(gaitfine, gaitbasis)
#y2cMap0 = solve(crossprod(gaitbasismat), t(gaitbasismat))
y2cMap = gaitSmooth$y2cMap

fRegressList1 = fRegress(kneefd, xfdlist, betalist,
                         y2cMap=y2cMap, SigmaE=SigmaE)

fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
betastderrlist = fRegressList2$betastderrlist
#titlelist = list("Intercept", "Hip coefficient")

op = par(mfrow=c(2,1))
plotbeta(betaestlist, betastderrlist, gaitfine)
par(op)

# Figure 10.9
# fRegress(deriv(kneefd, 2) ~ deriv(hipfd, 2))

xfdlist2 = list(const=rep(1,39), hip=deriv(hipfd, 2))
kneefd.accel = deriv(kneefd, 2)
gaitAccelRegr = fRegress(kneefd.accel, xfdlist2, betalist)

gaitt3 = seq(0, 20, length=401)
beta.hipFine = predict(gaitAccelRegr$betaestlist$hip$fd, gaitt3)

plot(gaitt3, beta.hipFine, type ='l', ylim=c(0, max(beta.hipFine)),
     xlab='', ylab='Hip acceleration and squared multiple correlation',
     lwd=2)
abline(v=c(7.5, 14.7), lty='dashed')

# Squared multiple correlation
# kneeAccel.R2 = var(gaitAccel.Regr$yhatfd) / var(kneefd.accel)

kneeAccel.pred = predict(gaitAccelRegr$yhatfd$fd, gaitt3)
kneeAccel.     = predict(kneefd.accel, gaitt3)

MS.pred = sd(t(kneeAccel.pred))^2
MS.accelfd = sd(t(kneeAccel.))^2
kneeAccel.R2 = (MS.pred / MS.accelfd)

lines(gaitt3, kneeAccel.R2, lty='dashed', lwd=2)

##
## Section 10.3 Beyond the Concurrent Model
##
#  (no computations in this section)

##
## Section 10.4 A Functional Linear Model for Swedish Mortality
##

# From Giles' Sweden.Rdata
Swede.Rdata = 'C:/Users/jim/FDA2008/Sweden.Rdata'
(mat0 = load(Swede.Rdata))
#Swede.Rdata = 'C:/Users/spencerg/fda/Rbook/Rbook/RCode/Sweden.Rdata'
#(mat0 = load(Swede.Rdata))
# SwedeMat Swede1920

# Download the contents of
# http://www.mortality.org/hmd/SWE/STATS/fltcoh_1x1.txt
# to a file 'Sweden Female Lifetable.txt

SwedeTab = read.table('Sweden Female Lifetable.txt')
SwedeTab = read.table('SwedeHazard.dat')
HazMat   = matrix(SwedeTab[,3],111,164,byrow=FALSE)
HazMat   = matrix(as.numeric(HazMat[1:81,]),81,164)
SwedeMat = log(HazMat)


# *****
#
# GILES:  Are these names correct?
#

dimnames(SwedeMat)[[2]] <- paste('b', 1751:1914, sep='')

# *****
#
# GILES:  Are these names correct?
#

# Figure 10.10

matplot(0:80, SwedeLogHazard[, c('b1780', 'b1820', 'b1880', 'b1885')],
        type='b')

SwedeLogHazard$b1900
# Huge spike in 1924 for the 1900 cohort

# Giles said the plot was mislabeled;  try 1760, 1820, 1880, 1920:
Swede4Lines = cbind(SwedeLogHazard[, c('b1760', 'b1820', 'b1880')],
    b1920=Swede1920)

matplot(0:80, Swede4Lines, type='l',lwd=2,xlab='age',ylab='log Hazard',col=1,
matplot(0:80, SwedeMat[, c('b1751', 'b1810', 'b1860', 'b1914')],
        type='b')

SwedeLogHazard$b1900
# Huge spike in 1924 for the 1900 cohort
# in this plot but not in Fig 10.10

matplot(0:80, SwedeMat[, c('b1751', 'b1810', 'b1860', 'b1914')], type='l',lwd=2,xlab='age',ylab='log Hazard',col=1,
  cex.lab=1.5,cex.axis=1.5)
legend(x='bottomright',legend=c('1751','1810', '1860','1914'),lwd=2,col=1,lty=1:4)

# Set up for 'linmod'


SwedeLogHazard <- SwedeMat[,1:144]

SwedeBasis = create.bspline.basis(c(0,80),23)

SwedeBeta0Par = fdPar(SwedeBasis, 2, 1e-5)
SwedeBeta1Par = fdPar(SwedeBasis, 2, 1e3)

SwedeBetaList = list(SwedeBeta0Par, SwedeBeta1Par, SwedeBeta1Par)

D2fdPar = fdPar(SwedeBasis, lambda=1e-7)

SwedeLogHazfd = smooth.basis(0:80, as.matrix(SwedeLogHazard), D2fdPar)$fd

# The following requires manually clicking on the plot
# for each of 144 birth year cohorts
#plotfit.fd(as.matrix(SwedeLogHazard),0:80,SwedeLogHazfd$fd)


NextYear = SwedeLogHazfd[2:144]
LastYear = SwedeLogHazfd[1:143]

Swede.linmodSmooth = linmod(NextYear, LastYear, SwedeBetaList)

# ***
#
# 'linmod' returns NULL
#
# ***????????

Swede.ages = seq(0, 80, 2)
Swede.beta1fd = eval.bifd(Swede.ages, Swede.ages,
    Swede.linmodSmooth$regfd)


# Figure 10.11

# ***?????

##
## Section 10.5 Permutation Tests of Functional Hypotheses
##
#  Section 10.5.1 Functional t-Tests

# Figure 10.12

ylim = with(growth, range(hgtm, hgtf))

with(growth, matplot(age, hgtm[, 1:10], type='l',
                     lty='dashed', ylab='height (cm)'))
with(growth, matlines(age, hgtf[, 1:10], lty='solid'))
legend('topleft', legend=c('girls', 'boys'),
       lty=c('solid', 'dashed'))

growthbasis = create.bspline.basis(breaks=growth$age, norder=6)
growfdPar = fdPar(growthbasis, 3, 10^(-0.5))

hgtffd = with(growth, smooth.basis(age,hgtf,growfdPar))
hgtmfd = with(growth, smooth.basis(age,hgtm,growfdPar))
tres = tperm.fd(hgtffd$fd,hgtmfd$fd)

# Figure 10.13

# Section 10.5.2 Functional F-Tests

# temp36fd, regionList, betaList from Section 10.1.1 above
F.res = Fperm.fd(temp36fd, regionList, betaList)

# Figure 10.14
# plot in black and white

with(F.res,{
            q = 0.95
           ylims = c(min(c(Fvals, qval, qvals.pts)), max(c(Fobs,
                qval)))
            plot(argvals, Fvals, type = "l", ylim = ylims, col = 1,
                lwd = 2, xlab = "day", ylab = "F-statistic", cex.lab=1.5,cex.axis=1.5)
            lines(argvals, qvals.pts, lty = 3, col = 1, lwd = 2)
            abline(h = qval, lty = 2, col = 1, lwd = 2)
            legendstr = c("Observed Statistic", paste("pointwise",
                1 - q, "critical value"), paste("maximum", 1 -
                q, "critical value"))
            legend(argvals[1], 1.2, legend = legendstr,
                col = c(1, 1, 1), lty = c(1, 3, 2), lwd = c(2,
                  2, 2))
        }
)

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
