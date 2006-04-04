
#  -----------------------------------------------------------------------
#                     Daily Weather Data
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  These analyses of the daily temperature and precipitation data for 35
#  Canadian weather stations are designed to be a tour of much of the
#  R and S-PLUS software.  The emphasis is on illustrating useful
#  techniques in a way that will help set up your own application.
#  Sophistication, computational efficiency and compression of code have
#  been ignored as considerations in favor of making the analyses as
#  transparent as possible.
#  The sections in this file are as follows:

#  As always, the first step is to inspect the data and to smooth the
#  data with the appropriate level of smoothing for the analyses being
#  considered.

#  The idea of a customized or "smart" roughness penalty
#  is introduced right away in the form of harmonic acceleration, which
#  is especially important for periodic data such as these with a
#  variation that is dominated by a sinusoid plus a constant signal, or
#  shifted harmonic variation.  However, in order to keep the analyses
#  simple and to economize on computational effort, we often compromise
#  the principle supported in the FDA book of using a saturated basis
#  capable of interpolating the data combined with a roughness penalty,
#  and instead opt for a Fourier series basis system with 65 basis
#  functions and no roughness penalty.

#  Nevertheless, there is a smoothing section below that uses 365
#  Fourier basis functions combined with a harmonic acceleration penalty
#  where we estimate the smoothing parameter by minimizing the
#  GCV or generalized cross-validation parameter.

#  Smoothing is followed by the display of various descriptive statistics,
#  including mean and standard deviation functions,  and covariance and
#  correlation surfaces.

#  A short section illustrates the principal components analysis of
#  temperature, and you may want to repeate these analyses for precipitation.

#  The functional linear model is used in various ways in the following
#  sections:
#  1.  Functional analysis of variance is used to show climate zone
#      effects for temperature and precipitation.
#  2.  The concurrent functional linear model is illustrated by
#      predicting log precipitation from climate zone and a functional
#      covariate constructed by removing climate effects from temperature.
#  3.  Log annual precipitation, a scalar dependent variable, is fit
#      by using temperature as a functional covariate.  Harmonic
#      acceleration roughness in the regression coefficient function
#      is penalized.
#  4.  The full log precipitation function is fit by the regressing
#      on the full temperature profile, and various levels of smoothing
#      are used to show the the effects of smoothing over both arguments
#      of the regression coefficient function.

#  The final section illustrates the smoothing of precipitation by a
#  curve that is constrained to be strictly positive, as is required
#  for a variable like precipitation.
#  -----------------------------------------------------------------------

#  Last modified 20 March 2006

#  ------------------------  input the data  -----------------------


#  set up the times of observation at noon

daytime   <- (1:365)-0.5
dayrange  <- c(0,365)
dayperiod <- 365

#  day values roughly in weeks

weeks <- seq(0,365,length=53)


#  define 1-character names for months
monthletter <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")


#  -------------  set up fourier basis  ---------------------------
#  Here it was decided that 65 basis functions captured enough of
#  the detail in the temperature data: about one basis function
#  per week.  However, see below for smoothing with a saturated
#  basis (365 basis functions) where smoothing is defined by the
#  GCV criterion.

nbasis     <- 65
daybasis65 <- create.fourier.basis(dayrange, nbasis, dayperiod)

#  -----------  set up the harmonic acceleration operator  ----------

harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), dayrange)

#  ---------  create fd objects for temp. and prec. ---------------

daytempfd <- data2fd(daily$tempav, daytime, daybasis65,
                     argnames=list("Day", "Station", "Deg C"))
dayprecfd <- data2fd(daily$precav, daytime, daybasis65,
                     argnames=list("Day", "Station", "Deg C"))

#  set up plotting arrangements for one and two panel displays allowing
#  for larger fonts

#  Plot temperature curves and values

par(mfrow=c(1,1), pty="m")
plotfit.fd(daily$tempav, daytime, daytempfd, titles=daily$place)

#  Plot residuals for three best fits and three worst fits

casenames <- daily$place
varnames  <- "Temperature"
rng       <- dayrange
index     <- c(1,2,3,33,34,35)
residual  <- TRUE
sortwrd   <- TRUE

plotfit.fd(daily$tempav, daytime, daytempfd, rng, index, 366,
           residual, sortwrd, titles=daily$place)

#  Plot precipitation curves and values

plotfit.fd(daily$precav, daytime, dayprecfd, titles=daily$place)

#  Assessment: the functions are definitely too rough with
#  this many basis functions, and especially for precip. which
#  has a much higher noise level.

#  These smoothing parameter values probably undersmooth the data,
#  but we can impose further smoothness on the results of analyses

#  Set up the functional parameter objects to define smoothing

templambda <- 1e1
preclambda <- 1e5
tempfdPar  <- fdPar(daybasis65, harmaccelLfd, templambda)
precfdPar  <- fdPar(daybasis65, harmaccelLfd, preclambda)

daytempfd <- smooth.fd(daytempfd, tempfdPar)
dayprecfd <- smooth.fd(dayprecfd, precfdPar)

#  Use function PLOTFIT.FD to plot the temperature data, fit and residuals
#    sorted by size of root mean square error

#  Plot temperature curves and values

plotfit.fd(daily$tempav, daytime, daytempfd, titles=daily$place)

#  Plot residuals for three best fits and three worst fits

plotfit.fd(daily$tempav, daytime, daytempfd, index=c(1,2,3,33,34,35),
          sortwrd=TRUE, residual=TRUE, titles=daily$place)

#  Plot curves and values only for January

plotfit.fd(daily$tempav, daytime, daytempfd, rng=c(0,31), titles=daily$place)

#  plot each pair of functions along with raw data

par(mfrow=c(1,2), pty="s")
index <- 1:35
for (i in index) {
    #par(ask=TRUE)
    plot(daytime,daily$tempav[,i], type="p", xlim=dayrange, col=1,
         xlab="Day", ylab="", main=paste(daily$place[i],"temperature"))
    lines.fd(daytempfd[i])
    plot(daytime,daily$precav[,i], type="p", xlim=dayrange, col=1,
         xlab="Day", ylab="", main=paste(daily$place[i],"precipitation"))
    lines.fd(dayprecfd[i])
}

#  plot all the functions

par(mfrow=c(1,2), pty="s")
plot(daytempfd, main="Temperature")
plot(dayprecfd, main="Precipitation")

#  -------------------------------------------------------------
#                 Choose level of smoothing using
#          the generalized cross-validation criterion
#              with smoothing function smooth.basis.
#  -------------------------------------------------------------

wtvec <- rep(1,365)

# set up a saturated basis capable of interpolating the data

nbasis      <- 365
daybasis365 <- create.fourier.basis(dayrange, nbasis)

#  --------------------  smooth temperature  ------------------

#  set up range of smoothing parameters in log_10 units

loglam <- -5:1
nlam   <- length(loglam)

dfsave  <- rep(0,nlam)
gcvsave <- rep(0,nlam)

#  loop through smoothing parameters

for (ilam in 1:nlam) {
	lambda <- 10^loglam[ilam]
	fdParobj <- fdPar(daybasis365, harmaccelLfd, lambda)
	smoothlist <- smooth.basis(daytime, daily$tempav, fdParobj)
	fdobj  <- smoothlist[[1]]
	df     <- smoothlist[[2]]
	gcv    <- smoothlist[[3]]
	dfsave[ilam]  <- df
	gcvsave[ilam] <- sum(gcv)
}

#  display and plot degrees of freedom and GCV criterion

cbind(loglam, dfsave, gcvsave)


cbind(loglam, dfsave, gcvsave)

par(mfrow=c(1,2), pty="s")
plot(loglam, gcvsave, type="b", cex=1,
     xlab="Log_10 lambda", ylab="GCV Criterion",
     main="Temperature Smoothing")

plot(loglam, dfsave, type="b",  cex=1,
     xlab="Log_10 lambda", ylab="Degrees of freedom",
     main="Temperature Smoothing")

#  Do final smooth with minimum GCV value

lambda   <- 0.01  #  minimum GCV estimate, corresponding to 255 df
fdParobj <- fdPar(daybasis365, harmaccelLfd, lambda)

smoothlist <- smooth.basis(daytime, daily$tempav, fdParobj)
daytempfd  <- smoothlist$fd
df         <- smoothlist$df
gcv        <- smoothlist$gcv
coef       <- smoothlist$coef
SSE        <- smoothlist$SSE

#  estimate standard error of fit

stderr <- sqrt(SSE/(35*(365-df)))  #  0.29 deg C

#  plot data and fit

par(mfrow=c(1,1), pty="m")
plotfit.fd(daily$tempav, daytime, daytempfd, titles=daily$place)

#  --------------------  smooth precipitation  ------------------

#  set up range of smoothing parameters in log.10 units

loglam <- 4:9
nlam <- length(loglam)

dfsave  <- rep(0,nlam)
gcvsave <- rep(0,nlam)

#  loop through smoothing parameters

for (ilam in 1:nlam) {
	lambda <- 10^loglam[ilam]
	cat(paste("lambda =",lambda,"\n"))
	fdParobj <- fdPar(daybasis365, harmaccelLfd, lambda)
	smoothlist <- smooth.basis(daytime, daily$precav, fdParobj)
	fdobj  <- smoothlist[[1]]
	df     <- smoothlist[[2]]
	gcv    <- smoothlist[[3]]
	dfsave[ilam]  <- df
	gcvsave[ilam] <- sum(gcv)
}

#  display and plot degrees of freedom and GCV criterion

cbind(loglam, dfsave, gcvsave)

par(mfrow=c(1,2), pty="m")
plot(loglam, gcvsave, type="b", cex=1,
     xlab="Log_10 lambda", ylab="GCV Criterion",
     main="Precipitation Smoothing")

plot(loglam, dfsave, type="b", cex=1,
     xlab="Log_10 lambda", ylab="Degrees of freedom",
     main="Precipitation Smoothing")

#  Do final smooth with minimum GCV value

lambda   <- 1e7  #  minimum GCV estimate, corresponding to 255 df
fdParobj <- fdPar(daybasis365, harmaccelLfd, lambda)

smoothlist <- smooth.basis(daytime, daily$precav, fdParobj)

dayprecfd  <- smoothlist$fd
df         <- smoothlist$df
gcv        <- smoothlist$gcv
coef       <- smoothlist$coef
SSE        <- smoothlist$SSE

#  estimate standard error of fit

stderr <- sqrt(SSE/(35*(365-df)))  #  0.94 mm

#  plot data and fit

par(mfrow=c(1,1), pty="m")
plotfit.fd(daily$precav, daytime, dayprecfd, titles=daily$place)

#  Assessment: the temperature curves are still pretty rough,
#  although the data themselves show that there are very
#  high frequency effects in the mean temperature, especially
#  early in the year.
#  The precip. curves may be oversmoothed for (some weather
#  stations.

#  smooth precipitation in Prince Rupert

PRprecfd <- smooth.basis(daytime, daily$precav[,29], fdParobj)$fd

PRprecvec <- eval.fd(daytime, PRprecfd)

par(ask=FALSE)
plot(daytime, daily$precav[,29], type="p", ylim=c(0,8),
     xlab="Day", ylab="Precipitation (mm)",
     main="Prince Rupert")
lines(daytime, PRprecvec, lwd=2)

#  ----------------------------------------------------------------------
#                Descriptive Statistics Functions
#  ----------------------------------------------------------------------

#  ---------  create fd objects for temp. and prec. ---------------

daytempfd <- data2fd(daily$tempav, daytime, daybasis65,
                     argnames=list("Day", "Station", "Deg C"))
dayprecfd <- data2fd(daily$precav, daytime, daybasis65,
                     argnames=list("Day", "Station", "Deg C"))

#  --  compute and plot mean and standard deviation of temperature -------

tempmeanfd  <- mean.fd(daytempfd)
tempstdvfd  <- stddev.fd(daytempfd)

par(mfrow=c(1,2), pty="s")
plot(tempmeanfd,               main="Mean")
plot(tempstdvfd, ylim=c(0,10), main="Standard Deviation")

#  --  plot the temperature variance-covariance bivariate function  ----

tempvarbifd <- var.fd(daytempfd)
tempvarmat  <- eval.bifd(weeks,weeks,tempvarbifd)

par(mfrow=c(1,2), pty="s")
contour(tempvarmat, xlab="Days", ylab="Days")
persp(tempvarmat,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Temperature Covariance", line=-4, outer=TRUE)

#  --  plot the temperature correlation function  ---

tempstddev <- sqrt(diag(tempvarmat))
tempcormat <- tempvarmat/outer(tempstddev,tempstddev)
tempcormat <- tempvarmat/(tempstddev %o% tempstddev)

contour(tempcormat, xlab="Days", ylab="Days")
persp(tempcormat,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Temperature Correlation", line=-4, outer=TRUE)

#  --  plot the precipitation variance-covariance bivariate function  ----

weeks <- seq(0,365,length=53)
precvarbifd <- var.fd(dayprecfd)
precvarmat  <- eval.bifd(weeks,weeks,precvarbifd)

contour(precvarmat, xlab="Days", ylab="Days")
persp(precvarmat,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Precipitation Covariance", line=-4, outer=TRUE)

#  --  plot the precipitation correlation function  ---

precstddev <- sqrt(diag(precvarmat))
preccormat <- precvarmat/(precstddev %o% precstddev)

contour(preccormat, xlab="Days", ylab="Days")
persp(preccormat,
      xlab="Days", ylab="Days", zlab="Covariance")
mtext("Precipitation Correlation", line=-4, outer=TRUE)

#  -----  compute and plot the covariance between temp. and prec.  --

covbifd <- var.fd(daytempfd, dayprecfd)
covmat  <- eval.bifd(weeks,weeks,covbifd)

contour(covmat, xlab="Weeks", ylab="Weeks")
persp(  covmat, xlab="Weeks", ylab="Weeks", zlab="Covariance")
mtext("Temperature-Precipitation Covariance", line=-4, outer=TRUE)

#  -----  compute and plot the correlation between temp. and prec.  --

cormat  <- covmat/(tempstddev %o% precstddev)

contour(cormat, xlab="Weeks", ylab="Weeks")
persp(cormat, xlab="Weeks", ylab="Weeks", zlab="Correlation")
mtext("Temperature-Precipitation Correlation", line=-4, outer=TRUE)

#  -----------------------------------------------------------------------
#               PCA of temperatures with varimax rotation
#  -----------------------------------------------------------------------

harmfdPar     <- fdPar(daybasis65, harmaccelLfd, 1e5)
daytemppcaobj <- pca.fd(daytempfd, nharm=4, harmfdPar)

daytemppcaobj <- varmx.pca.fd(daytemppcaobj)

#  plot harmonics

par(mfrow=c(1,1), pty="m")
plot.pca.fd(daytemppcaobj)

#  plot log eigenvalues

daytempeigvals <- daytemppcaobj[[2]]
par(ask=FALSE)
plot(1:20, log10(daytempeigvals[1:20]), type="b",
     xlab="Eigenvalue Number", ylab="Log 10 Eigenvalue")
abline(lsfit(5:20, log10(daytempeigvals[5:20])), lty=2)

#  plot factor scores

harmscr <- daytemppcaobj[[3]]

plot(harmscr[,1], harmscr[,2], xlab="Harmonic 1", ylab="Harmonic 2")
text(harmscr[,1], harmscr[,2], daily$place, col=4)

plot(harmscr[,3], harmscr[,4], xlab="Harmonic 3", ylab="Harmonic 4")
text(harmscr[,3], harmscr[,4], daily$place, col=4)

#  ------------------------------------------------------------------
#               Functional linear models
#  ------------------------------------------------------------------

#  ---------------------------------------------------------------
#             Predicting temperature from climate zone
#  ---------------------------------------------------------------

#  return data to original ordering

daily$tempav <- matrix(scan("../data/dailtemp.txt",0), 365, 35)
daily$precav <- matrix(scan("../data/dailprec.txt",0), 365, 35)

#  set up a smaller basis using only 65 Fourier basis functions
#  to save some computation time

smallnbasis <- 65
smallbasis  <- create.fourier.basis(dayrange, smallnbasis)

tempfd      <- data2fd(daily$tempav, daytime, smallbasis)

smallbasismat <- eval.basis(daytime, smallbasis)
y2cMap <- solve(crossprod(smallbasismat)) %*% t(smallbasismat)

#  names for (climate zones

zonenames <- c("Canada  ",
               "Atlantic", "Pacific ", "Contintl", "Arctic  ")

#  indices for (weather stations in each of four climate zones

atlindex <- c(1,2,4,8,9,13,14,15,19,22,23,24,25,28,34)
pacindex <- c(12,17,18,30,31)
conindex <- c(3,5,6,7,16,20,26,27,29,32,33,35)
artindex <- c(10,11,21)

#  Set up a design matrix having a column for (the grand mean, and
#    a column for (each climate zone effect. Add a dummy contraint
#    observation

zmat <- matrix(0,35,5)
zmat[        ,1] <- 1
zmat[atlindex,2] <- 1
zmat[pacindex,3] <- 1
zmat[conindex,4] <- 1
zmat[artindex,5] <- 1

#  attach a row of 0, 1, 1, 1, 1 to force zone
#  effects to sum to zero, and define first regression
#  function as grand mean for (all stations

z36    <- matrix(1,1,5)
z36[1] <- 0
zmat   <- rbind(zmat, z36)

#  revise YFDOBJ by adding a zero function

coef   <- tempfd$coefs
coef36 <- cbind(coef,matrix(0,smallnbasis,1))
tempfd$coefs <- coef36

p <- 5
xfdlist <- vector("list",p)
for (j in 1:p) xfdlist[[j]] <- zmat[,j]

#  set up the basis for (the regression functions

nbetabasis <- 11
betabasis  <- create.fourier.basis(dayrange, nbetabasis)

#  set up the functional parameter object for (the regression fns.

betafd    <- fd(matrix(0,nbetabasis,1), betabasis)
estimate  <- TRUE
lambda    <- 0
betafdPar <- fdPar(betafd, harmaccelLfd, lambda, estimate)

betalist <- vector("list",p)
for (j in 1:p) betalist[[j]] <- betafdPar

#  compute regression coefficient functions and
#  predicted functions

fRegressList <- fRegress(tempfd, xfdlist, betalist)

#  plot regression functions

betaestlist <- fRegressList$betaestlist
par(mfrow=c(3,2))
for (j in 1:p) {
	betaestParfdj <- betaestlist[[j]]
	plot(betaestParfdj$fd, xlab="Day", ylab="Temp.",
	     main=zonenames[j])
}

#  plot predicted functions

yhatfdobj <- fRegressList$yhatfdobj
plot(yhatfdobj)

#  compute residual matrix and get covariance of residuals

yhatmat  <- eval.fd(daytime, yhatfdobj)
ymat     <- eval.fd(daytime, tempfd)
temprmat <- ymat[,1:35] - yhatmat[,1:35]
SigmaE   <- var(t(temprmat))

#  plot covariance surface for errors

par(mfrow=c(1,1))
contour(SigmaE, xlab="Day", ylab="Day")
lines(dayrange,dayrange,lty=4)

#  plot standard deviation of errors

par(mfrow=c(1,1), mar=c(5,5,3,2), pty="m")
stddevE <- sqrt(diag(SigmaE))
plot(daytime, stddevE, type="l",
     xlab="Day", ylab="Standard error (deg C)")

#  Repeat regression, this time outputting results for
#  confidence intervals

stderrList <- fRegress.stderr(fRegressList, y2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  plot regression function standard errors

par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
	betastderrj <- eval.fd(daytime, betastderrlist[[j]])
	plot(daytime, betastderrj,
	        type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
	        main=zonenames[j])
}

#  plot regression functions with confidence limits

par(mfrow=c(2,3))
for (j in 1:p) {
	betafdParj  <- betaestlist[[j]]
	betafdj     <- betafdParj$fd
	betaj       <- eval.fd(daytime, betafdj)
	betastderrj <- eval.fd(daytime, betastderrlist[[j]])
	matplot(daytime, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
	        type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
	        main=zonenames[j])
}

#  -----------------------------------------------------------------------
#         predict log precipitation from climate zone and temperature
#  -----------------------------------------------------------------------

#  Be sure to run previous analysis predicting temperature from
#  climate zone before running this example.

#  set up functional data object for log precipitation

precfd     <- data2fd(daily$precav, daytime, smallbasis)

logprecmat <- log10(eval.fd(daytime, precfd))

lnprecfd <- data2fd(logprecmat, daytime, smallbasis)
lnprecfd$fdnames[[1]] <- "Days"
lnprecfd$fdnames[[2]] <- "Station"
lnprecfd$fdnames[[3]] <- "log.{10} mm"

#  plot log precipitation functions

par(mfrow=c(1,1), pty="m")
plot(lnprecfd)
title("Log Precipitation Functions")

#  revise LOGPREDFD by adding a zero function

coef   <- lnprecfd$coefs
nbasis <- smallbasis$nbasis
coef36 <- cbind(coef,matrix(0,nbasis,1))
lnprecfd$coefs <- coef36

#  set up the XFDLIST list

p <- 6
xfdlist <- vector("list",p)

#  load first five members with columns of design matrix

for (j in 1:5) xfdlist[[j]] <- zmat[,j]

#  set up a FD object for (temperature residuals

lambda     <- 1e5
fdParobj   <- fdPar(smallbasis, harmaccelLfd, lambda)
smoothList <- smooth.basis(daytime, temprmat, fdParobj)
temprfdobj <- smoothList$fd

#  plot temperature residuals

par(mfrow=c(1,1), pty="m")
plot(temprfdobj)

#  extend temperature residual functions to include
#  zero function

coef   <- temprfdobj$coefs
nbasis <- dim(coef)[1]
coef36 <- cbind(coef,matrix(0,nbasis,1))
temprfdobj$coefs <- coef36

#  add TEMPRFDOBJ to the set of predictors

xfdlist[[6]]  <- temprfdobj
betalist[[6]] <- betafdPar

#  set up the basis for (the regression functions

nbetabasis <- 13
betabasis  <- create.fourier.basis(dayrange, nbetabasis)

#  set up the functional parameter object for (the regression fns.

betafd    <- fd(matrix(0,nbetabasis,p), betabasis)
estimate  <- TRUE
lambda    <- 0
betafdPar <- fdPar(betafd, harmaccelLfd, lambda, estimate)
for (j in 1:p) betalist[[j]] <- betafdPar

#  compute regression coefficient functions and
#  predicted functions

fRegressList <- fRegress(lnprecfd, xfdlist, betalist)

betaestlist <- fRegressList$betaestlist
yhatfdobj   <- fRegressList$yhatfdobj

#  plot regression functions

prednames <- c(zonenames, "tempres ")
par(mfrow=c(2,3),pty="s")
for (j in 1:p) {
	betaParfdj <- betaestlist[[j]]
	betafdj    <- betaParfdj$fd
    plot(betafdj)
    title(prednames[j])
}

#  plot predicted functions

par(mfrow=c(1,1), pty="m")
plot(yhatfdobj)

#  compute residual matrix and get covariance of residuals

yhatmat    <- eval.fd(daytime, yhatfdobj)
ymat       <- eval.fd(daytime, lnprecfd)
lnprecrmat <- ymat[,1:35] - yhatmat[,1:35]
SigmaE     <- var(t(lnprecrmat))

contour(SigmaE)

#  repeat regression analysis to get confidence intervals

stderrList <- fRegress.stderr(fRegressList, y2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  plot regression functions

prednames <- c(zonenames, "tempres ")

#  plot regression function standard errors

par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
	betastderrfdj <- betastderrlist[[j]]
	betastderrj <- eval.fd(daytime, betastderrfdj)
	plot(daytime, betastderrj,
	        type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
	        main=prednames[j])
}

#  plot regression functions with confidence limits

par(mfrow=c(2,3), pty="s")
for (j in 1:p) {
	betafdParj  <- betaestlist[[j]]
	betafdj     <- betafdParj$fd
	betaj       <- eval.fd(daytime, betafdj)
	betastderrfdj <- betastderrlist[[j]]
	betastderrj   <- eval.fd(daytime, betastderrfdj)
	matplot(daytime, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
	        type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
	        main=prednames[j])
}

#  ---------------------------------------------------------------
#      log annual precipitation predicted by temperature profile
#  ---------------------------------------------------------------

#  set up log10 total precipitation

annualprec <- log10(apply(daily$precav,2,sum))

#  set up a smaller basis using only 65 Fourier basis functions
#  to save some computation time

smallnbasis <- 65
smallbasis  <- create.fourier.basis(dayrange, smallnbasis)
tempfd      <- data2fd(daily$tempav, daytime, smallbasis)

smallbasismat <- eval.basis(daytime, smallbasis)
y2cMap <- solve(crossprod(smallbasismat)) %*% t(smallbasismat)

#  set up the covariates, the first the constant, and the second
#  temperature

p <- 2
constantfd <- fd(matrix(1,1,35), create.constant.basis(dayrange))

xfdlist <- vector("list",2)
xfdlist[[1]] <- constantfd
xfdlist[[2]] <- tempfd[1:35]

#  set up the functional parameter object for (the regression fns.
#  the smoothing parameter for (the temperature function
#  is obviously too small here, and will be revised below by
#  using cross-validation.

betalist   <- vector("list",2)
#  set up the first regression function as a constant
betabasis1 <- create.constant.basis(dayrange)
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)
betalist[[1]] <- betafdPar1
#  set up the second with same basis as for (temperature
#  35 basis functions would permit a perfect fit to the data
nbetabasis  <- 35
betabasis2  <- create.fourier.basis(dayrange, nbetabasis)
betafd2     <- fd(matrix(0,nbetabasis,1), betabasis2)
lambda      <- 10
betafdPar2  <- fdPar(betafd2, harmaccelLfd, lambda)
betalist[[2]] <- betafdPar2

#  carry out the regression analysis

fRegressList <- fRegress(annualprec, xfdlist, betalist)

betaestlist   <- fRegressList$betaestlist

annualprechat <- fRegressList$yhatfdobj

#  constant term

alphafdPar <- betaestlist[[1]]

alphafdPar$fd$coefs

#  plot the coefficient function for (temperature

betafdPar <- betaestlist[[2]]
betafd    <- betafdPar$fd
par(mfrow=c(1,1), pty="m")
plot(betafd)
title("Regression coefficient for temperature")

#  plot the fit

plot (annualprechat, annualprec, type="p")
lines(annualprechat, annualprechat, lty=4)

#  compute cross-validated SSE"s for a range of smoothing parameters

loglam <- seq(5,15,0.5)
nlam   <- length(loglam)
SSE.CV <- matrix(0,nlam,1)
for (ilam in 1:nlam) {
    lambda       <- 10^loglam[ilam]
    betalisti    <- betalist
    betafdPar2   <- betalisti[[2]]
    betafdPar2$lambda <- lambda
    betalisti[[2]] <- betafdPar2
    SSE.CV[ilam]   <- fRegress.CV(annualprec, xfdlist, betalisti)
    print(c(ilam, loglam[ilam], SSE.CV[ilam]))
}

plot(loglam, SSE.CV, type="b",
	  xlab="log_{10} smoothing parameter lambda",
	  ylab="Cross-validation score")

#  analysis with minimum CV smoothing

lambda        <- 10^12.5
betafdPar2    <- fdPar(betafd2, harmaccelLfd, lambda)
betalist[[2]] <- betafdPar2

#  carry out the regression analysis

fRegressList <- fRegress(annualprec, xfdlist, betalist)

betaestlist   <- fRegressList$betaestlist

annualprechat <- fRegressList$yhatfdobj

#  constant term

alphafdPar <- betaestlist[[1]]

alphafdPar$fd$coefs

#  plot the coefficient function for (temperature

betafdPar <- betaestlist[[2]]
betafd    <- betafdPar$fd

par(mfrow=c(1,1), pty="m")
plot(betafd)
title("Regression coefficient for temperature")

#  plot the fit

par(mfrow=c(1,1), pty="m")
plot (annualprechat, annualprec, type="p")
lines(annualprechat, annualprechat, lty=2)

#  compute squared multiple correlation

covmat <- var(cbind(annualprec, annualprechat))
Rsqrd <- covmat[1,2]^2/(covmat[1,1]*covmat[2,2])
#   0.7540

#  compute SigmaE

resid  <- annualprec - annualprechat
SigmaE <- mean(resid^2)
SigmaE <- SigmaE*diag(rep(1,35))

#  recompute the analysis to get confidence limits

stderrList <- fRegress.stderr(fRegressList, NULL, SigmaE)

betastderrlist <- stderrList$betastderrlist

#  constant  coefficient standard error:

betastderr1 <- betastderrlist[[1]]
betastderr1$coefs

#  plot the temperature coefficient function

betafdParj      <- betaestlist[[2]]
betafd          <- betafdParj$fd
betastderrfd    <- betastderrlist[[2]]
betavec         <- eval.fd(daytime, betafd)
betastderrvec   <- eval.fd(daytime, betastderrfd)

betaplotmat <- cbind(betavec, betavec+2*betastderrvec,
                              betavec-2*betastderrvec)

matplot(daytime, betaplotmat, type="l", lty=c(1,4,4),
        xlab="Day", ylab="Temperature Reg. Coeff.")
lines(dayrange,c(0,0),lty=2)

#  ---------------------------------------------------------------
#         predict log precipitation from temperature
#  ---------------------------------------------------------------

#  set up a smaller basis using only 65 Fourier basis functions
#  to save some computation time

smallnbasis <- 65
smallbasis  <- create.fourier.basis(dayrange, smallnbasis)

tempfd      <- data2fd(daily$tempav, daytime, smallbasis)

#  change 0's to 0.05 mm in precipitation data

prectmp <- daily$precav
for (j in 1:35) {
    index <- prectmp[,j] == 0
    prectmp[index,j] <- 0.05
}

#  work with log base 10 precipitation

logprecmat <- log10(prectmp)

#  set up functional data object for log precipitation

lnprecfd <- data2fd(logprecmat, daytime, smallbasis)
lnprecfdnames <- vector("list",3)
lnprecfdnames[[1]] <- "Days"
lnprecfdnames[[2]] <- "Station"
lnprecfdnames[[3]] <- "log.{10} mm"
lnprecfd$fdnames <- lnprecfdnames

#  plot precipitation functions

par(mfrow=c(1,1), pty="m")
plot(lnprecfd)
title("Log Precipitation Functions")

#  set up smoothing levels for (s (xLfd) and for (t (yLfd)

xLfdobj <- harmaccelLfd
yLfdobj <- harmaccelLfd
xlambda <- 1e9
ylambda <- 1e7

#  compute the linear model

wtvec <- matrix(1,35,1)

linmodstr <- linmod(daytempfd, lnprecfd, wtvec,
                    xLfdobj, yLfdobj, xlambda, ylambda)

afd <- linmodstr$alphafd   #  The intercept function
bfd <- linmodstr$regfd     #  The bivariate regression function

#  plot the intercept function

plot(afd, xlab="Day", ylab="Intercept function")

#  plot the regression function as a surface

bfdmat <- eval.bifd(weeks, weeks, bfd)

persp(weeks, weeks, bfdmat, xlab="Day(t)", ylab="Day(s)")

#  Get fitted functions

lnprechatfd <- linmodstr$yhatfd

# Compute mean function as a benchmark for (comparison

lnprecmeanfd <- mean(lnprecfd)
lnprechat0   <- eval.fd(weeks, lnprecmeanfd)

#  Plot actual observed, fitted, and mean log precipitation for
#      each weather station,

lnprecmat    <- eval.fd(weeks, lnprecfd)
lnprechatmat <- eval.fd(weeks, lnprechatfd)
plotrange    <- c(min(lnprecmat),max(lnprecmat))
#par(ask=TRUE)
for (i in 1:35) {
    lnpreci    <- eval.fd(lnprecfd[i],    weeks)
    lnprechati <- eval.fd(lnprechatfd[i], weeks)
    SSE <- sum((lnpreci-lnprechati)^2)
    SSY <- sum((lnpreci-lnprechat0)^2)
    RSQ <- (SSY-SSE)/SSY
    plot(weeks, lnprecmat[,i], type="l", lty=4, ylim=plotrange,
         xlab="Day", ylab="Log Precipitation",
         main=paste(daily$place[i],"  R^2 =",round(RSQ,3)))
    lines(weeks, lnprechatmat[,i])
}
par(ask=FALSE)

#  -------------------------------------------------------------------
#              Smooth Vancouver's precipitation with a
#                       positive function.
#  -------------------------------------------------------------------

#  select Vancouver's precipitation

index <- (1:35)[daily$place == "Vancouver  "] 
VanPrec  <- daily$precav[,index]

#  smooth the data using 65 basis functions

lambda    <- 1e4
fdParobj  <- fdPar(smallbasis, harmaccelLfd, lambda)
VanPrecfd <- smooth.basis(daytime, VanPrec, fdParobj)$fd

#  Plot temperature curves and values

plotfit.fd(VanPrec, daytime, VanPrecfd, titles=daily$place[26])

#  set up the functional parameter object for (positive smoothing

dayfdPar <- fdPar(smallbasis, harmaccelLfd, lambda)

#  smooth the data with a positive function

Wfd1 <- smooth.pos(daytime, VanPrec, dayfdPar)$Wfdobj

#  plot both the original smooth and positive smooth

VanPrecvec    <- eval.fd(daytime, VanPrecfd)
VanPrecposvec <- eval.posfd(daytime, Wfd1)

par(ask=FALSE)
plot(daytime,  VanPrec, type="p")
lines(daytime, VanPrecposvec, lwd=2, col=4)
lines(daytime, VanPrecvec, lwd=2, col=3)
legend(100, 8, c("Positive smooth", "Unrestricted smooth"), 
       lty=c(1,1), col=c(4,3))

#  plot the residuals

VanPrecres <- VanPrec - VanPrecposvec
plot(daytime, VanPrecres^2, type="p")
title("Squared residuals from positive fit")

#  compute a postive smooth of the squared residuals

lambda <- 1e3
dayfdPar <- fdPar(smallbasis, harmaccelLfd, lambda)

Wfd <- smooth.pos(daytime, VanPrecres^2, dayfdPar)$Wfdobj

#  plot the square root of this smooth along with the residuals

VanPrecvarhat <- eval.posfd(daytime, Wfd)
VanPrecstdhat <- sqrt(VanPrecvarhat)

plot(daytime, VanPrecres^2, type="p")
lines(daytime, VanPrecvarhat, lwd=2)

#  set up a weight function for (revised smoothing

wtvec <- 1/VanPrecvarhat

lambda   <- 1e3
dayfdPar <- fdPar(smallbasis, harmaccelLfd, lambda)

Wfd2 <- smooth.pos(daytime, VanPrec, dayfdPar, wtvec)$Wfdobj

#  plot the two smooths, one with weighting, one without

VanPrecposvec2 <- eval.posfd(daytime, Wfd2)

plot(daytime,  VanPrec, type="p")
lines(daytime, VanPrecposvec2, lwd=2, col=4)
lines(daytime, VanPrecposvec, lwd=2, col=3)
legend(100, 8, c("Weighted", "Unweighted"), 
       lty=c(1,1), col=c(4,3))