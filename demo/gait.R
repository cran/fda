#  --------------------------------------------------------------------
#                        Gait data
#  --------------------------------------------------------------------

#  --------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  The gait data were chosen for these sample analyses because they are
#  bivariate:  consisting of both hip and knee angles observed over a
#  gait cycle for 39 children.  The bivariate nature of the data implies
#  certain displays and analyses that are not usually considered, and
#  especially the use of canonical correlation analysis.
#
#  As with the daily weather data, the harmonic acceleration roughness
#  penalty is used throughout since the data are periodic with a strong
#  sinusoidal component of variation.
#
#  After setting up the data, smoothing the data using GCV (generalized
#  cross-validation) to select a smoothing parameter, and displaying
#  various descriptive results, the data are subjected to a principal
#  components analysis, followed by a canonical correlation analysis of
#  thejoint variation of hip and knee angle, and finally a registration
#  of the curves.  The registration is included here especially because
#  the registering of periodic data requires the estimation of a phase
#  shift constant for each curve in addition to possible nonlinear
#  transformations of time.
#
#  --------------------------------------------------------------------

#  Last modified 21 November 2008 by Jim Ramsay
#  Previously modified 25 February 2007 by Spencer Graves

#  attach the FDA functions

library(fda)

#  Set up the argument values: equally spaced over circle of
#  circumference 20.  Earlier  analyses of the gait data used time
#  values over [0,1], but led to singularity problems in the use of
#  function fRegress.  In general, it is better use a time interval
#  that assigns roughly one time unit to each inter-knot interval.

(gaittime <- as.numeric(dimnames(gait)[[1]])*20)
gaitrange <- c(0,20)

#  set up a three-dimensional array of function values

apply(gait, 3, range)

# -----------  set up the harmonic acceleration operator  ----------

harmaccelLfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)

#  Set up basis for representing gait data.  The basis is saturated
#  since there are 20 data points per curve, and this set up defines
#  21 basis functions.  Recall that a fourier basis has an odd number
#  of basis functions.

gaitbasis <- create.fourier.basis(gaitrange, nbasis=21)

#  -------------------------------------------------------------------
#                 Choose level of smoothing using
#          the generalized cross-validation criterion
#  -------------------------------------------------------------------

#  set up range of smoothing parameters in log_10 units

gaitLoglam <- seq(-4,0,0.25)
nglam   <- length(gaitLoglam)

gaitSmoothStats <- array(NA, dim=c(nglam, 3),
      dimnames=list(gaitLoglam, c("log10.lambda", "df", "gcv") ) )
gaitSmoothStats[, 1] <- gaitLoglam

#  loop through smoothing parameters

for (ilam in 1:nglam) {
  gaitSmooth <- smooth.basisPar(gaittime, gait, gaitbasis,
                   Lfdobj=harmaccelLfd, lambda=10^gaitLoglam[ilam])
  gaitSmoothStats[ilam, "df"]  <- gaitSmooth$df
  gaitSmoothStats[ilam, "gcv"] <- sum(gaitSmooth$gcv)
  # note: gcv is a matrix in this case
}

#  display and plot GCV criterion and degrees of freedom

gaitSmoothStats
plot(gaitSmoothStats[, 1], gaitSmoothStats[, 3])

#  set up plotting arrangements for one and two panel displays
#  allowing for larger fonts

op <- par(mfrow=c(2,1))
plot(gaitLoglam, gaitSmoothStats[, "gcv"], type="b",
     xlab="Log_10 lambda", ylab="GCV Criterion",
     main="Gait Smoothing", log="y")

plot(gaitLoglam, gaitSmoothStats[, "df"], type="b",
     xlab="Log_10 lambda", ylab="Degrees of freedom",
     main="Gait Smoothing")
par(op)

# With gaittime <- (1:20)/21,
#    GCV is minimized with lambda = 10^(-2).

str(gait)
gaitfd <- smooth.basisPar(gaittime, gait,
       gaitbasis, Lfdobj=harmaccelLfd, lambda=1e-2)$fd

str(gaitfd)
names(gaitfd$fdnames) <- c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] <- c("Hip", "Knee")

#  --------  plot curves and their first derivatives  ----------------

#par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
op <- par(mfrow=c(2,1))
plot(gaitfd, cex=1.2)
par(op)

#  plot each pair of curves interactively

#par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
op <- par(mfrow=c(2,1))
plotfit.fd(gait, gaittime, gaitfd, cex=1.2)
# Problem:  does not work properly;
# ask=TRUE with appropriate choices for lty, etc.,
# might make it stop after each one,
# but would not fix the overplotting
# Need to modify plotfit.fd to accept ylim = array,
# not just a vector of length 2.
par(op)

#  plot the residuals, sorting cases by residual sum of squares

#par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
plotfit.fd(gait, gaittime, gaitfd, residual=TRUE, sort=TRUE, cex=1.2)

#  plot first derivative of all curves

#par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
plot(gaitfd, Lfdob=1)

#  -----------------------------------------------------------------
#            Display the mean, variance and covariance functions
#  -----------------------------------------------------------------

#  ------------  compute the mean functions  --------------------

gaitmeanfd <- mean.fd(gaitfd)

#  plot these functions and their first two derivatives

#par(mfcol=c(2,3),pty="s")
#op <- par(mfrow=c(3,2))
op <- par(mfcol=2:3)
plot(gaitmeanfd)
plot(gaitmeanfd, Lfdobj=1)
plot(gaitmeanfd, Lfdobj=2)
par(op)

#  --------------  Compute the variance functions  -------------

gaitvarbifd <- var.fd(gaitfd)
str(gaitvarbifd)

gaitvararray <- eval.bifd(gaittime, gaittime, gaitvarbifd)

#par(mfrow=c(1,1), mar=c(3,4,2,1), pty="m")

#  plot variance and covariance functions as contours

filled.contour(gaittime, gaittime, gaitvararray[,,1,1], cex=1.2)
title("Knee - Knee")

filled.contour(gaittime, gaittime, gaitvararray[,,1,2], cex=1.2)
title("Knee - Hip")

filled.contour(gaittime, gaittime, gaitvararray[,,1,3], cex=1.2)
title("Hip - Hip")

#  plot variance and covariance functions as surfaces

persp(gaittime, gaittime, gaitvararray[,,1,1], cex=1.2)
title("Knee - Knee")

persp(gaittime, gaittime, gaitvararray[,,1,2], cex=1.2)
title("Knee - Hip")

persp(gaittime, gaittime, gaitvararray[,,1,3], cex=1.2)
title("Hip - Hip")

#par(mfrow=c(1,1), mar=c(3,4,2,1), pty="m")

#  plot correlation functions as contours

gaitCorArray <- cor.fd(gaittime, gaitfd)

quantile(gaitCorArray)

contour(gaittime, gaittime, gaitCorArray[,,1,1], cex=1.2)
title("Knee - Knee")

contour(gaittime, gaittime, gaitCorArray[,,1,2], cex=1.2)
title("Knee - Hip")

contour(gaittime, gaittime, gaitCorArray[,,1,3], cex=1.2)
title("Hip - Hip")


#  --------------------------------------------------------------
#            Principal components analysis
#  --------------------------------------------------------------

#  do the PCA with varimax rotation

# Smooth with lambda as determined above
gaitfdPar  <- fdPar(gaitbasis, harmaccelLfd, lambda=1e-2)
gaitpca.fd <- pca.fd(gaitfd, nharm=4, gaitfdPar)

gaitpca.fd <- varmx.pca.fd(gaitpca.fd)

#  plot harmonics using cycle plots

#par(mfrow=c(1,1), mar=c(3,4,2,1), pty="s")
op <- par(mfrow=c(2,2))
plot.pca.fd(gaitpca.fd, cycle=TRUE)
par(op)

#  --------------------------------------------------------------
#           Canonical correlation analysis
#  --------------------------------------------------------------

hipfd  <- gaitfd[,1]
kneefd <- gaitfd[,2]

hipfdPar  <- fdPar(hipfd,  harmaccelLfd, 1e2)
kneefdPar <- fdPar(kneefd, harmaccelLfd, 1e2)

ccafd    <- cca.fd(hipfd, kneefd, ncan=3, hipfdPar, kneefdPar)

#  plot the canonical weight functions

op <- par(mfrow=c(2,1), mar=c(3,4,2,1), pty="m")
plot.cca.fd(ccafd, cex=1.2)
par(op)

#  display the canonical correlations

round(ccafd$ccacorr[1:6],3)
plot(1:6, ccafd$ccacorr[1:6], type="b")


#  --------------------------------------------------------------
#         Register the angular acceleration of the gait data
#  --------------------------------------------------------------

#  compute the acceleration and mean acceleration

D2gaitfd <- deriv.fd(gaitfd,2)
D2gaitmeanfd  <- mean.fd(D2gaitfd)

#  set up basis for warping function

nwbasis   <- 7
wbasis    <- create.bspline.basis(gaitrange,nwbasis,3)
Warpfd    <- fd(matrix(0,nwbasis,39),wbasis)
WarpfdPar <- fdPar(Warpfd)

#  register the functions

regstr <- register.fd(D2gaitmeanfd, D2gaitfd, WarpfdPar, periodic=TRUE)

xfine        <- seq(0,1,len=101)
D2gaitregfd  <- regstr$regfd
D2gaitregmat <- eval.fd(xfine, D2gaitregfd)
warpfd       <- regstr$Wfd
shift        <- regstr$shift
warpmat      <- eval.monfd(xfine, warpfd)
warpmat      <- warpmat/outer(rep(1,101),warpmat[101,])

print(round(shift,1))
hist(shift)

#  plot warping functions

op <- par(mfrow=c(1,1), mar=c(4,4,2,1), pty="m")
matplot(xfine, warpmat, type="l", xlab="t", ylab="h(t)",
        main="Time warping functions", cex=1.2)
par(op)

#  plot the deformation functions, def(t) = warp(t) - t

defmat <- warpmat
for (i in 1:39)
	defmat[,i] <- warpmat[,i] - xfine - shift[i]

op <- par(mfrow=c(1,1), mar=c(4,4,2,1), pty="m")
matplot(xfine, defmat, type="l", xlab="t", ylab="h(t) - t",
        main="Deformation functions", cex=1.2)
par(op)

#  plot both the unregistered and registered versions of the curves

op <- par(mfrow=c(2,2))
plot(D2gaitfd,    ask=FALSE)
plot(D2gaitregfd, ask=FALSE)
par(op)

#  --------------------------------------------------------------
#              Predict knee angle from hip angle
#             for angle and angular acceleration
#  --------------------------------------------------------------

#  set up the data

hipfd  <- gaitfd[,1]
kneefd <- gaitfd[,2]
ncurve <- dim(kneefd$coefs)[2]

kneemeanfd <- mean(kneefd)

#  define the functional parameter object for regression functions

betafdPar <- fdPar(gaitbasis, harmaccelLfd)
betalist  <- list(betafdPar,betafdPar)

#  ----------  predict knee angle from hip angle --------

conbasis <- create.constant.basis(c(0,20))
constfd  <- fd(matrix(1,1,ncurve), conbasis)

#  set up the list of covariate objects

xfdlist  <- list(constfd, hipfd)

#  fit the current functional linear model

fRegressout <- fRegress(kneefd, xfdlist, betalist)

#  set up and plot the fit functions and the regression functions

kneehatfd   <- fRegressout$yhatfd
betaestlist <- fRegressout$betaestlist

alphafd   <- betaestlist[[1]]$fd
hipbetafd <- betaestlist[[2]]$fd

op <- par(mfrow=c(2,1), ask=FALSE)
plot(alphafd,   ylab="Intercept")
plot(hipbetafd, ylab="Hip coefficient")
par(op)

#  compute and plot squared multiple correlation function

gaitfine    <- seq(0,20,len=101)
kneemat     <- eval.fd(gaitfine, kneefd)
#kneehatmat  <- eval.fd(gaitfine, kneehatfd)
kneehatmat  <- predict(kneehatfd, gaitfine)
kneemeanvec <- as.vector(eval.fd(gaitfine, kneemeanfd))

SSE0 <- apply((kneemat - outer(kneemeanvec, rep(1,ncurve)))^2, 1, sum)
SSE1 <- apply((kneemat - kneehatmat)^2, 1, sum)
Rsqr <- (SSE0-SSE1)/SSE0

op <- par(mfrow=c(1,1),ask=FALSE)
plot(gaitfine, Rsqr, type="l", ylim=c(0,0.4))

#  for each case plot the function being fit, the fit,
#                     and the mean function

op <- par(mfrow=c(1,1),ask=TRUE)
for (i in 1:ncurve) {
  plot( gaitfine, kneemat[,i], type="l", lty=1, col=4, ylim=c(0,80))
  lines(gaitfine, kneemeanvec,           lty=2, col=2)
  lines(gaitfine, kneehatmat[,i],        lty=3, col=4)
  title(paste("Case",i))
}
par(op)

#  ----------  predict knee acceleration from hip acceleration --------

D2kneefd     <- deriv(kneefd, 2)
D2hipfd      <- deriv(hipfd, 2)
D2kneemeanfd <- mean(D2kneefd)

#  set up the list of covariate objects

D2xfdlist  <- list(constfd,D2hipfd)

#  fit the current functional linear model

D2fRegressout <- fRegress(D2kneefd, D2xfdlist, betalist)

#  set up and plot the fit functions and the regression functions

D2kneehatfd   <- D2fRegressout$yhatfd
D2betaestlist <- D2fRegressout$betaestlist

D2alphafd   <- D2betaestlist[[1]]$fd
D2hipbetafd <- D2betaestlist[[2]]$fd

op <- par(mfrow=c(2,1), ask=FALSE)
plot(D2alphafd,   ylab="D2Intercept")
plot(D2hipbetafd, ylab="D2Hip coefficient")
par(op)

#  compute and plot squared multiple correlation function

D2kneemat     <- eval.fd(gaitfine, D2kneefd)
#D2kneehatmat  <- eval.fd(gaitfine, D2kneehatfd)
D2kneehatmat  <- predict(D2kneehatfd, gaitfine)
D2kneemeanvec <- as.vector(eval.fd(gaitfine, D2kneemeanfd))

D2SSE0 <- apply((D2kneemat - outer(D2kneemeanvec, rep(1,ncurve)))^2, 1, sum)
D2SSE1 <- apply((D2kneemat - D2kneehatmat)^2, 1, sum)
D2Rsqr <- (D2SSE0-D2SSE1)/D2SSE0

par(mfrow=c(1,1),ask=FALSE)
plot(gaitfine, D2Rsqr, type="l", ylim=c(0,0.5))

#  for each case plot the function being fit, the fit, and the mean function

op <- par(mfrow=c(1,1),ask=TRUE)
for (i in 1:ncurve) {
  plot( gaitfine, D2kneemat[,i], type="l", lty=1, col=4, ylim=c(-20,20))
  lines(gaitfine, D2kneemeanvec,           lty=2, col=2)
  lines(gaitfine, D2kneehatmat[,i],        lty=3, col=4)
  lines(c(0,20), c(0,0), lty=2, col=2)
  title(paste("Case",i))
}
par(op)





