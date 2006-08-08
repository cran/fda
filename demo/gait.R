#  -----------------------------------------------------------------------
#                        Gait data
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  The gait data were chosen for these sample analyses because they are
#  bivariate:  consisting of both hip and knee angles observed over a
#  gait cycle for 39 children.  The bivariate nature of the data implies
#  certain displays and analyses that are not usually considered, and
#  especially the use of canonical correlation analysis.
#  As with the daily weather data, the harmonic acceleration roughness
#  penalty is used throughout since the data are periodic with a strong
#  sinusoidal component of variation.
#  After setting up the data, smoothing the data using GCV or generalized
#  cross-validation to select a smoothing parameter, and displaying
#  various descriptive results, the data are subjected to a principal
#  components analysis, followed by a canonical correlation analysis of the
#  joint variation of hip and knee angle, and finally a registration of the
#  curves.  The registration is included here especially because the 
#  registering of periodic data requires the estimation of a phase shift
#  constant for each curve in addition to possible nonlinear 
#  transformations of time.
#
#  -----------------------------------------------------------------------

#  Last modified 27 March 2006

#  attach the FDA functions

#  Windows ... R

attach("c:\\Program Files\\R\\R-2.2.0\\fdaR\\R\\.RData")

#  Windows ... S-PLUS

attach("c:\\Program Files\\Insightful\\Splus70\\fdaS\\functions\\.Data")


#  -------------  input the data for the two measures  ---------------

hip  <- matrix(scan("../data/hip.txt", 0), 20, 39)
knee <- matrix(scan("../data/knee.txt",0), 20, 39)

#  set up the argument values

gaittime  <- (1:20)/21
gaitrange <- c(0,1)

#  set up a three-dimensional array of function values

gaitarray <- array(0, c(20, 39, 2))
dimnames(gaitarray) <- list(NULL, NULL, c("Hip Angle", "Knee Angle"))
gaitarray[,,1] <- hip
gaitarray[,,2] <- knee

# -----------  set up the harmonic acceleration operator  ----------

harmaccelLfd <- vec2Lfd(c(0, 0, (2*pi)^2, 0))

#  Set up basis for representing gait data.  The basis is saturated since
#  there are 20 data points per curve, and this set up defines 21 basis
#  functions. Recall that a fourier basis has an odd number of basis functions.

gaitbasis <- create.fourier.basis(gaitrange,21)

#  -----------  create the fd object   -----------

#                 Choose level of smoothing using
#          the generalized cross-validation criterion

#  set up range of smoothing parameters in log_10 units

loglam <- seq(-13,-10,0.5)
nlam   <- length(loglam)

dfsave  <- rep(0,nlam)
gcvsave <- rep(0,nlam)

#  loop through smoothing parameters

for (ilam in 1:nlam) {
	lambda <- 10^loglam[ilam]
	fdParobj <- fdPar(gaitbasis, harmaccelLfd, lambda)
	smoothlist <- smooth.basis(gaittime, gaitarray, fdParobj)
	df  <- smoothlist$df
	gcv <- smoothlist$gcv
	dfsave[ilam]  <- df
	gcvsave[ilam] <- sum(gcv)  # note: gcv is a matrix in this case
}

#  display and plot GCV criterion and degrees of freedom

cbind(loglam, dfsave, gcvsave)

#  set up plotting arrangements for one and two panel displays allowing
#  for larger fonts

par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
plot(loglam, gcvsave, type="b",
     xlab="Log_10 lambda", ylab="GCV Criterion", 
     main="Gait Smoothing")

plot(loglam, dfsave, type="b",
     xlab="Log_10 lambda", ylab="Degrees of freedom", 
     main="Gait Smoothing")

#  lambda value minimizing GCV is 10^(-11.5)

lambda    <- 10^(-11.5)
gaitfdPar <- fdPar(gaitbasis, harmaccelLfd, lambda)

gaitfd <- smooth.basis(gaittime, gaitarray, gaitfdPar)$fd

names(gaitfd$fdnames) = c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] = c("Hip", "Knee")

#  --------  plot curves and their first derivatives  ----------------

par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
plot(gaitfd, cex=1.2)

#  plot each pair of curves interactively

par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
plotfit.fd(gaitarray, gaittime, gaitfd, cex=1.2)

#  plot the residuals, sorting cases by residual sum of squares

par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
plotfit.fd(gaitarray, gaittime, gaitfd, residual=T, sort=T, cex=1.2)

#  plot first derivative of all curves

par(mfrow=c(1,2), mar=c(3,4,2,1), pty="s")
plot(gaitfd, Lfdob=1)

#  ------------  compute the mean functions  --------------------

gaitmeanfd <- mean.fd(gaitfd)

#  plot these functions and their first two derivatives

par(mfcol=c(2,3),pty="s")
plot(gaitmeanfd)
plot(gaitmeanfd, Lfdobj=1)
plot(gaitmeanfd, Lfdobj=2)

#  --------------  PCA of gait data  --------------------

#  do the PCA with varimax rotation

gaitpca.fd <- pca.fd(gaitfd, nharm=4, gaitfdPar)

gaitpca.fd <- varmx.pca.fd(gaitpca.fd)

#  plot harmonics using cycle plots

par(mfrow=c(1,1), mar=c(3,4,2,1), pty="s")
plot.pca.fd(gaitpca.fd, cycle=T)

#  ------  Canonical correlation analysis of knee-hip curves  ------

hipfd  <- gaitfd[,1]
kneefd <- gaitfd[,2]

hipfdPar  <- fdPar(hipfd,  harmaccelLfd, 1e-8)
kneefdPar <- fdPar(kneefd, harmaccelLfd, 1e-8)

ccafd    <- cca.fd(hipfd, kneefd, ncan=3, hipfdPar, kneefdPar)

#  plot the canonical weight functions

#  each harmonic in turn

par(mfrow=c(2,1), mar=c(3,4,2,1), pty="m")
plot.cca.fd(ccafd, cex=1.2)

#  display the canonical correlations

round(ccafd$ccacorr[1:6],3)

%  ----------  compute the variance and covariance functions  -------

gaitvarbifd = var.fd(gaitfd)

gaitvararray = eval.bifd(gaittime, gaittime, gaitvarbifd)

par(mfrow=c(1,1), mar=c(3,4,2,1), pty="m")

#  plot variance and covariance functions as contours

contour(gaittime, gaittime, gaitvararray[,,1,1], cex=1.2)
title("Knee - Knee")

contour(gaittime, gaittime, gaitvararray[,,1,2], cex=1.2)
title("Knee - Hip")

contour(gaittime, gaittime, gaitvararray[,,1,3], cex=1.2)
title("Hip - Hip")

#  plot variance and covariance functions as surfaces

persp(gaittime, gaittime, gaitvararray[,,1,1], cex=1.2)
title("Knee - Knee")

persp(gaittime, gaittime, gaitvararray[,,1,2], cex=1.2)
title("Knee - Hip")

persp(gaittime, gaittime, gaitvararray[,,1,3], cex=1.2)
title("Hip - Hip")

#  ----  Register the first derivative of the gait data  --------

#  set up basis for warping function

nwbasis <- 7
wbasis  <- create.bspline.basis(gaitrange,nwbasis,5)

Dgaitfd <- deriv.fd(gaitfd)

Dgaitmeanfd  <- mean.fd(Dgaitfd)

lambda <- 1e-2
Wfd    <- fd(matrix(0,nwbasis,39),wbasis)
WfdPar <- fdPar(Wfd, 3, lambda)

regstr <- register.fd(Dgaitmeanfd, Dgaitfd, WfdPar, periodic=T)

xfine       <- seq(0,1,len=101)
Dgaitregfd  <- regstr$regfd
Dgaitregmat <- eval.fd(xfine, Dgaitregfd)
warpfd      <- regstr$Wfd
shift       <- regstr$shift
warpmat     <- eval.monfd(xfine, warpfd)
warpmat     <- warpmat/outer(rep(1,101),warpmat[101,])

par(mfrow=c(1,1), mar=c(4,4,2,1), pty="m")
matplot(xfine, warpmat, type="l", xlab="t", ylab="h(t)", 
        main="Time warping functions", cex=1.2)

#  plot both the unregistered and registered versions of the curves

par(mfrow=c(2,2))
plot(Dgaitfd,    ask=F)
plot(Dgaitregfd, ask=F)

#  plot the deformation functions, def(t) = warp(t) - t

defmat <- warpmat
for (i in 1:39) 
	defmat[,i] <- warpmat[,i] - xfine - shift[i]

par(mfrow=c(1,1), mar=c(4,4,2,1), pty="m")
matplot(xfine, defmat, type="l", xlab="t", ylab="h(t) - t", 
        main="Deformation functions", cex=1.2)



