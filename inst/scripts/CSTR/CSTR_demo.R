##
##%%  1.  Introduction
##
#%  Estimate the four parameters defining two differential equations
#%  that describe the behavior of the output concentration
#%  and temperature for a nonisothermal continuously stirred
#%  tank reactor, abbreviated nonisothermal CSTR.
#%  The system of two nonlinear equations has five forcing or
#%  input functions.
#%  These equations are taken from
#%  Marlin, T. E. (2000) Process Control, 2nd Edition, McGraw Hill,
#%  pages 899-902.

#%  Last modified 2009.05.31 by Spencer Graves
#%  Matlab version previously modified 26 October 2005

# The R version of this script was translated from Matlab
# and still retains some of the Matlab commands in comments.
# It also includes commands in comments that can be used
# to transfer data between Matlab and R for comparison.

library(fda)
set.seed(1) # so the random numbers will replicate
##
##%%  2.  Set up the problem
##
#%  set up the system of equations with known constants
#%  and values of the four parameters

#%  store known system constants in a struct variable fitstruct

fitstruct <- list(V    = 1.0,#;       %  volume in m^3
                  Cp   = 1.0,#;       %  concentration in cal/(g.K)
                  rho  = 1.0,#;       %  density in g/m^3
                  delH = -130.0,#;    %  cal/kmol
                  Cpc  = 1.0,#;       %  concentration in cal/(g.K)
                  rhoc = 1.0,#;       %  cal/kmol
                  Tref = 350)#;       % reference temperature

#% store true values of known parameters

EoverRtru = 0.83301#;     %  E/R in units K/1e4
kreftru   = 0.4610 #      %  reference value
atru      = 1.678#;       %  a in units (cal/min)/K/1e6
btru      = 0.5#;         %  dimensionless exponent

#% enter these parameter values into fitstruct

fitstruct$kref   = kreftru#;
fitstruct$EoverR = EoverRtru#;   %  kref = 0.4610
fitstruct$a      = atru#;        %  a in units (cal/min)/K/1e6
fitstruct$b      = btru#;        %  dimensionless exponent

##
##%% 3.  Choose the input design
##
#%  Outputs have oscillating responses to lowering
#%  the temperature of the coolant when the tank is
#%  operating at a high temperature, around 365 deg K.
#%  Outputs respond smoothly to changes when the tank
#%  is operating at a cooler temperature, around 330 deg K.

#%  These two designs differ in terms of the temperature
#%  at which the reactor is operated.  The "cool" design
#%  using a coolant temperature of 330 deg K, and the
#%  outputs respond smoothly to step changes in input.
#%  The "hot" design uses a coolant temperature of 365 deg K,
#%  and both concentration and temperature show sharp
#%  oscillations when an input is changed, and especially
#%  when coolant temperature is dropped, at which point
#%  the reactor is nearly unstable.

#%  set up time values at which output values are to be
#%  computed in the numerical solution to the equation.

Tlim  = 64#;   % reaction observed over interval [0, Tlim]
delta = 1/12#; % observe every five seconds
#tspan = (0:delta:Tlim)#;
tspan = seq(0, Tlim, delta)#;
nspan = length(tspan)#;

##
##%%  4.  Display the two input conditions
##
#%  set up input function data for both cool and hot
#%  operations, and plot inputs

coolStepInput <- CSTR2in(tspan, 'all.cool.step')
hotStepInput <- CSTR2in(tspan, 'all.hot.step')

#%  plot inputs
# figure(1)
op <- par(mfrow=c(5, 1), mar=c(2, 4, 2, 1)+.1)

plot(tspan, coolStepInput[, "F."], type="l",
     main="Input flow rate", xlab="", ylab=expression(F(t)),
     axes=FALSE)
axis(1)
axis(2, c(.5, 1., 1.5), las=1)

plot(tspan, coolStepInput[, "CA0"], type="l",
     main="Input concentration", xlab="",
     ylab=expression(C[0](t)), axes=FALSE)
axis(1)
axis(2, c(1.8, 2, 2.2), las=1)

plot(tspan, coolStepInput[, "T0"], type="l",
     main="Input temperature", xlab="",
     ylab=expression(T[0](t)), axes=FALSE)
axis(1)
axis(2, c(303, 323, 343), las=1)

T.cin.lvls = table(c(coolStepInput[, "Tcin"], hotStepInput[, "Tcin"]))
T.cin.rng = range(coolStepInput[, "Tcin"], hotStepInput[, "Tcin"])

plot(range(tspan), T.cin.rng, type="n",
     main="Coolant temparature (red = hot, blue = cool)",
     xlab="", ylab=expression(T[cin](t)), axes=FALSE)
lines(tspan, coolStepInput[, "Tcin"], col="blue")
lines(tspan, hotStepInput[, "Tcin"], col="red")
axis(1)
axis(2, c(330, 350, 370), las=1)

plot(tspan, coolStepInput[, "Fc"], type="l",
     main="Coolant flow rate", xlab="",
     ylab=expression(F[c](t)), axes=FALSE)
axis(1)
axis(2, c(10, 15, 20), las=1)

par(op)

# Matlab:  32 lines;  R:  31 lines

##
##%%  5.  Solve the equations for both the hot and cool conditions
##
#%  5.0.  set constants for ODE solver

#odeoptions = odeset('RelTol',1e-7,'AbsTol',1e-7);

#%  5.1.  cool condition solution

#%  initial conditions

Cinit_cool = 1.5965#;    %  initial concentration in kmol/m^3
Tinit_cool = 341.3754#;  %  initial temperature in deg K
#yinit = [Cinit_cool, Tinit_cool];
yinitCool = c(Conc = Cinit_cool, Temp=Tinit_cool)

#%%  load cool input into fitstruct

fitStrCool <- fitstruct
fitStrCool$Tcin = coolStepInput[, "Tcin"]

#%  5.2.  solve  differential equation with true parameter values

#h = waitbar(0,'Simulating Cool Output...');
#[t, y_cool] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
#                            fitstruct, 'all_cool_step', Tlim);
#close(h)

#library(odesolve)
library(deSolve)
coolStepSoln <- lsoda(y=yinitCool, times=tspan, func=CSTR2, parms=
  list(fitstruct=fitStrCool, condition='all.cool.step', Tlim=Tlim) )

#%  5.3.  set up separate variables for concentration and temperature

C_cool = coolStepSoln[, "Conc"]
plot(tspan, C_cool)

T_cool = coolStepSoln[, "Temp"]
plot(tspan, T_cool)
plot(tspan[-1], diff(T_cool))

#%  5.4.  hot condition solution

#%  initial conditions

Cinit_hot  = 0.2651#;    %  initial concentration in kmol/m^3
Tinit_hot  = 394.0532#;  %  initial temperature in deg K
yinitHot = c(Conc=Cinit_hot, Temp=Tinit_hot);

#%  load hot input into fitstruct

fitStrHot <- fitstruct
fitStrHot$Tcin = hotStepInput[, 'Tcin']

#%  solve  differential equation with true parameter values
#h = waitbar(0,'Simulating Hot Output...');
#[t, y_hot] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
#                            fitstruct, 'all_hot_step', Tlim);
#close(h)
#library(odesolve)
hotStepSoln <- lsoda(y=yinitHot, times=tspan, CSTR2, parms=
  list(fitstruct=fitStrHot, condition='all.hot.step', Tlim=Tlim) )

#%  set up separate variables for concentration and temperature

C_hot = hotStepSoln[, "Conc"]
T_hot = hotStepSoln[, "Temp"]

#% 5.5.  plot deterministic model outputs
# figure(2)
op <- par(mfrow=c(2,1))
t.rng <- range(tspan)
C.rng <- range(C_cool, C_hot)
plot(t.rng, C.rng, type="n", bty='n',
     main="Output concentration (red = hot, blue = cool)",
     xlab="", ylab="C(t)", las=1)
lines(tspan, C_cool, col="blue")
lines(tspan, C_hot, col="red")
abline(v=c(44, 48), lty="longdash")

T.rng <- range(T_cool, T_hot)
plot(t.rng, T.rng, type="n", bty='n',
     main="Output temperature (red = hot, blue = cool)",
     xlab="", ylab="T(t)", las=1)
lines(tspan, T_cool, col="blue")
lines(tspan, T_hot, col="red")
abline(v=c(44, 48), lty="longdash")

par(op)
##
##%%  6.  Set up the B-spline basis for representing solutions.
##
#%  The splines must have discontinuous derivatives and
#%  continuous function values at points where there are
#%  step changes in an input.

#%  this set up is for the cool operation.  The hot
#%  operation requires twice as many knots.

#brk    = (4:4:60);  %  times of step changes
brk <- seq(4, 60, 4)
#knots  = 0:1/3:Tlim;  %  knots at every 20 secs.
knots0 = seq(0, Tlim, 1/3)
knots  = sort(c(knots0, brk, brk));
nknots = length(knots)#;  %  223 knots

#%  set up the basis for the expansion as an order 4
#%  B-spline basis with these knots

#library(fda)
norder = 4;
nbasis = length(knots) + norder - 2;
CSTRbasis0 = create.bspline.basis(c(0,Tlim),nbasis,norder,knots);

#%  set up quadrature points and weights for estimating
#%  the integrals in the roughness penalties

nquadint = 5;
#[CSTRbasis, quadpts, quadwts] = quadset(nquadint, CSTRbasis);
CSTRbasis <- quadset(nquadint, CSTRbasis0)
#str(CSTRbasis)
#all.equal(quadpts, CSTRbasis$quadvals[, "quadpts"])

nquad <- dim(CSTRbasis$quadvals)[1]

#nquad = length(quadpts);
# 960;  someone once recored 4032 quadrature points,
# but I got 960 with both R and Matlab

##
##%%  7.  Compute fits to error less data and display results
##
#%  7.1.  compute the expansion for the errorless data

fdParobj  = fdPar(CSTRbasis, 2, 1e-8);

Cfdtru = smooth.basis(tspan, C_cool, fdParobj);
Tfdtru = smooth.basis(tspan, T_cool, fdParobj);

# Compare with Matlab
# library(R.matlab)
# fdtru <- readMat("CSTRfdtru.mat")
# str(fdtru)
# str(Cfdtru)
# str(Tfdtru)
# sqrt(mean(Cfdtru$coef^2)) = 1.57
# quantile(Cfdtru$coef - fdtru$Cfdtru$coef)
#      0%      25%      50%      75%     100%
#-1.4e-05 -1.3e-06  9.4e-08  1.5e-06  1.3e-05
# sqrt(mean(Tfdtru$coef^2)) = 342
# quantile(Tfdtru$coef - fdtru$Tfdtru$coef)
#      0%      25%      50%      75%     100%
#-1.6e-03 -7.4e-05 -2.2e-06  4.7e-05  1.3e-03

#%  7.2.  plot the fit to the solution using the expansion

#figure(3)
#subplot(2,1,1)
op <- par(mfrow=c(2,1))

plotfit.fd(C_cool, tspan, Cfdtru$fd, titles="Concentration",
           xlab="", sub="", ylab="C(t)")
plotfit.fd(T_cool, tspan, Tfdtru$fd, titles="Temperature",
           xlab="", sub="", ylab="T(t)")
par(op)

#%  7.3.  compute weights to correct for differences in
#%  scale for the two variables

CvectruSpan = eval.fd(tspan, Cfdtru$fd);

Cwt     = var(CvectruSpan);

TvectruSpan = eval.fd(tspan, Tfdtru$fd);

Twt     = var(TvectruSpan);

wt <- c(Conc=Cwt, Temp=Twt)

##
##%%   8.  Set up some data with noise
##
#tobs = 1:1/3:Tlim;  %  observations at every 20 secs.
tobs = seq(1, Tlim, 1/3);#  %  observations at every 20 secs.
N = length(tobs);#   %  190 observations

#%  8.1.  fix standard errors

stderrfac = 0.2#;  %  standard error of measurements
stderrC   = stderrfac*sqrt(Cwt)#;
stderrT   = stderrfac*sqrt(Twt)#;

#%  8.2.  compute errorless values

Cvec = eval.fd(tobs, Cfdtru$fd);
Tvec = eval.fd(tobs, Tfdtru$fd);

#%  8.3.  add errors

Cobs    = Cvec + rnorm(N)*stderrC;
Tobs    = Tvec + rnorm(N)*stderrT;

##****
##**** IF YOU WANT TO COMPARE ANSWERS WITH MATLAB,
##**** TRANSFER MATLAB SIMULATIONS TO ENSURE COMPARABILITY
##****
## You can still get good answers from the R code,
## but you won't be able to compare them with Matlab.
##
#library(R.matlab)
#confirm <- readMat("~R\library\fda\scripts\CSTR\CSTRsim.mat")
#confirm <- readMat("CSTRsim.mat")
#Cobs <- confirm$Cobs
#Tobs <- confirm$Tobs

#%  8.4.  smooth the noisy data using smoothing splines

lambdaC = 1e0;
CfdPar  = fdPar(CSTRbasis, 2, lambdaC);
Cfdsmth = smooth.basis(tobs, Cobs, CfdPar);

lambdaT = 1e-1;
TfdPar  = fdPar(CSTRbasis, 2, lambdaT);
Tfdsmth = smooth.basis(tobs, Tobs, TfdPar);

##
## save.image("CSTR-section8-Matlab sims.Rdata")
## load("CSTR-section8-Matlab sims.Rdata")
##

# TRANSFER MATLAB smooth_basis answers to check
#
# library(R.matlab)
# fdsmth <- readMat("CSTRsmth.mat")
# str(fdsmth)
# str(Cfdsmth)
# str(Tfdsmth)
# sqrt(mean(Cfdsmth$coef^2)) = 1.57
# quantile(Cfdsmth$coef - fdsmth$Cfdsmth$coef)
#      0%      25%      50%      75%     100%
#-4.2e-13 -3.0e-14 -6.0e-15  1.7e-14  9.5e-14
# sqrt(mean(Tfdsmth$coef^2)) = 342
# quantile(Tfdsmth$coef - fdsmth$Tfdsmth$coef)
#      0%      25%      50%      75%     100%
#-9.7e-12 -6.3e-13 -5.7e-14  5.1e-13  8.1e-12

#%  8.5.  plot the data

#figure(4)
#subplot(2,1,1)
op <- par(mfrow=c(2,1))

plot(tobs, Cobs, ylab="C(t)", xlab='',
     main="Concentration");
plot(tobs, Tobs, xlab='', ylab="T(t)",
     main="Temperature");
par(op)

save(file="CSTR_fig4.Rdata")
# load("CSTR_fig4.Rdata")

#%  8.6.  plot the data and the smooths

#figure(5)
#subplot(2,1,1)

op <- par(mfrow=c(2,1))
plotfit.fd(Cobs, tobs, Cfdsmth$fd, xlab='', ylab='C(t)',
           titles="Concentration", sub='');
plotfit.fd(Tobs, tobs, Tfdsmth$fd, xlab='', ylab='T(t)',
           titles="Temperature", sub='');
par(op)

##
##%%  9.  Load data into struct variable datstruct
##

#%  this information only depends on the
#%  data and sampling design

#%  9.1.  weights for variables

datstruct0 <- list(Cwt=Cwt, Twt=Twt)

#%  9.2.  data

datstruct0$y = cbind(Conc=Cobs, Temp=Tobs);

#%  9.3.  basis values at sampling points

basismat            = eval.basis(tobs, CSTRbasis);
Dbasismat           = eval.basis(tobs, CSTRbasis, 1);
datstruct0$basismat  = basismat;
datstruct0$Dbasismat = Dbasismat;

#%  9.4.  forcing function values at quadrature points

quadpts <- CSTRbasis$quadvals[, "quadpts"]
coolStepQuadInput <- CSTR2in(quadpts, 'all.cool.step');

datstruct <- c(datstruct0, as.data.frame(coolStepQuadInput))

#%  9.5.  basis values at quadrature points

quadbasismat            = eval.basis(quadpts, CSTRbasis);
Dquadbasismat           = eval.basis(quadpts, CSTRbasis, 1);

datstruct$quadpts       = CSTRbasis$quadvals[, "quadpts"]
datstruct$quadwts       = CSTRbasis$quadvals[, "quadwts"]
datstruct$quadbasismat  = quadbasismat;
datstruct$Dquadbasismat = Dquadbasismat;

##
##  10.  Analyze the noisy data to estimate
##       * a CSTRbasis representation (section 11) and
##       * parameters of the CSTR differential equations (sect. 12-14)
##
#%  This is a very difficult problem because all four
#%  parameters affect the speed of the reaction in nearly
#%  identical ways.  Convergence is slow.

##
##%%  11.  Preliminary explorations of sampling variances
##
#%  Compute a solution using the true parameters as a
#%  starting point and zero iterations

#%  11.1.  Specify that both output variables are measured

fitStrHot$fit <- c(Conc=1,Temp=1)

#%  11.2.  Use coefficients from spline smoothing as initial
#%  values for coefficient vectors

Ccoef <- Cfdsmth$fd$coef
Tcoef <- Tfdsmth$fd$coef

fitStrHot$coef0 <- data.frame(Conc=Ccoef, Temp=Tcoef)
all.equal(fitStrHot$coef0, data.frame(Conc=Ccoef, Temp=Tcoef))
# TRUE

#%  11.3.  Specify that all parameters will be estimated

fitStrHot$estimate <- rep(1,4)

estind <- which(fitStrHot$estimate==1)

#%  11.4.  set up the true parameter values

parvec0 <- c(kref=kreftru, EoverR=EoverRtru,
             a=atru, b=btru);

#%  11.5.  specify smoothing parameter values

lambda = 1e1*c(lambdaC, lambdaT);

#%  11.6.  Obtain coefficients for a solution in CSTRbasis

fitStrHot$fit <- c(Conc=1, Temp=1)
(et.1 <- system.time(
CSTR.1 <- CSTRfn(parvec0, datstruct, fitStrHot, CSTRbasis, lambda) ))
# 27 sec. 2009.05.31

fitStrHot$fit <- c(Conc=0, Temp=1)
(et.01 <- system.time(
CSTR1.01 <- CSTRfn(parvec0, datstruct, fitStrHot, CSTRbasis, lambda) ))
# 25 sec. 2009.05.31

fitStrHot$fit <- c(Conc=1, Temp=0)
(et.10 <- system.time(
CSTR1.10 <- CSTRfn(parvec0, datstruct, fitStrHot, CSTRbasis, lambda) ))
# 38 sec. 2009.05.31

#library(R.matlab)
#CSTR.1Mat <- readMat('CSTR1.mat')
#str(CSTR.1Mat)
#sqrt(mean(CSTR.1Mat$res^2)) # 0.168
#sqrt(mean(CSTR.1$res^2)) # 0.168
# d.res <- (CSTR.1$res - as.vector(CSTR.1Mat$res))
#sqrt(mean(d.res^2)) # 5.46e-7
# quantile(d.res)
#      0%      25%      50%      75%     100%
#-1.9e-06 -2.3e-07  2.3e-08  2.7e-07  1.8e-06
#**GREAT

# sqrt(mean(CSTR.1$Dres^2)) # 3.56
# d.Dres <- (CSTR.1$Dres - CSTR.1Mat$jacobian)
# sqrt(mean(d.Dres^2))# 1.1e-5
# quantile(d.Dres)
#      0%      25%      50%      75%     100%
#-3.7e-05 -5.9e-07  6.3e-08  1.5e-06  1.9e-05

# Check CSTR.1$Dres

fitStrHot$fit <- c(Conc=1, Temp=1)

ncoef2 <- length(CSTR.1$res)
parNames <- dimnames(CSTR.1$Dres)[[2]]
k2 <- length(parNames)
Dres. <- array(NA, dim=c(ncoef2, k2), dimnames=
               list(NULL, parNames) )
delta = 1e-4;
for(i in 1:k2){
  pvi <- parvec0
  pvi[i] <- pvi[i]+delta
  CSTR.i <- CSTRfn(pvi, datstruct, fitStrHot, CSTRbasis, lambda)
  Dres.[,i] <- as.vector(CSTR.i$res-CSTR.1$res)/delta
  cat(i, "")
}
mean(abs(CSTR.1$Dres-Dres.))
#[1] .907
apply(abs(CSTR.1$Dres-Dres.), 2, mean)
#       kref      EoverR           a           b
#0.016452836 0.005014252 0.006450557 3.598695014

#*** CONCLUSION:
#*** The partials with respect to kref, EoverR, and 'a' seem OK
#*** but the partial with respect to 'b' seems mostly error.

apply(abs(CSTR.1$Dres), 2, mean)
#      kref     EoverR          a          b
#6.14653067 1.84677871 0.80659077 0.04626569

# The first 3 columns of CSTR.1$Dres
# are close to the first difference approximation to the
# partial derivatives.
# However, the last one is very different.

#%    get the singular values of the jacobian

#s = svd(full(jacobian));
CSTR.1.Dsvd <- svd(CSTR.1$Dres)
s = CSTR.1.Dsvd$d

#%  log 10 of condition number

log10(s[1]/s[4])
#   3.02
# Matches Matlab 2007.05.29 s. graves

#%  Sampling variances will vary by over 6 orders of magnitude!

#%  11.7.  Take only the first two parameters from parvec,
#     and include derivatives for only those in the Jacobian, Dres

fitStrH12 = fitStrHot
fitStrH12$estimate = c(1, 1, 0, 0)
estind12 = which(fitStrH12$estimate == 1);

parvec12 = c(kreftru, EoverRtru)

#%    calculate the estimates

#[res, jacobian] = CSTRfn(parvec0, datstruct, fitstruct, ...
#                         CSTRbasis, lambda);
(et.2 <- system.time(
CSTR.2 <- CSTRfn(parvec12, datstruct, fitStrH12, CSTRbasis, lambda) ))
# 26 sec.  2009.05.31

#str(CSTR.1)
#str(CSTR.2)
all.equal(CSTR.1$res, CSTR.2$res)
# TRUE
all.equal(CSTR.1$Dres[, 1:2], CSTR.2$Dres)
# TRUE
# library(R.matlab)
#CSTR.2Mat <- readMat('CSTR12.mat')
#str(CSTR.2Mat)
#sqrt(mean(CSTR.2$res^2)) # 0.168
# d.res.2 <- (CSTR.2$res - as.vector(CSTR.2Mat$res))
# sqrt(mean((d.res.2^2) )) # 5.46e-7
# quantile(d.res.2)
#           0%           25%           50%           75%          100%
#-1.902187e-06 -2.333586e-07  2.302885e-08  2.729940e-07  1.844674e-06
#**GREAT

# sqrt(mean(CSTR.2$Dres^2)) # 5.0
# d.Dres.2 <- (CSTR.2$Dres - CSTR.2Mat$jacobian)
# sqrt(mean(d.Dres.2^2))# 1.6e-5
# quantile(d.Dres.2)
#         0%         25%         50%         75%        100%
#-3.706184e-05 -1.981728e-05 -2.700406e-07  2.673825e-06  1.853554e-05

#%  get the singular values of the jacobian

s12 = svd(CSTR.2$Dres)$d

log10(s12[1]/s12[2])
# R = 0.826  same as Matlab

#%  11.8.  Use only the first and third parameters
#     and include derivatives for only those in the Jacobian, Dres

fitStrH13 = fitStrHot
fitStrH13$estimate = c(1, 0, 1, 0)
estind13 = which(fitStrH13$estimate == 1);

parvec13 = c(kreftru, atru)

#%  calculate the estimates

# CSTRfn use 3
#[res, jacobian] = CSTRfn(parvec0, datstruct, fitstruct, ...
#                         CSTRbasis, lambda);

(et.13 <- system.time(
CSTR.13 <- CSTRfn(parvec13, datstruct, fitStrH13, CSTRbasis, lambda) ))
# 26 sec.  2009.05.31

#%  get the singular values of the jacobian

s13 = svd(CSTR.13$Dres)$d

#%  log 10 of condition number
log10(s13[1]/s13[2])
#  1.1789 on 2007.05.30  Matches Matlab

#%  11.9.  Use only the third and fourth parameters
#     and include derivatives for only those in the Jacobian, Dres

fitStrH34 = fitStrHot
fitStrH34$estimate = c(0, 0, 1, 1)
estind = which(fitstruct$estimate == 1);

parvec34 = c(atru, btru)

#%  calculate the estimates

# CSTRfn use 4
#[res, jacobian] = CSTRfn(parvec0, datstruct, fitstruct, ...
#                         CSTRbasis, lambda);
(et.34 <- system.time(
CSTR.34 <- CSTRfn(parvec34, datstruct, fitStrH34, CSTRbasis, lambda) ))
# 26 sec. 2009.05.31

#%  get the singular values of the jacobian

s34 = svd(CSTR.34$Dres)$d

#%  log 10 of condition number
log10(s34[1]/s34[2])
# 2.1200 on 2007.05.30  matches Matlab

#%  11.10.  Use only the first, second and third parameters
#     and include derivatives for only those in the Jacobian, Dres

fitStrH123 = fitStrHot
fitStrH123$estimate = c(1, 1, 1, 0)
estind123 = which(fitstruct$estimate == 1);

parvec123 = c(kreftru, EoverRtru, atru)

#%  calculate the estimates

# CSTRfn use 5
#[res, jacobian] = CSTRfn(parvec0, datstruct, fitstruct, ...
#                         CSTRbasis, lambda);
(et.123 <- system.time(
CSTR.123 <- CSTRfn(parvec123, datstruct, fitStrH123, CSTRbasis, lambda) ))
# 26 sec. 2009.05.31

#CSTR.123Mat <- readMat('CSTR123.mat')
#str(CSTR.123Mat)
#sqrt(mean(CSTR.123$res^2)) # 0.168
# d.res.123 <- (CSTR.123$res - as.vector(CSTR.123Mat$res))
#sqrt(mean(d.res.123^2)) # 5.46e-7
# quantile(d.res.123)
#           0%           25%           50%           75%          100%
#-1.902187e-06 -2.333586e-07  2.302885e-08  2.729940e-07  1.844674e-06
#**GREAT

# sqrt(mean(CSTR.123$Dres^2)) # 4.12
# d.Dres.123 <- (CSTR.123$Dres - CSTR.123Mat$jacobian)
# sqrt(mean(d.Dres.123^2))# 1.3e-5
# quantile(d.Dres.123)
#         0%         25%         50%         75%        100%
#-3.706184e-05 -1.981728e-05 -2.700406e-07  2.673825e-06  1.853554e-05
# Good

#%  get the singular values of the jacobian

s123 = svd(CSTR.123$Dres)$d

#%  log 10 of condition number
log10(s123[1]/s123[3])
# 1.236 on 2007.05.30
# matches Matlab

##
##%%  12.  Estimate all four parameters
##

#%  set the optimization parameters for the outer optimization
#%  only function values are used at this point

tolval = 1e-8;

#%  use coefficients from spline smoothing as initial
#%  values for coefficient vectors

fitStrH1234 = fitStrHot ;
fitStrH1234$coef0 = data.frame(Conc=Ccoef, Temp=Tcoef)

#%  Use only the first, second and third parameters

fitStrH1234$fit = c(1,1)

fitStrH1234$estimate = c(1, 1, 1, 1)
estind1234 = which(fitStrH1234$estimate == 1)

#%  set up some initial values for parameters

parvec1234 = c(kref=0.4, EoverR=0.8, a=1.7, b=0.5)

#%  set up the true parameter values

parvectru = c(kreftru, EoverRtru, atru, btru)

#%  specify smoothing parameter values

lambda = 10.*c(lambdaC, lambdaT)

#%  calculate the estimates

# CSTRfn use 6 (indirectly via CSTRres0)
#tic;
#[parvec, resnorm, residual, exitflag, output, lambdaout, jacobian] = ...
#    lsqnonlin(@CSTRfn, parvec0, [], [], optionsCSTRfn, ...
#              datstruct, fitstruct, CSTRbasis, lambda);
#toc
#library(R.matlab)
#CSTR.mat <- readMat("CSTRfnTst0.mat")

save(file="CSTR-section12-nls-prep.Rdata",
     list=c("parvec0", "parvec1234", "datstruct", "fitStrCool", "fitStrHot",
		"fitStrH1234", "CSTRbasis", "lambda", "btru"))
# load("CSTR-section12-nls-prep.Rdata")
# or
# fdaPath <- system.file(package = 'fda')
# fdaDir <- dir(fdaPath)
# fdaScriptsPath <- dir(fdaPath, full.names=TRUE)[fdaDir == 'scripts']
# fdaScriptsDir <- dir(fdaScriptsPath)
# fdaCSTRpath <- dir(fdaScriptsPath, full.names=TRUE)[fdaScriptsDir=='CSTR']
# fdaCSTRdir <- dir(fdaCSTRpath)
# fdaCSTRdat <- dir(fdaCSTRpath, full.names=TRUE)[grep('\\.Rdata$', fdaCSTRdir)]
# load(fdaCSTRdat)

.datstruct <- datstruct
.fitstruct <- fitStrH1234
.CSTRbasis <- CSTRbasis
.lambda <- lambda

.CSTRres0.trace <- NULL

et.nls0.0 <- system.time(
nlsFit0.0 <- nls(formula=
     ~CSTRres0(kref=kref, EoverR=EoverR, a=a, b=b,gradwrd=FALSE),
    start=as.list(parvec1234),
    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE) )
#

#et.nls0.0 <- system.time(
#nlsFit0.0 <- NLS(formula=
#     ~CSTRres0(kref=kref, EoverR=EoverR, a=a, b=b,gradwrd=FALSE),
#    start=as.list(parvec1234),
#    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE) )
# warnOnly=TRUE to force output:
# CSTR computations are so intense that 'nls' overestimates
#   the numerical precision feasible and often
#   terminates with an error & returns nothing
#   when it can't achieve that level of precision.
#......70.23466 :  0.4 0.8 1.7 0.5 7
#  It.   1, fac=           1, eval (no.,total): ( 1,  1):..... new dev = 10.9019
#10.90186 :  0.4672504 0.8296264 1.9357038 0.4600686
#  It.   2, fac=           1, eval (no.,total): ( 1,  2):..... new dev = 10.6761
#10.67609 :  0.4663919 0.8403480 1.7992337 0.4784262
#  It.   3, fac=           1, eval (no.,total): ( 1,  3):..... new dev = 10.6673
#10.66733 :  0.4661362 0.8397804 1.7173975 0.4961549
#  It.   4, fac=           1, eval (no.,total): ( 1,  4):..... new dev = 10.6667
#10.66668 :  0.4661972 0.8398665 1.7431566 0.4910539
#  It.   5, fac=           1, eval (no.,total): ( 1,  5):..... new dev = 10.6664
#10.66637 :  0.4662310 0.8399787 1.7192274 0.4961931
#  It.   6, fac=           1, eval (no.,total): ( 1,  6):..... new dev = 10.6667
#  It.   6, fac=         0.5, eval (no.,total): ( 2,  7):..... new dev = 10.6665
# ...
#  It.   6, fac=     0.03125, eval (no.,total): ( 6, 11):..... new dev = 10.6664
#10.66637 :  0.4662309 0.8399774 1.7200228 0.4960231
#  It.   7, fac=      0.0625, eval (no.,total): ( 1, 12):..... new dev = 10.6664
#10.66636 :  0.4662251 0.8399652 1.7199096 0.4960453
#  It.   8, fac=       0.125, eval (no.,total): ( 1, 13):..... new dev = 10.6664
#  ...
#  It.   8, fac= 0.000976562, eval (no.,total): ( 8, 20):..... new dev = 10.6664
#Warning messages:
#1: Stepsize reduced below the minimum with parvec = 0.467, 0.83, 1.94, 0.46 on iteration 9 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstr, fitstruct = fitstr, CSTRbasis = CSTRb,
#...
#9: Stepsize reduced below the minimum with parvec = 0.466, 0.84, 1.72, 0.496 on iteration 10 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstr, fitstruct = fitstr, CSTRbasis = CSTRb,
#10: step factor 0.000488281 reduced below 'minFactor' of 0.000976562 in: nls(formula = ~CSTRres0(kref = kref, EoverR = EoverR, a = a,

# NOTE:  The above trace includes one string of 6 periods '.' and 20 strings of 5.
# CSTRres0 produces one period '.'  each time it is called with an object named
# '.CSTRres0.trace' available.  The resulting matrix '.CSTRres0.trace'
# has 106 rows = 6+5*20.

# Why so much trouble with step size reduction, especially
# the final 'nls' "step factor ... reduced below 'minFactor'"
# message?

dim(.CSTRres0.trace)
# 106 385
str(.datstruct)
str(.fitstruct)
str(.CSTRbasis)

dimnames(.CSTRres0.trace) <- list(NULL,
       c("kref", "EoverR", "a", "b", "SSE",
         t(outer(c("C", "T"), 1:190, paste, sep=""))) )

save(file="CSTR-nlsFit0.0.Rdata",
     list=c("nlsFit0.0", ".CSTRres0.trace") )
# load("CSTR-nlsFit0.0.Rdata")

deviance(nlsFit0.0)
#     10.66636
# vs. 10.6666 reported by Matlab

CSTRres0.trace <- cbind(
   iter = rep(0:8, c(6, 5,5,5,5,5, 30, 5, 40)),
   stepFactor = c(1, rep(c(rep(1, 6), 1/2^c(0:5, 4, 3:10)), each=5)),
   evalInGroup= c(1, rep(c(1,1,1,1,1,1, 1:6, 1, 1:8), each=5)),
   totalEvals = c(0, rep(0:20, each=5)),
   .CSTRres0.trace)
dimnames(CSTRres0.trace)

CSTRres0.trace[1:9, 1:11]

plot(CSTRres0.trace[, 1])
lines(8*CSTRres0.trace[, 2])

plot(CSTRres0.trace[, 3])
plot(CSTRres0.trace[, 4])

write.table(CSTRres0.trace, "CSTRres0-trace.csv", sep=",")
# CSTRres0.trace<-as.matrix(read.table("CSTRres0-trace.csv",header=TRUE,sep=","))
str(CSTRres0.trace)

(n1 <- 105/5)
(k <- dim(CSTRres0.trace)[2])
CSTRres0.t3 <- CSTRres0.trace[-1,]
dim(CSTRres0.t3) <- c(5, n1, k)
str(CSTRres0.tr3)

CSTRres0.tr3 <- aperm(CSTRres0.t3, c(2:3, 1) )
dimnames(CSTRres0.tr3) <- list(
      NULL, dimnames(CSTRres0.trace)[[2]], NULL)
str(CSTRres0.tr3)

dimnames(CSTRres0.tr3)[[2]][1:22]

k0 <- 9
N1 <- dim(CSTRres0.tr3)[2]
resi <- (k0+1):N1
Theta <- c("kref", "EoverR", "a", "b")

# Manually check 'nls' computations

RalstonJenrichStep <- function(x=CSTRres0.tr3[1,,],
                    theta=Theta, res=resi){
# Single step computation per
# Ralston and Jenrich (1978) algorithm, as described in
# Bates and Watts (1988)
# Nonlinear Regression Analysis and Its Applications
# (Springer, sec. 3.5.4, p. 84-85)
  H1 <- x[, -1]-x[, 1]
  p <- length(theta)
  d.th <- rep(NA, p)
  for(i in 1:p)
    d.th[i] <- (x[theta[i], i+1]-x[theta[i], 1])
  fit <- lm.fit(H1[res,], x[res,1])
  Dth <- (-d.th*coef(fit))
  names(Dth) <- paste("d", theta, sep=".")
  Dth
}

Theta0 <- paste(Theta, "0", sep="")
Theta.best <- paste(Theta, "best", sep=".")

d.Theta <- paste("d", Theta, sep=".")

k00 <- 4
names0 <- dimnames(CSTRres0.tr3)[[2]]
(names2 <- c(names0[1:k00], Theta0, "SSE0",
             Theta.best, "best.SSE", "best.i", d.Theta))
k1 <- length(names2)
CSTRres0.trSum <- array(NA, dim=c(n1, k1), dimnames = list(
          NULL, names2) )

k.best <- k0+(1:5)
RalstonJenrichStep()
k1. <- max(k.best)+(2:5)

for(i1 in 1:n1){
  CSTRres0.trSum[i1, 1:k0] <- CSTRres0.tr3[i1, 1:k0, 1]
#
  SSEi <- CSTRres0.tr3[i1, "SSE", ]
  bestSSE <- which(SSEi==min(SSEi))[1]
  CSTRres0.trSum[i1, k.best] <- CSTRres0.tr3[i1, 5:9, bestSSE]
  CSTRres0.trSum[i1, "best.i"] <- bestSSE
  CSTRres0.trSum[i1, k1.] <- RalstonJenrichStep(CSTRres0.tr3[i1,,])
}

CSTRres0.trSum[7, ]
j1 <- 9
j1ref <- 6
CSTRres0.trSum[j1ref, Theta0]+.25*CSTRres0.trSum[j1ref,d.Theta]
(th9 <- CSTRres0.trSum[j1ref, Theta0]-.25*CSTRres0.trSum[j1ref,d.Theta])
names(th9) <- Theta
tst <- do.call(CSTRres0, as.list(th9))
sum(tst^2)

j1 <- 11
j1ref <- 6
(CSTRres0.trSum[j1ref, Theta0]+
         CSTRres0.trSum[j1, "stepFactor"]*CSTRres0.trSum[j1ref,d.Theta])
(th11 <- (CSTRres0.trSum[j1ref, Theta0]-
         CSTRres0.trSum[j1, "stepFactor"]*CSTRres0.trSum[j1ref,d.Theta]))
names(th11) <- Theta
res11 <- do.call(CSTRres0, as.list(th11))
sum(res11^2)
sum(res11^2)-10.666


j1 <- 15
j1ref <- 13
(CSTRres0.trSum[j1ref, Theta0]+
         CSTRres0.trSum[j1, "stepFactor"]*CSTRres0.trSum[j1ref,d.Theta])
(th11 <- (CSTRres0.trSum[j1ref, Theta0]-
         CSTRres0.trSum[j1, "stepFactor"]*CSTRres0.trSum[j1ref,d.Theta]))
names(th11) <- Theta
res11 <- do.call(CSTRres0, as.list(th11))
sum(res11^2)
sum(res11^2)-10.666

j1 <- 16
j1ref <- 13
(CSTRres0.trSum[j1ref, Theta0]+
         CSTRres0.trSum[j1, "stepFactor"]*CSTRres0.trSum[j1ref,d.Theta])
(th11 <- (CSTRres0.trSum[j1ref, Theta0]-
         CSTRres0.trSum[j1, "stepFactor"]*CSTRres0.trSum[j1ref,d.Theta]))
names(th11) <- Theta
res11 <- do.call(CSTRres0, as.list(th11))
sum(res11^2)
sum(res11^2)-10.666

theta7 <- CSTRres0.trSum[j1ref, Theta0]
names(theta7) <- Theta
theta7
d.theta7 <- CSTRres0.trSum[j1ref, d.Theta]
theta7+d.theta7
stepFac <- outer(c(1:3, 5, 7), 10^((-3):0))
stepFac. <- sort(c(0, -stepFac, stepFac))

nStep <- length(stepFac.)
CSTRres0.trFinalStep <- array(NA, dim=c(nStep, 6), dimnames=
       list(NULL, c("stepFac", Theta, "SSE")))
CSTRres0.trFinalStep[, "stepFac"] <- stepFac.
for(i in 1:nStep){
  th.i <- (theta7+stepFac.[i]*d.theta7)
  CSTRres0.trFinalStep[i, Theta] <- th.i
  res.i <- do.call(CSTRres0, as.list(th.i))
  CSTRres0.trFinalStep[i, "SSE"] <- sum(res.i^2)
  cat(i, "")
}

apply(CSTRres0.trFinalStep, 2, min)
plot(CSTRres0.trFinalStep[, c(1, 6)])
abline(v=0)
plot(CSTRres0.trFinalStep[6:36, c(1, 6)])
abline(v=0)
plot(CSTRres0.trFinalStep[11:31, c(1, 6)])
abline(v=0)
plot(CSTRres0.trFinalStep[11:27, c(1, 6)])
abline(v=0)

write.table(CSTRres0.trFinalStep, "negative stepFactor required.csv",
            sep=",", row.names=FALSE)

write.table(CSTRres0.trSum, "CSTRres0.trSum.csv", sep=",",
            row.names=FALSE)

stepCompare <- rep(NA, n1-1)
i0 <- 1
for(i1 in 2:n1){
  if(diff(CSTRres0.trSum[c(i0, i1), "SSE"])<0)
    i0 <- i1-1
  stepCompare[i1-1] <- sum(abs(CSTRres0.trSum[i1, Theta] -
         (CSTRres0.trSum[i1-1, Theta]+
         CSTRres0.trSum[i1-1, 'stepFactor']*
                           CSTRres0.trSum[i1-1, d.Theta])) )
}

quantile(stepCompare)
plot(stepCompare)

write.table(
  cbind(CSTRres0.trSum, c(NA, sum.abs.chg.pred.error=stepCompare)),
            "CSTRres0.trSum.csv", sep=",", row.names=FALSE)


write.table(RmatComp, "R-Matlab convergence comparison.csv",
            sep=",", row.names=FALSE)
RmatComp <- read.table("R-Matlab convergence comparison.csv",
       header=TRUE, sep=",", stringsAsFactors=FALSE)
dimnames(RmatComp)[[2]] <- RmatCnames

# CONCLUSION:
# The 'nls' code does NOT check negative step sizes,
# as noted in Bates and Watts (1988)
# Nonlinear Regression Analysis and Its Applications
# (Springer, sec. 3.5.4, p. 84-85)
# citing Ralston and Jenrich (1978).
# The problem is that a second plan does NOT adequately
# approximate the nonlinear manifold.


# More comparisons with Matlab in this example:

library(R.matlab)
# (1) Open 'CSTR_demo.m' in Matlab,
# (2)  Create a connection to Matlab
(matlab <- Matlab())
# (3) In Matlab, addpath to 'MatlabServer'
#     at system.file("externals", package="R.matlab")
# (4) In Matlab > MatlabServer
#     This locks up Matlab until released by close(matlab) below
# (5) Open the connection
(isOpenMatlab <- open(matlab))
# (6) Transfer '.CSTRres0.trace' to Matlab
str(.CSTRres0.trace)
setVariable(matlab, R_CSTRres0trace = .CSTRres0.trace)
# (7) Close the connection
close(matlab)

# (8) In Matlab, create
#     RmatlabCSTRres0trace = [R_CSTRres0trace, NaN]
# (9) Fill the  with CSTRfn(R_CSTRres0trace(i, 1:4), ...)

#(10) In Matlab > MatlabServer % to transfer control back to R.
#(11) Reopen the connection
(isOpenMatlab <- open(matlab))
#(12) Get the results
CSTRres0.R.Matlab.trace <-  getVariable(matlab, 'RmatlabCSTRres0trace')
#(13) Allow further work in Matlab
close(matlab)

dMatR <- (CSTRres0.R.Matlab.trace$RmatlabCSTRres0trace[, 6] -
  CSTRres0.R.Matlab.trace$RmatlabCSTRres0trace[, 5])

RmatComp <- cbind(
     rep(0:8, c(6, 5,5,5,5,5, 30, 5, 40)),
     c(1, rep(c(rep(1, 6), 1/2^c(0:5, 4, 3:10)), each=5)),
     c(0, rep(0:20, each=5)), c(0, rep(1:5, 21)),
     CSTRres0.R.Matlab.trace$RmatlabCSTRres0trace, dMatR,
     dMatR / CSTRres0.R.Matlab.trace$RmatlabCSTRres0trace[, 5] )

RmatCnames <- c('iter', 'stepFactor', 'group', 'evalInGroup',
    "kref", "EoverR", "a", "b", "R.SSE",
    'Matlab.SSE', 'Matlab-R.SSE', 'rel.dSSE((M-R)/R)')
dimnames(RmatComp) <- list(NULL, RmatCmanes)
class(RmatComp)

write.table(RmatComp, "R-Matlab convergence comparison.csv",
            sep=",", row.names=FALSE)
RmatComp <- read.table("R-Matlab convergence comparison.csv",
       header=TRUE, sep=",", stringsAsFactors=FALSE)
dimnames(RmatComp)[[2]] <- RmatCnames

plot(RmatComp[,'Matlab-R.SSE'])
plot(RmatComp[-(1:6),'Matlab-R.SSE'])
plot(RmatComp[-(1:16),'Matlab-R.SSE'])
plot(RmatComp[, 'rel.dSSE((M-R)/R)'] )
plot(RmatComp[-(1:6), 'rel.dSSE((M-R)/R)'] )

# Conclusion:  The R code seems as good as the Matlab.

RmCdim <- dim(RmatComp)
# 106 12
(n1 <- 105/5)
RmatC3 <- array(NA, dim=c(n1, 12, 5), dimnames=list(
      NULL, RmatCnames, c("center", "d.kref", "d.EoverR", "d.a", "d.b")))
for(i in 1:n1)
  for(j in 1:5){
    i1 <- 1+j+5*(i-1)
    RmatC3[i,,j] <- unlist(RmatComp[i1,])
}

RmatComp[1:8, ]
RmatC3[1:2,,]

table(RmatComp[,1])
# 0  1  2  3  4  5  6  7  8
# 6  5  5  5  5  5 30  5 40

RmatComp8 <- RmatComp[67:106, -1]
RmatComp[66:68,]
RmatComp8[1:2,]
dim(RmatComp8)
# 40 11

plot(RmatComp8[, "kref"])
plot(RmatComp8[, "EoverR"])
plot(RmatComp8[, "a"])
plot(RmatComp8[, "b"])

dRmatComp8 <- RmatComp8-rep(c(0,0,0, unlist(RmatComp[62,-(1:4)])), each=40)
plot(kref~stepFactor, dRmatComp8)
plot(EoverR~stepFactor, dRmatComp8)
plot(a~stepFactor, dRmatComp8)
plot(b~stepFactor, dRmatComp8)

plot(R.SSE~stepFactor, dRmatComp8)
plot(Matlab.SSE~stepFactor, dRmatComp8)

plot(dRmatComp8[, c("stepFactor", "Matlab-R.SSE")])
plot(dRmatComp8[, c("stepFactor", "rel.dSSE((M-R)/R)")])

plot(RmatComp[c( 67:106), c("kref", "EoverR")])
plot(RmatComp[c(62, 67:106), c("kref", "EoverR")])
plot(RmatComp[c(62, 67:106), c("kref", "R.SSE")])
plot(RmatComp[c(62, 67:106), c("kref", "Matlab.SSE")])
plot(RmatComp8[ , c("stepFactor", "Matlab.SSE")])

plot(RmatComp8[ , c("stepFactor", "R.SSE")], bty="n")

# Return to nlsFit0.0 example
coef(nlsFit0.0)
#     kref    EoverR         a         b
#0.4661415 0.8397917 1.7182484 0.4963714
# vs.
#0.4662    0.8396    1.7014    0.5002 per Matlab

# nls does not return the Jacobian, so call CSTFfn directly

nlsFit0.0a <- CSTRfn(parvec=coef(nlsFit0.0), datstruct, fitStrH1234,
                     CSTRbasis, lambda)
fitStrH1234.10 <- fitStrH1234
fitStrH1234.10$fit <- c(Conc=1, Temp=0)
nlsFit0.0a10 <- CSTRfn(parvec=coef(nlsFit0.0), datstruct, fitStrH1234.10,
                     CSTRbasis, lambda)

fitStrH1234.01 <- fitStrH1234
fitStrH1234.01it <- c(Conc=0, Temp=1)
nlsFit0.0a01 <- CSTRfn(parvec=coef(nlsFit0.0), datstruct, fitStrH1234.01,
                     CSTRbasis, lambda)

(s1234 <- svd(nlsFit0.0a$Dres)$d)
# log10 of the condition number
log10(s1234[1]/s1234[4])
# 3.0173 vs. 3.0146 for the Matlab solution

###################################################
#
#  There seems to be a problem with the gradient,
#  because 'nls' won't converge when I try to use it.
#  Moreover, when I compare it with a first difference approximation,
#  the first three columns match, but the fourth does not.

#nlsFit0g <- nls(formula=
#     ~CSTRres0(kref=kref, EoverR=EoverR, a=a, b=b,gradwrd=FALSE),
#    start=as.list(parvec0),
#    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE)

###################################################
# Try passing arguments via a list 'data'

et.nlsFit <- system.time(
nlsFit1234o <- NLS(formula=~CSTRres(kref=kref, EoverR=EoverR, a=a, b=b,
   datstruct=datstr, fitstruct=fitstr, CSTRbasis=Cb, lambda=lam),
              data=list(datstr=datstruct, fitstr=fitStrH1234,
                             Cb=CSTRbasis, lam=lambda),
    start=as.list(parvec1234),
    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE) )
#70.23466 :  0.4 0.8 1.7 0.5
#  It.   1, fac=           1, eval (no.,total): ( 1,  1): new dev = 10.9019
#10.90186 :  0.4672504 0.8296264 1.9357045 0.4600684
#  It.   2, fac=           1, eval (no.,total): ( 1,  2): new dev = 10.6765
#10.67654 :  0.4666348 0.8408360 1.8546879 0.4679637
#  It.   3, fac=           1, eval (no.,total): ( 1,  3): new dev = 10.6739
#10.67385 :  0.4661300 0.8397689 1.7155038 0.4958567
#  It.   4, fac=           1, eval (no.,total): ( 1,  4): new dev = 10.6663
#10.66634 :  0.4661415 0.8397918 1.7182263 0.4963761
#  It.   5, fac=           1, eval (no.,total): ( 1,  5): new dev = 10.6666
#...
#  It.   5, fac= 0.000976562, eval (no.,total): (11, 15): new dev = 10.6663
#> warnings()
#Warning messages:
#1: Stepsize reduced below the minimum with parvec = 0.467, 0.83, 1.94, 0.46 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
# ...
#13: Stepsize reduced below the minimum with parvec = 0.466, 0.84, 1.72, 0.496 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#14: step factor 0.000488281 reduced below 'minFactor' of 0.000976562 in: nls(formula = ~CSTRres(kref = kref, EoverR = EoverR, a = a, b = b,   ...

et.nlsFit/60
# 16.23 minutes

nlsFit1234a <- CSTRfn(parvec=coef(nlsFit1234o), datstruct, fitStrH1234,
                     CSTRbasis, lambda)

save(file="CSTR-section12-nls-soln.Rdata",
     list=c("nlsFit1234o", "nlsFit1234a") )
# load(file="CSTR-section12-nls-soln.Rdata")

#%  get the singular values of the jacobian
(s1234a <- svd(nlsFit1234a$Dres)$d)
#%  log 10 of condition number

log10(s1234a[1]/s1234a[4])
# 3.0173 as before

#%  get final solution values

#[res, Dres, fitstruct] = CSTRfn(parvec, datstruct, fitstruct, ...
#                                         CSTRbasis, lambda, 0);

#%  display the parameter estimates

#disp(['Initial   values: ', num2str(parvec0)])
#disp(['Estimated values: ', num2str(parvec)])
#disp(['True      values: ', num2str(parvectru)])
# Matlab answers
#Initial   values: 0.4       0.8         1.7        0.5
#Estimated values: 0.46617   0.83961     1.7014     0.50021
#True      values: 0.461     0.83301     1.678      0.5

rbind("Initial values"=parvec1234,
      "Estimated values"=coef(nlsFit1234o),
      "True valules"=parvectru)
#                      kref    EoverR        a         b
#Initial values   0.4000000 0.8000000 1.700000 0.5000000
#Estimated values 0.4661344 0.8397756 1.718018 0.4964178
#True valules     0.4610000 0.8330100 1.678000 0.5000000

#%  get the final coefficient values for the solutions

coef1234o <- nlsFit1234a$fitstruct$coef0
Ccoef1234o = coef1234o[1:nbasis]
Tcoef1234o = coef1234o[(nbasis+1):(2*nbasis)]

Cfd1234o = fd(Ccoef1234o, CSTRbasis);
Tfd1234o = fd(Tcoef1234o, CSTRbasis);

Cvec1234o = eval.fd(tobs, Cfd1234o);
Tvec1234o = eval.fd(tobs, Tfd1234o);

Cvectru = eval.fd(tobs, Cfdtru$fd);
Tvectru = eval.fd(tobs, Tfdtru$fd);

#%  display maximum absolute errors

100.*quantile(abs(Cvec1234o - Cvectru)/Cvectru)
#          0%          25%          50%          75%         100%
#0.0009704469 0.0535586085 0.1144641708 0.1820555419 1.6127594082
(Cerrpct1234o = 100.*median(abs(Cvec1234o - Cvectru)/Cvectru))
# 0.1145 vs. 0.1149 for Matlab

100.*quantile(abs(Tvec1234o - Tvectru)/Tvectru)
#          0%          25%          50%          75%         100%
#5.638475e-05 1.377171e-02 2.636143e-02 4.506094e-02 3.984905e-01
(Terrpct1234o = 100.*median(abs(Tvec1234o - Tvectru)/Tvectru))
# 0.02640 vs. 0.02636 for Matlab

#%  plot the estimated and true solutions

#figure(6)
#subplot(2,1,1)
op <- par(mfrow=c(2,1))
plot(tobs, Cvectru, type="l", col="red", xlab="", ylab="C(t)",
     main="Concentration (red = true, blue = estimate)",
     ylim=c(1.2, 1.8))
lines(tobs, Cvec1234o, col="blue")
points(tobs, Cobs, col="blue")
plot(tobs, Tvectru, type="l", col="red", xlab="", ylab="C(t)",
     main="Temperature", ylim=c(330, 360))
lines(tobs, Tvec1234o, col="blue")
par(op)

##
##%%  13.  Fit the hot solution with the previous parameter estimates
##  (NOTE:  This script earlier said "Fit ... with the cool parameter";
##   however, that version appeared to actually use the hot parameter
##   estimates.  That version of this script used a much smaller
##   number of variable names with many variables being redefined
##   in the course of the script.  This made it difficult to
##   understand.  Therefore, many such variables have been given
##   different names for different uses.  If R and Matlab do not agree
##   on something, they might be using different variables.)
##
#%  replace true by estimated parameter values

fitStrHC <- fitStrHot
parvec1234o <- coef(nlsFit1234o)
fitStrHC$kref <- parvec1234o[1]
fitStrHC$EoverR <- parvec1234o[2]
fitStrHC$a <- parvec1234o[3]
fitStrHC$b <- parvec1234o[4]

#yinit = [Cinit_hot, Tinit_hot];
# yinitHot from 5.4 above

#%  Matlab:  load hot input into fitstruct
#   R:  Retained from 5.4 above in fitStrHot
#fitstruct.Tcin = Tcin_hot;
all.equal(fitStrHot$Tcin, hotStepInput[, "Tcin"])
#TRUE
plot(fitStrHot$Tcin, type='l')

#%  solve  differential equation with true parameter values

#h = waitbar(0,'Simulating Actual Model Output...');
#[t, yest_hot] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
#                            fitstruct, 'all_hot_step', Tlim);
#close(h)

library(odesolve)
hotStepSolHC <- lsoda(y=yinitHot, times=tspan, func=CSTR2, parms=
  list(fitstruct=fitStrHC, condition='all.hot.step', Tlim=Tlim) )

#%  set up separate variables for concentration and temperature

C.HC <- hotStepSolHC[, "Conc"]
T.HC <- hotStepSolHC[, "Temp"]

#Cerrpct = 100.*median(abs(C_hotest - C_hot)./Cinit_hot)
#Terrpct = 100.*median(abs(T_hotest - T_hot)./Tinit_hot)

100.*quantile(abs(C.HC - C_hot)/Cinit_hot)
#       0%       25%       50%       75%      100%
# 0.000000  1.462733  1.827320  2.535828 11.598890
CerrpctHC = 100.*median(abs(C.HC-C_hot)/Cinit_hot)
# 1.832445  vs. 0.1089 for Matlab ... ?????

100.*quantile(abs(T.HC - T_hot)/Tinit_hot)
#        0%        25%        50%        75%       100%
#0.00000000 0.03326384 0.05696738 0.08813524 0.53601556
TerrpctHC = 100.*median(abs(T.HC-T_hot)/Tinit_hot)
# 0.05692276 vs. 0.060908 from Matlab

#% plot the estimated and true hot solutions

#figure(8)
#subplot(2,1,1)

op <- par(mfrow=c(2,1))
plot(tspan, C_hot, type="l", col="red", ylab="C(t)", xlab="",
     main="Concentration (red = true, blue = estimates)")
lines(tspan, C.HC, col="blue")

plot(tspan, T_hot, type="l", col="red", ylab="T(t)", xlab="",
     main="Temperature")
lines(tspan, T.HC, col="blue")

par(op)
# Matlab plot suggests it got a better fit for this.

# Methods for 'nls' objects:

#anova(nlsFit) # only defined to compare nls fits
coef(nlsFit0.0)
confintNlsFit <- confint(nlsFit0.0)
deviance(nlsFit0.0)
df.residual(nlsFit0.0)
#fitted(nlsFit1234o)
formula(nlsFit0.0)
logLik(nlsFit0.0)
#predict(nlsFit1234o)
profileNlsFit <- profile(nlsFit0.0)
#residuals(nlsFit1234o)
summary(nlsFit0.0)
var.kEab <- vcov(nlsFit0.0)
s.kEab <- sqrt(diag(var.kEab))
cor.kEab <- var.kEab/outer(s.kEab,s.kEab)
round(cor.kEab, 2)
#        kref EoverR     a     b
#kref    1.00   0.79  0.18 -0.11
#EoverR  0.79   1.00  0.18 -0.14
#a       0.18   0.18  1.00 -1.00
#b      -0.11  -0.14 -1.00  1.00

# weights(nlsFit1234o)

##
##%%  14.  Estimate kref, EoverR and a with only
##         one response variable observed
##

# 14.1.  with only temperature observed.

#fitstruct.fit = [1,0];

fitStrHTemp <- fitStrHot
fitStrHTemp$fit <- c(Conc=0, Temp=1)

#%  use coefficients from spline smoothing as initial
#%  values for coefficient vectors

#fitstruct.coef0 = [Ccoef; Tcoef];
all.equal(fitStrHTemp$coef0, data.frame(Conc=Ccoef, Temp=Tcoef))
# TRUE

#%  set up the true parameter values

#parvectru = [kreftru, EoverRtru, atru, btru];
#fitstruct.b = btru;
fitStrHTemp$b <- btru

#%  calculate the estimates

#tic;
#[parvec, resnorm, residual, exitflag] = ...
#    lsqnonlin(@CSTRfn, parvec0, [], [], optionsCSTRfn, ...
#              datstruct, fitstruct, CSTRbasis, lambda);
#toc

et.nlsFitT <- system.time(
nlsFitHTemp <- NLS(formula=~CSTRres(kref=kref, EoverR=EoverR, a=a, b=b,
   datstruct=datstr, fitstruct=fitstr, CSTRbasis=Cb, lambda=lam),
              data=list(datstr=datstruct, fitstr=fitStrHTemp,
                             Cb=CSTRbasis, lam=lambda),
    start=as.list(parvec1234),
    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE) )
#24.25368 :  0.4 0.8 1.7 0.5
#  It.   1, fac=           1, eval (no.,total): ( 1,  1): new dev = 5.1049
#5.104895 :  0.4669959 0.8442459 1.8177825 0.4746155
#  It.   2, fac=           1, eval (no.,total): ( 1,  2): new dev = 5.09656
#5.096559 :  0.4674112 0.8432472 1.7790100 0.4843426
#  It.   3, fac=           1, eval (no.,total): ( 1,  3): new dev = 5.09652
#5.09652 :  0.4673539 0.8431802 1.7783980 0.4845128
#  It.   4, fac=           1, eval (no.,total): ( 1,  4): new dev = 5.09652
#5.09652 :  0.4673512 0.8431759 1.7783556 0.4845199

et.optim <- system.time(
fitHTempOptim <- optim(par=parvec1234, CSTRsse, control=list(trace=1),
      hessian=TRUE, datstruct=datstruct, fitstruct=fitStrHTemp,
      CSTRbasis=CSTRbasis, lambda=lambda) )
#  Nelder-Mead direct search function minimizer
#function value for initial parameters = 24.253679
#  Scaled convergence tolerance is 3.61408e-07
#Stepsize computed as 0.170000
#BUILD              5 151.892559 24.253679
#LO-REDUCTION       7 75.404476 23.414726
#HI-REDUCTION       9 57.803851 11.658206
#LO-REDUCTION      11 46.027880 10.927218
#HI-REDUCTION      13 24.253679 10.927218
#HI-REDUCTION      15 23.414726 7.367107
# ...
#HI-REDUCTION     371 5.096522 5.096521
#HI-REDUCTION     373 5.096522 5.096520
#HI-REDUCTION     375 5.096521 5.096520
#LO-REDUCTION     377 5.096521 5.096520
#REFLECTION       379 5.096521 5.096520
#Exiting from Nelder Mead minimizer
#    381 function evaluations used
#Warning messages:
#1: Stepsize reduced below the minimum with parvec = 0.45, 0.824, 1.76, 0.454 on iteration 8 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,
#2: Stepsize reduced below the minimum with parvec = 0.476, 0.785, 1.8, 0.478 on iteration 7 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,
save(fitHTempOptim, file="fitHTempOptim.Rdata")

et.nlminb <- system.time(
fitHTemp.nlminb <- nlminb(start=parvec1234, CSTRsse, control=list(trace=1),
      datstruct=datstruct, fitstruct=fitStrHTemp,
      CSTRbasis=CSTRbasis, lambda=lambda) )
#  0:     24.253679: 0.400000 0.800000  1.70000 0.500000
#  1:     5.4717994: 0.430618 0.793938  1.69470 0.475947
#  ...
# 19:     5.0965195: 0.467351 0.843176  1.77835 0.484520

et.BFGS <- system.time(
fitHTempBFGS <- optim(par=parvec1234, CSTRsse, method="BFGS",
      control=list(trace=1), hessian=TRUE,
      datstruct=datstruct, fitstruct=fitStrHTemp,
      CSTRbasis=CSTRbasis, lambda=lambda) )
 et.BFGS <- system.time(
+ fitHTempBFGS <- optim(par=parvec1234, CSTRsse, method="BFGS",
+       control=list(trace=1), hessian=TRUE,
+       datstruct=datstruct, fitstruct=fitStrHTemp,
+       CSTRbasis=CSTRbasis, lambda=lambda) )
#initial  value 24.253679
#initial  value 24.253679
#iter  10 value 5.096521
#iter  10 value 5.096521
#final  value 5.096521
#converged
#There were 43 warnings (use warnings() to see them)
#1: Dres0 has rank 36 < dim(Dres0) = 2110, 450 in iteration 1.  Ignoring singular subspace.  First (rank+1) singular values = 2.62663415956034e+87, 7.9543348596211e+86, 1.34554413837184e+86, 6.88290776803602e+85, 2.37720043427008e+84, 9.03623954728286e+83, 4.85409854039128e+82, 3.68743192927867e+82, 8.96301017461952e+81, 1.35075068899481e+81, 7.7956101085923e+80, 3.11569652390776e+80, 3.14754356858988e+79, 1.27127026656233e+79, 9.66842764582914e+78, 5.4421729448804e+77, 3.30528955432854e+77, 1.59771662666156e+77, 2.29305235111279e+76, 2.20616992908930e+76, 8.56537521877822e+75, 6.1740807002733e+75, 3.83213890257276e+75, 1.17533452113321e+75, 9.08537244727635e+74, 3.12969169628295e+74, 1.22550633680182e+74, 1.21353024155092e+74, 9.03380429389963e+72, 7.7912156066903e+72, 1.67487850658530e+72, 1.45346805994433e+72, 9.26781751419885e+71, 6.71232207220219e+71, 5.89967951140504e+71, 5.87353759710936e+71, 4.19176007081319e+71 in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#2: 5 of 960 values of log(abs(betaCC)) exceed the max = 236.594237631128;  thresholding. in: CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
# ...
#14: 10 of 960 values of log(abs(betaCC)) exceed the max = 236.594237631128;  thresholding. in: CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
#15: Stepsize reduced below the minimum with parvec = 661, -130, -113, -518 on iteration 1 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#16: Dres0 has rank 33 < dim(Dres0) = 2110, 450 in iteration 2.  Ignoring singular subspace.  First (rank+1) singular values = 8.27266427524798e+89, 7.20259593841941e+88, 2.62631439599031e+87, 9.83165427623727e+86, 7.95336649362295e+86, 1.34538033692328e+86, 6.88206982505691e+85, 3.35672793688527e+84, 2.37691104083859e+84, 9.03514169141088e+83, 4.85351552538531e+82, 3.68742373674273e+82, 8.96191924684419e+81, 4.17659180448618e+81, 1.350749532483e+81, 7.79466466725767e+80, 3.11569547020656e+80, 3.14754241174991e+79, 1.27112328389124e+79, 9.66836539526564e+78, 5.44217115328192e+77, 3.30488731649415e+77, 1.59771612283330e+77, 9.14885501802162e+76, 2.2927725045991e+76, 2.20590126818276e+76, 8.56429187892855e+75, 6.17397745002856e+75, 3.83186707427973e+75, 1.18139890393627e+75, 9.08562728588156e+74, 3.67846861044427e+74, 3.22574925264745e+74, 1.66914902035077e+74 in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#17: Dres0 has rank 45 < dim(Dres0) = 2110, 450 in iteration 3.  Ignoring singular subspace.  First (rank+1) singular values = 1.36835942223115e+86, 1.16793342243662e+86, 4.13198219419799e+85, 1.54183817464470e+85, 7.03121864766749e+84, 3.56351528831994e+84, 1.86500736788785e+83, 1.20889765298304e+83, 4.80041869285615e+82, 4.61351031711264e+81, 3.47383544394199e+81, 2.41158066653168e+81, 4.45803536245087e+80, 1.29910091761474e+80, 3.97532231554287e+79, 2.93800215781369e+79, 7.6542651045543e+78, 3.7288220048164e+78, 3.56851636797138e+78, 2.39592116234976e+78, 1.67074908085537e+78, 7.749521743393e+77, 2.94723826149075e+77, 1.61759305551093e+77, 1.11375965504932e+76, 8.43735739864358e+75, 4.91811178497923e+75, 3.15659239434470e+75, 1.42090745842028e+75, 5.52610607689986e+74, 3.19579361315506e+74, 2.12289397293076e+74, 1.26454610737948e+74, 3.50143949903395e+73, 2.55830798820949e+73, 5.13491036399458e+72, 4.18537704764922e+72, 2.90131538849364e+72, 9.03014598329064e+71, 7.97007997672222e+71, 6. [... truncated] in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#18: 2 of 960 values of log(abs(betaCC)) exceed the max = 236.594237631128;  thresholding. in: CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
#19: 1 of 960 values of log(abs(betaCC)) exceed the max = 236.594237631128;  thresholding. in: CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
#20: Dres0 has rank 45 < dim(Dres0) = 2110, 450 in iteration 4.  Ignoring singular subspace.  First (rank+1) singular values = 1.36671900312776e+86, 1.03061722203722e+86, 4.12700023131503e+85, 1.36056435052334e+85, 7.0228244917128e+84, 3.5591887150713e+84, 1.64567308243172e+83, 1.20731021294271e+83, 4.79160375004388e+82, 4.07109993498621e+81, 3.0502210176552e+81, 2.40471461249827e+81, 4.45331276744087e+80, 1.13846736057371e+80, 3.94808725491752e+79, 2.57455964234459e+79, 6.72106411386494e+78, 3.29018107563867e+78, 3.13533305702428e+78, 2.0990890273364e+78, 1.46418252288113e+78, 6.79326817789413e+77, 2.57970720194838e+77, 1.41935713134249e+77, 1.11096677762739e+76, 7.44550113534982e+75, 4.31275336100792e+75, 2.78551184259285e+75, 1.25387092170516e+75, 4.81278395059157e+74, 2.82009663375723e+74, 2.11657664499561e+74, 1.11589880504595e+74, 3.06941290449734e+73, 2.25756047204379e+73, 5.12623420270095e+72, 3.69328260266361e+72, 2.54109284553373e+72, 8.01762882692504e+71, 7.02375291963679e+71, 5. [... truncated] in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#21: 1 of 960 values of log(abs(betaCC)) exceed the max = 236.594237631128;  thresholding. in: CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
#22: Dres0 has rank 43 < dim(Dres0) = 2110, 450 in iteration 5.  Ignoring singular subspace.  First (rank+1) singular values = 1.36436907141415e+86, 8.02546317670237e+85, 4.1198634367384e+85, 1.05947904419045e+85, 7.01079961692214e+84, 3.55299073738439e+84, 1.28146429619685e+83, 1.20503553035439e+83, 4.77895498848676e+82, 3.17021327683460e+81, 2.44754939781046e+81, 2.24969695411365e+81, 4.44658713433275e+80, 8.55908976646808e+79, 3.90878447697525e+79, 1.93477401709288e+79, 5.07504403676784e+78, 2.56195687236048e+78, 2.36548705241738e+78, 1.57718106272791e+78, 1.10059610190817e+78, 5.11289212680266e+77, 1.93623371521044e+77, 1.07206101839586e+77, 1.10471963156292e+76, 5.79792389520512e+75, 3.25727151009611e+75, 2.16911991377158e+75, 9.76403537637953e+74, 3.71760882366843e+74, 2.19601891217486e+74, 2.09733792618382e+74, 8.6893688500941e+73, 2.30886515183495e+73, 1.75792446993994e+73, 5.09767297839446e+72, 2.87599937697087e+72, 1.91358013231849e+72, 6.21549015685329e+71, 5.47313489431808e+71,  [... truncated] in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#23: 1 of 960 values of log(abs(betaCC)) exceed the max = 236.594237631128;  thresholding. in: CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
#24: 1 of 960 values of log(abs(betaCC)) exceed the max = 236.594237631128;  thresholding. in: CSTRfitLS(coef1, datstruct, fitstruct, lambda, 1)
#25: Dres0 has rank 45 < dim(Dres0) = 2110, 450 in iteration 6.  Ignoring singular subspace.  First (rank+1) singular values = 1.36365613834975e+86, 7.08202056159313e+85, 4.11769829919760e+85, 9.34930756884983e+84, 7.00715135595220e+84, 3.55111047242630e+84, 4.53732391910521e+83, 1.20434558422072e+83, 1.13082023101592e+83, 4.77511930484497e+82, 2.79753486994937e+81, 2.42294608712413e+81, 1.98902859352469e+81, 4.44454593253216e+80, 7.49599749442163e+79, 3.89627860181026e+79, 1.69391799046411e+79, 4.45501581053059e+78, 2.26078318186119e+78, 2.07353969050831e+78, 1.38079233703420e+78, 9.6367770026885e+77, 4.48127806209146e+77, 1.69489937684374e+77, 9.4204707319405e+76, 1.10126299184138e+76, 5.11635288313078e+75, 2.86379880877699e+75, 1.91413018648104e+75, 8.61622920692052e+74, 3.28099211849122e+74, 2.08455772099957e+74, 1.93786377413672e+74, 7.66790733625262e+73, 2.02251134992934e+73, 1.55129670760345e+73, 5.07786853025171e+72, 2.53845951792551e+72, 1.67635887876641e+72, 5.48442292648315e+71,  [... truncated] in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#...
#29: Dres0 has rank 39 < dim(Dres0) = 2110, 450 in iteration 10.  Ignoring singular subspace.  First (rank+1) singular values = 2.46962460470095e+85, 7.40725922943645e+84, 1.29536977314064e+84, 1.27512414151221e+84, 6.32054069245095e+83, 1.71008178402774e+83, 2.05635674725325e+82, 1.75658295007180e+82, 8.90752641033155e+81, 2.06839722395062e+81, 3.80596071436606e+80, 6.6397254382806e+79, 5.11734083061055e+79, 5.74910302776551e+78, 8.80747425894892e+76, 4.13540230384686e+76, 3.61538340661968e+75, 5.3689845917628e+74, 1.08690227870789e+74, 9.36127135652053e+73, 3.50215525097826e+73, 2.23150954142171e+73, 1.57641391450363e+73, 5.84368024946491e+72, 4.66966428614869e+72, 3.54559342788578e+72, 1.40297458042960e+72, 4.1912916894818e+71, 2.83869096551342e+71, 1.42324562203236e+71, 7.88760835411115e+70, 4.64541907192365e+70, 3.35438445076084e+70, 2.25502531572203e+70, 1.01544687184171e+70, 9.3337612754605e+69, 6.98680853775722e+69, 6.62605507822225e+69, 5.64116275769335e+69, 4.01117979884277e+69 in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#30: Z'Z+R'R is not positive definite.  rank = 400;  dim = 450; increasing the 50smallest eigenvalues  to improve numeric stability. in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#31: Z'Z+R'R is ill conditioned.  Increasing the 436 smallest eigenvalues  to improve numeric stability. in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#32: Dres0 has rank 401 < dim(Dres0) = 2110, 450 in iteration 1.  Ignoring singular subspace.  First (rank+1) singular values = 99769355232869203978, 53529110408723652618, 28480258155647963136, 16323078033250699264, 13390094302987558912, 9306692922815381504, 5437135806167479296, 2439023121435212800, 1373249988599398656, 741278833621818112, 644138862371225728, 482838374869996992, 324463929132468608, 235069953251189056, 160174765077720928, 151647654981736512, 104882569139425632, 88574384869123712, 61707221012508336, 52251631614340744, 45305755063135216, 27282145025379336, 16456738040027608, 10612479669907044, 10494583111145418, 10201968007566054, 6583707701504801, 5011196943303142, 4093237621031927, 1621402650750643, 1000429154505716, 491664044104954, 451566825225660, 348869330227222, 271529847778930, 142293402355207, 125893232801337, 120407009324784, 103971780756566, 103485881862785, 91842961227864.3, 75277614582267.4, 61077675141905.5, 58791585175972.4, 54211454872604.2, 42241058828068.6, 4 [... truncated] in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#...
#41: Dres0 has rank 321 < dim(Dres0) = 2110, 450 in iteration 10.  Ignoring singular subspace.  First (rank+1) singular values = 7477824946949989376, 4001445968819990528, 2101415777382548736, 1219643884762974208, 974776163908151424, 656979357089595904, 383965583010148928, 177186728246312032, 106381689263043792, 58101768316896032, 49460724506784552, 38228880226336160, 25921216839516184, 18981625346241828, 13621708891123060, 11658234364469878, 9338585369129394, 6972301450747438, 5609476229732162, 3870017842949299, 2527921836226676, 1231220541529113, 954783320023108, 816593533807939, 489961409308487, 459193036900121, 100300846707915, 33119326494657.1, 22894498190903.1, 12015131332740.4, 2089072519322.06, 1683800801206.65, 1173167827797.79, 1115221151039.82, 884800687589.179, 872378569289.6, 816584909722.852, 623564165665.158, 586328328582.332, 482671372698.702, 399012631914.728, 325669594079.731, 320405113993.014, 248386618891.291, 238593817034.184, 233831182093.588, 204811983414.304, 18849940 [... truncated] in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#42: Z'Z+R'R is not positive definite.  rank = 431;  dim = 450; increasing the 19smallest eigenvalues  to improve numeric stability. in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#43: Z'Z+R'R is ill conditioned.  Increasing the 397 smallest eigenvalues  to improve numeric stability. in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#> et.BFGS /60
#    user   system  elapsed
#136.3972   3.0905 139.5402 minutes

et.CG <- system.time(
fitHTempCG <- optim(par=parvec1234, CSTRsse, method="CG",
      control=list(trace=1), hessian=TRUE,
      datstruct=datstruct, fitstruct=fitStrHTemp,
      CSTRbasis=CSTRbasis, lambda=lambda) )
#  Conjugate gradients function minimizer
#Method: Fletcher Reeves
#tolerance used in gradient test=7.27596e-12
#0 1 24.253679
#parameters    0.40000    0.80000    1.70000    0.50000
#****** i> 1 9 7.521463
#parameters    0.44228    0.79163    1.69268    0.46679
# i< 2 11 5.489809
# ...
# i> 99 360 5.099004
#parameters    0.46545    0.84147    1.70066    0.49957
# i> 100 362 5.098996
#parameters    0.46547    0.84146    1.70079    0.49955
#There were 43 warnings (use warnings() to see them)

#str(fitHTemp.nlminb)
#fitHTemp.nlminb$par-parvec0
#names(fitHTemp.nlminb$par) <- names(parvec0)

#et.nlsFitHTemp2 <- system.time(
#nlsFitHTemp2 <- nls(formula=~CSTRres(kref=kref, EoverR=EoverR, a=a, b=b,
#   datstruct=datstr, fitstruct=fitstr, CSTRbasis=Cb, lambda=lam),
#              data=list(datstr=datstruct, fitstr=fitStrHTemp,
#                             Cb=CSTRbasis, lam=lambda),
#    start=as.list(fitHTemp.nlminb$par),
#    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE) )

et.SANN <- system.time(
fitHTempSANN <- optim(par=parvec1234, CSTRsse, method="SANN",
      control=list(trace=1), hessian=TRUE,
      datstruct=datstruct, fitstruct=fitStrHTemp,
      CSTRbasis=CSTRbasis, lambda=lambda) )
)
#sann objective function values
#initial       value 24.253679
#iter     1000 value 5.337875
#iter     2000 value 5.337875
#Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
#In addition: There were 50 or more warnings (use warnings() to see the first 50)
#iter     3000 value 5.337875
#iter     4000 value 5.175040
#Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
#In addition: Warning messages:
#1: svd failed using LINPACK = FALSE with n = 2110and p = 450;  x stored in '.svd.LINPACK.error.matrix' in: svd2(Dres0)
#2: Stepsize reduced below the minimum with parvec = 0.382, 0.643, 2.42, 0.303 on iteration 9 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,
#3: Stepsize reduced below the minimum with parvec = 0.296, 0.515, 1.95, 0.287 on iteration 8 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,
#4: Stepsize reduced below the minimum with parvec = 0.526, 0.786, 1.89, 0.508 on iteration 9 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,
#...
#Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
#In addition: Warning message:
#svd failed using LINPACK = FALSE with n = 2110and p = 450;  x stored in '.svd.LINPACK.error.matrix' in: svd2(Dres0)
#iter     9999 value 5.175040
#final         value 5.175040
#sann stopped after 9999 iterations
#Warning message:
#svd failed using LINPACK = FALSE with n = 2110and p = 450;  x stored in '.svd.LINPACK.error.matrix' in: svd2(Dres0)
#
#outer(et.SANN, c(1, 60, 60^2, 24*60^2), "/")
#            seconds or minutes  or hours  or days
#user.self  376365.11 6272.7518 104.545864 4.3560777
#sys.self     7895.69  131.5948   2.193247 0.0913853
#elapsed    385924.48 6432.0747 107.201244 4.4667185

#%  get final solution values

#[res, Dres, fitstruct] = CSTRfn(parvec, datstruct, fitstruct, ...
#                                         CSTRbasis, lambda, 0);

#%  display the parameter estimates

#disp(['Initial   values: ', num2str(parvec0)])


#disp(['Estimated values: ', num2str(parvec)])
#disp(['True      values: ', num2str(parvectru)])
# Matlab:
#Initial   values: 0.4         0.8         1.7         0.5
#Estimated values: 0.46327     0.83906      1.7033     0.49754
#True      values: 0.461     0.83301       1.678         0.5



#%  get the final coefficient values for the solutions

#coef     = fitstruct.coef0;
#Ccoefest = coef(1:nbasis);
#Tcoefest = coef(nbasis+1:2*nbasis);

#Cfdest = fd(Ccoefest, CSTRbasis);
#Tfdest = fd(Tcoefest, CSTRbasis);

#Cvecest = eval_fd(tobs, Cfdest);
#Tvecest = eval_fd(tobs, Tfdest);

#Cvectru = eval_fd(tobs, Cfdtru);
#Tvectru = eval_fd(tobs, Tfdtru);

#%  display maximum absolute errors

#Cerrpct = 100.*median(abs(Cvecest - Cvectru)./Cvectru)
#Terrpct = 100.*median(abs(Tvecest - Tvectru)./Tvectru)

#%  plot the estimated and true solutions

# Matlab answers for figure(9)
# fit9Mat <- read.csv('fit9Mat.csv')
# to compare Matlab with R ...
# ... when and if we get answers for fitHTemp ...

#figure(9)
#subplot(2,1,1)
#plot(tobs, Cvectru, 'r-', tobs, Cvecest, 'b-')
#ylabel('\fontsize{16} C(t)')
#title('\fontsize{16} Concentration (red = true, blue = estimated)')
#axis([0, Tlim, 1.2, 1.8])
#subplot(2,1,2)
#plot(tobs, Tvectru, 'r-', tobs, Tvecest, 'b-')
#ylabel('\fontsize{16} T(t)')
#title('\fontsize{16} Temperature')
#axis([0, Tlim, 330, 360])


# 14.2.  with only concentration observed.

#fitstruct.fit = [1,0];

fitStrHConc <- fitStrHot
fitStrHConc$fit <- c(Conc=1, Temp=0)

#%  use coefficients from spline smoothing as initial
#%  values for coefficient vectors

#fitstruct.coef0 = [Ccoef; Tcoef];
all.equal(fitStrHConc$coef0, data.frame(Conc=Ccoef, Temp=Tcoef))
# TRUE

#%  set up the true parameter values

#parvectru = [kreftru, EoverRtru, atru, btru];
#fitstruct.b = btru;
fitStrHConc$b <- btru

#%  calculate the estimates

#tic;
#[parvec, resnorm, residual, exitflag] = ...
#    lsqnonlin(@CSTRfn, parvec0, [], [], optionsCSTRfn, ...
#              datstruct, fitstruct, CSTRbasis, lambda);
#toc

et.Conc <- system.time(
nlsFitHConc <- NLS(formula=~CSTRres(kref=kref, EoverR=EoverR, a=a, b=b,
   datstruct=datstr, fitstruct=fitstr, CSTRbasis=Cb, lambda=lam),
              data=list(datstr=datstruct, fitstr=fitStrHConc,
                             Cb=CSTRbasis, lam=lambda),
    start=as.list(parvec1234),
    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE) )
#5.275455 :  0.46100 0.83301 1.67800 0.50000
#  It.   1, fac=           1, eval (no.,total): ( 1,  1): new dev = 5.25505
#5.255047 :  0.4629298 0.8368685 1.6220946 0.5127860
#  It.   2, fac=           1, eval (no.,total): ( 1,  2): new dev = 5.25488
#5.254876 :  0.4629855 0.8369806 1.6236225 0.5126904
#  It.   3, fac=           1, eval (no.,total): ( 1,  3): new dev = 5.25488
#5.254876 :  0.4629772 0.8369686 1.6235186 0.5127036

#%  get final solution values

#[res, Dres, fitstruct] = CSTRfn(parvec, datstruct, fitstruct, ...
#                                         CSTRbasis, lambda, 0);
parvecHConc <- coef(nlsFitHConc)
nlsFitHConc. <- CSTRfn(parvec=parvecHConc, datstruct,
                      fitStrHConc, CSTRbasis, lambda)

#%  display the parameter estimates

#disp(['Initial   values: ', num2str(parvec0)])
#disp(['Estimated values: ', num2str(parvec)])
#disp(['True      values: ', num2str(parvectru)])
rbind("Initial values"=parvec1234,
      "Estimated values"=parvecHConc,
      "True valules"=parvectru)

#%  get the final coefficient values for the solutions

#coef     = fitstruct.coef0;
coefHConc <- nlsFitHConc.$fitstruct$coef0
#Ccoefest = coef(1:nbasis);
CcoefestHConc <- coefHConc[1:nbasis]
#Tcoefest = coef(nbasis+1:2*nbasis);
TcoefestHConc <- coefHConc[(nbasis+1):(2*nbasis)]

#Cfdest = fd(Ccoefest, CSTRbasis);
CfdestHConc <- fd(CcoefestHConc, CSTRbasis)
#Tfdest = fd(Tcoefest, CSTRbasis);
TfdestHConc <- fd(TcoefestHConc, CSTRbasis)

#Cvecest = eval_fd(tobs, Cfdest);
CvecestHConc <- eval.fd(tobs, CfdestHConc)
#Tvecest = eval_fd(tobs, Tfdest);
TvecestHConc <- eval.fd(tobs, TfdestHConc)

#Cvectru = eval_fd(tobs, Cfdtru)
#Tvectru = eval_fd(tobs, Tfdtru);
# Conputed above following 'nlsFit' computation

#%  display maximum absolute errors

quantile((CvecestHConc-Cvectru)/Cvectru)
#           0%           25%           50%           75%          100%
#-0.0109490852 -0.0014348313 -0.0006694950  0.0003116575  0.0029172501
#Cerrpct = 100.*median(abs(Cvecest - Cvectru)./Cvectru)
CerrpctHConc <- 100*median(abs((CvecestHConc-Cvectru)/Cvectru))
# 0.0934

quantile((TvecestHConc-Tvectru)/Tvectru)
#           0%           25%           50%           75%          100%
#-4.597968e-04 -9.886245e-05  1.483907e-07  1.618749e-04  1.755848e-03
#Terrpct = 100.*median(abs(Tvecest - Tvectru)./Tvectru)
TerrpctHConc <- 100*median(abs((TvecestHConc-Tvectru)/Tvectru))
# 0.0124

#%  plot the estimated and true solutions

#figure(10)% was figure(9): second copy for obs. Conc only
#subplot(2,1,1)
op <- par(mfrow=c(2,1))
plot(tobs, Cvectru, type="l", col="red", ylab="C(t)", xlab="",
     main="Concentration (red = true, blue = estimates)")
lines(tobs, CvecestHConc, col="blue")

plot(tobs, Tvectru, type="l", col="red", ylab="T(t)", xlab="",
     main="Temperature")
lines(tobs, TvecestHConc, col="blue")

par(op)

save(file="CSTR_fig10.Rdata")

##
##%%  15.  Fit the hot solution with the cool parameter estimates
##
#%  replace true by estimated parameter values

fitStrHTC <- fitStrHTemp
fitStrHTC$kref <- parvecHConc[1]
fitStrHTC$EoverR <- parvecHConc[2]
fitStrHTC$a <- parvecHConc[3]
fitStrHTC$b <- parvecHConc[4]

#yinit = [Cinit_hot, Tinit_hot];

#%  load hot input into fitstruct

#fitstruct.Tcin = Tcin_hot;

#%  solve  differential equation with true parameter values

#h = waitbar(0,'Simulating Actual Model Output...');
#[t, yest_hot] = ode45(@CSTR2, tspan, yinitHot, odeoptions, ...
#                            fitstruct, 'all_hot_step', Tlim);
#close(h)

library(odesolve)

hotStepSolHTC <- lsoda(y=yinitHot, times=tspan, func=CSTR2, parms=
  list(fitstruct=fitStrHTC, condition='all.hot.step', Tlim=Tlim) )
#%  set up separate variables for concentration and temperature

C.HTC <- hotStepSolHTC[, "Conc"]
T.HTC <- hotStepSolHTC[, "Temp"]

#% plot the estimated and true hot solutions

#figure(11) % was figure(10)
#subplot(2,1,1)
op <- par(mfrow=c(2,1))
plot(tspan, C_hot, type="l", col="red", ylab="C(t)", xlab="",
     main="Concentration (red = true, blue = estimates)")
lines(tspan, C.HTC, col="blue")

plot(tspan, T_hot, type="l", col="red", ylab="T(t)", xlab="",
     main="Temperature")
lines(tspan, T.HTC, col="blue")

par(op)

save(file="CSTR_fig11.Rdata")
