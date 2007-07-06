##
##%%  1.  Introduction
##
#%  Estimate the four parameters defining two differential equations
#%  that describe the behavior of the output concentration
#%  and temperature for a nonisotherman continuously stirred 
#%  tank reactor, abbreviated nonisothermal CSTR.  
#%  The system of two nonlinear equations has five forcing or  
#  input functions.
#  These equations are taken from
#%  Marlin, T. E. (2000) Process Control, 2nd Edition, McGraw Hill,
#%  pages 899-902.

#%  Last modified 16 April 2007 by Spencer Graves
#%  Matlab version previously modified 26 October 2005

#%addpath('c:/matlab7/fdaM')
#addpath('d:/spencerg/statmtds/fda/Matlab/fdaM')
library(fda) 
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

#%[F,CA0,T0,Tcin_cool,Fc] = CSTR2in(tspan, 'all_cool_step');
#[F,CA0,T0,Tcin_cool] = CSTR2in(tspan, 'all_cool_step');
#[F,CA0,T0,Tcin_hot, Fc] = CSTR2in(tspan, 'all_hot_step');
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
#yinit = c(Conc = Cinit_cool, Temp=Tinit_cool)
yinit = c(Conc = Cinit_cool, Temp=Tinit_cool)

#%%  load cool input into fitstruct

fitStrCool <- fitstruct
fitStrCool$Tcin = coolStepInput[, "Tcin"]

#%  5.2.  solve  differential equation with true parameter values

#h = waitbar(0,'Simulating Cool Output...');
#[t, y_cool] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
#                            fitstruct, 'all_cool_step', Tlim);
#close(h)

coolStepSoln <- lsoda(y=yinit, times=tspan, func=CSTR2, parms=
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
library(odesolve)
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

#Cfdtru = smooth_basis(t, C_cool, fdParobj);
#Tfdtru = smooth_basis(t, T_cool, fdParobj);
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

#plotfit_fd(C_cool, t, Cfdtru)
plotfit.fd(C_cool, tspan, Cfdtru$fd, titles="Concentration",
           xlab="", sub="", ylab="C(t)")

#xlabel('')
#ylabel('\fontsize{16} C(t)')
#title('\fontsize{16} Concentration')
#axis([0, Tlim, 1.2, 1.8])
#subplot(2,1,2)

#plotfit_fd(T_cool, t, Tfdtru)
plotfit.fd(T_cool, tspan, Tfdtru$fd, titles="Temperature",
           xlab="", sub="", ylab="T(t)")

#xlabel('')
#ylabel('\fontsize{16} T(t)')
#title('\fontsize{16} Temperature')
#axis([0, Tlim, 330, 360])

par(op)

#%  7.3.  compute weights to correct for differences in
#%  scale for the two variables

#Cvectru = eval_fd(t, Cfdtru);
#str(Cfdtru)
CvectruSpan = eval.fd(tspan, Cfdtru$fd);

Cwt     = var(CvectruSpan);

#Tvectru = eval_fd(t, Tfdtru);
TvectruSpan = eval.fd(tspan, Tfdtru$fd);

Twt     = var(TvectruSpan);

#wt   = [Cwt, Twt];
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

#Cvec = eval_fd(tobs, Cfdtru);
Cvec = eval.fd(tobs, Cfdtru$fd);
#Tvec = eval_fd(tobs, Tfdtru);
Tvec = eval.fd(tobs, Tfdtru$fd);

#%  8.3.  add errors

#Cobs    = Cvec + randn(N,1).*stderrC;
Cobs    = Cvec + rnorm(N)*stderrC;
#Tobs    = Tvec + randn(N,1).*stderrT;
Tobs    = Tvec + rnorm(N)*stderrT;
#Cobs.R <- Cobs
#Tobs.R <- Tobs
# rm(list=c("Cobs","Tobs"))

##
##**** IF YOU WANT TO COMPARE ANSWERS WITH MATLAB,
## TRANSFER MATLAB SIMULATIONS TO ENSURE COMPARABILITY
##
## You can still get good answers from the R code,
## but you won't be able to compare them with Matlab.
##
#library(R.matlab)
#confirm <- readMat("CSTRsim.mat")
#Cobs <- confirm$Cobs
#Tobs <- confirm$Tobs

#%  8.4.  smooth the noisy data using smoothing splines

lambdaC = 1e0;
CfdPar  = fdPar(CSTRbasis, 2, lambdaC);
#Cfdsmth = smooth_basis(tobs, Cobs, CfdPar);
Cfdsmth = smooth.basis(tobs, Cobs, CfdPar);
#Cfdsmth.old = smooth.basis(tobs, Cobs, CfdPar);
#all.equal(Cfdsmth, Cfdsmth.old)


lambdaT = 1e-1;
TfdPar  = fdPar(CSTRbasis, 2, lambdaT);
#Tfdsmth = smooth_basis(tobs, Tobs, TfdPar);
Tfdsmth = smooth.basis(tobs, Tobs, TfdPar);

##
##***** save.image("CSTR-section8-Matlab sims.Rdata")
##

# TRANSFER MATLAB smooth_basis answers TO check 
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

#plot(tobs, Cobs, '.');
plot(tobs, Cobs, ylab="C(t)", xlab='',
     main="Concentration");

#xlabel('')
#ylabel('\fontsize{16} C(t)')
#title('\fontsize{16} Concentration')
#axis([0, Tlim, 1.2, 1.8])
#subplot(2,1,2)
#plot(tobs, Tobs, '.');
plot(tobs, Tobs, xlab='', ylab="T(t)",
     main="Temperature");
#xlabel('')
#ylabel('\fontsize{16} T(t)')
#title('\fontsize{16} Temperature')
#axis([0, Tlim, 330, 360])
par(op)

#%  8.6.  plot the data and the smooths

#figure(5)
#subplot(2,1,1)

op <- par(mfrow=c(2,1))
#plotfit_fd(Cobs, tobs, Cfdsmth);
plotfit.fd(Cobs, tobs, Cfdsmth$fd, xlab='', ylab='C(t)',
           titles="Concentration", sub='');
#xlabel('')
#ylabel('\fontsize{16} C(t)')
#title('\fontsize{16} Concentration')
#axis([0, Tlim, 1.2, 1.8])
#subplot(2,1,2)

#plotfit_fd(Tobs, tobs, Tfdsmth);
#plotfit_fd(Tobs, tobs, Tfdsmth);
plotfit.fd(Tobs, tobs, Tfdsmth$fd, xlab='', ylab='T(t)',
           titles="Temperature", sub='');
#xlabel('')
#ylabel('\fontsize{16} T(t)')
#title('\fontsize{16} Temperature')
#axis([0, Tlim, 330, 360])

par(op)

##
##%%  9.  Load data into struct variable datstruct
##

#%  this is information that only depends on the
#%  data and sampling design

#%  9.1.  weights for variables

#datstruct.Cwt = Cwt;
#datstruct.Twt = Twt;
datstruct0 <- list(Cwt=Cwt, Twt=Twt)

#%  9.2.  data

#datstruct0$y = data.frame(Conc=Cobs, Temp=Tobs);
datstruct0$y = cbind(Conc=Cobs, Temp=Tobs);

#%  9.3.  basis values at sampling points

#basismat            = eval_basis(tobs, CSTRbasis);
basismat            = eval.basis(tobs, CSTRbasis);
#Dbasismat           = eval_basis(tobs, CSTRbasis, 1);
Dbasismat           = eval.basis(tobs, CSTRbasis, 1);
datstruct0$basismat  = basismat;
datstruct0$Dbasismat = Dbasismat;

#%  9.4.  forcing function values at quadrature points

#[Fquad, CA0quad, T0quad, Tcinquad, Fcquad] = ...
#                     CSTR2in(quadpts, 'all_cool_step');
#str(CSTRbasis)
quadpts <- CSTRbasis$quadvals[, "quadpts"]
coolStepQuadInput <- CSTR2in(quadpts, 'all.cool.step');

#datstruct.F    = Fquad;
#datstruct.CA0  = CA0quad;
#datstruct.T0   = T0quad;
#datstruct.Tcin = Tcinquad;
#datstruct.Fc   = Fcquad;

datstruct <- c(datstruct0, as.data.frame(coolStepQuadInput))

#%  9.5.  basis values at quadrature points

#quadbasismat            = eval_basis(quadpts, CSTRbasis);
quadbasismat            = eval.basis(quadpts, CSTRbasis);
#Dquadbasismat           = eval_basis(quadpts, CSTRbasis, 1);
Dquadbasismat           = eval.basis(quadpts, CSTRbasis, 1);

#datstruct.quadpts       = quadpts;
#datstruct.quadwts       = quadwts;
#datstruct.quadbasismat  = quadbasismat;
#datstruct.Dquadbasismat = Dquadbasismat;

datstruct$quadpts       = CSTRbasis$quadvals[, "quadpts"]
datstruct$quadwts       = CSTRbasis$quadvals[, "quadwts"]
datstruct$quadbasismat  = quadbasismat;
datstruct$Dquadbasismat = Dquadbasismat;

##
##  10.  Analyze the noisy data to estimate
##       * a CSTRbasis representation (section 11) and
##       * parameters of the CSTR differential equations (sect. 12) 
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

#fitstruct.fit = [1,1];
fitStrHot$fit <- c(Conc=1,Temp=1)

#%  11.2.  Use coefficients from spline smoothing as initial
#%  values for coefficient vectors

#Ccoef = getcoef(Cfdsmth);
#Tcoef = getcoef(Tfdsmth);
#fitStrHot.coef0 = [Ccoef; Tcoef];

#Ccoef = getcoef(Cfdsmth);
#Tcoef = getcoef(Tfdsmth);

Ccoef <- Cfdsmth$coef
Tcoef <- Tfdsmth$coef

#fitstruct.coef0 = [Ccoef; Tcoef];
fitStrHot$coef0 <- data.frame(Conc=Ccoef, Temp=Tcoef)

#%  11.3.  Specify that all parameters will be estimated

#fitstruct.estimate = [1, 1, 1, 1];
fitStrHot$estimate <- rep(1,4)

#estind = find(fitstruct.estimate == 1);
estind <- which(fitStrHot$estimate==1)
  
#%  11.4.  set up the true parameter values

#parvec0 = [kreftru, EoverRtru, atru, btru];
parvec0 <- c(kref=kreftru, EoverR=EoverRtru,
             a=atru, b=btru);

#%  11.5.  specify smoothing parameter values

#lambda = 1e1.*[lambdaC, lambdaT];
lambda = 1e1*c(lambdaC, lambdaT);

#%  11.6.  Obtain coefficients for a solution in CSTRbasis 

#[res, jacobian] = CSTRfn(parvec0, datstruct, fitstruct, ...
#                         CSTRbasis, lambda);
start.time <- proc.time()
CSTR.1 <- CSTRfn(parvec0, datstruct, fitStrHot, CSTRbasis, lambda)
(et.1 <- proc.time()-start.time)

#library(R.matlab)
#CSTR.1Mat <- readMat('CSTR1.mat')
#str(CSTR.1Mat)
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
start.time <- proc.time()
CSTR.2 <- CSTRfn(parvec12, datstruct, fitStrH12, CSTRbasis, lambda)
(et.2 <- proc.time()-start.time)

str(CSTR.1)
str(CSTR.2)
all.equal(CSTR.1$res, CSTR.2$res)
# TRUE 
all.equal(CSTR.1$Dres[, 1:2], CSTR.2$Dres)
# TRUE  

#CSTR.2Mat <- readMat('CSTR12.mat')
#str(CSTR.2Mat)
#sqrt(mean(CSTR.2$res^2)) # 0.168
# d.res.2 <- (CSTR.2$res - as.vector(CSTR.2Mat$res))
#sqrt(mean(d.res.2^2)) # 0.021





# was :  5.46e-7
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

#s = svd(full(jacobian));
s12 = svd(CSTR.2$Dres)$d

#%  log 10 of condition number
log10(s12[1]/s12[2])
# R = 0.826
#%  same as Matlab 

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
start.time <- proc.time()
CSTR.13 <- CSTRfn(parvec13, datstruct, fitStrH13, CSTRbasis, lambda)
(elapsed.time.13 <- proc.time()-start.time)
                                     
#%  get the singular values of the jacobian

#s = svd(full(jacobian));
s13 = svd(CSTR.13$Dres)$d

#%  log 10 of condition number
log10(s13[1]/s13[2])
#  1.1789 on 2007.05.30
# Matches Matlab 

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
start.time <- proc.time()
CSTR.34 <- CSTRfn(parvec34, datstruct, fitStrH34, CSTRbasis, lambda)
(et.34 <- proc.time()-start.time)

#%  get the singular values of the jacobian

s34 = svd(CSTR.34$Dres)$d

#%  log 10 of condition number
log10(s34[1]/s34[2])
# 2.1200 on 2007.05.30 
#%  matches Matlab 

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
start.time <- proc.time()
CSTR.123 <- CSTRfn(parvec123, datstruct, fitStrH123, CSTRbasis, lambda)
(et.123 <- proc.time()-start.time)

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
#optionsCSTRfn = optimset('LargeScale', 'iter', 'MaxIter', 50, 'TolCon',
#     tolval, 'TolFun', tolval, 'TolX',   tolval, 'TolPCG', tolval, 
#                          'Jacobian', 'on')

#%  use coefficients from spline smoothing as initial
#%  values for coefficient vectors

#Ccoef = getcoef(Cfdsmth);
#Tcoef = getcoef(Tfdsmth);

#%fitstruct.coef0 = [Ccoef; Tcoef];
fitStrH1234 = fitStrHot ; 
fitStrH1234$coef0 = data.frame(Conc=Ccoef, Temp=Tcoef)

#%  Use only the first, second and third parameters

fitStrH1234$fit = c(1,1)

fitStrH1234$estimate = c(1, 1, 1, 1)  
estind1234 = which(fitStrH1234$estimate == 1)

#%  set up some initial values for parameters

parvec1234 = c(0.4, 0.8, 1.7, 0.5)

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

.datstruct <- datstruct
.fitstruct <- fitStrH1234
.CSTRbasis <- CSTRbasis
.lambda <- lambda

save(file="CSTR-section12-nls-prep.Rdata",
     list=c("parvec0", "datstruct", "fitStrCool", "fitStrHot", 
		"fitStrH1234", "CSTRbasis", "lambda"))
#load("CSTR-section12-nls-prep.Rdata")

start.time <- proc.time()
nlsFit0.0 <- nls(formula=
     ~CSTRres0(kref=kref, EoverR=EoverR, a=a, b=b,gradwrd=FALSE),
    start=as.list(parvec0),
    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE)
(et.nls0.0 <- proc.time()-start.time)
# warnOnly=TRUE to force output:
# CSTR computations are so intense that 'nls' overestimates
#   the numerical precision feasible and often
#   terminates with an error & returns nothing
#   when it can't achieve that level of precision.

save(file="CSTR-nlsFit0.0.Rdata", list="nlsFit0.0")

deviance(nlsFit0.0)
#  10.66635
# vs. 10.6666 reported by Matlab

coef(nlsFit0.0)
#     kref    EoverR         a         b 
#0.4661271 0.8397104 1.7185853 0.4963070 
# vs.
#0.4662    0.8396    1.7014    0.5002 per Matlab

# nls does not return the Jacobian, so call CSTFfn directly

nlsFit0.0a <- CSTRfn(parvec=coef(nlsFit0.0), datstruct, fitStrH1234,
                     CSTRbasis, lambda)

#s <- svd(jacobian)
(s1234 <- svd(nlsFit0.0a$Dres)$d)
# log10 of the condition number 
log10(s1234[1]/s1234[4])
# 3.017316 vs. 3.0146 for the Matlab solution 

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
start.time <- proc.time()
nlsFit <- nls(formula=~CSTRres(kref=kref, EoverR=EoverR, a=a, b=b,
   datstruct=datstr, fitstruct=fitstr, CSTRbasis=Cb, lambda=lam),
              data=list(datstr=datstruct, fitstr=fitStrH1234,
                             Cb=CSTRbasis, lam=lambda), 
    start=as.list(parvec0),
    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE)
(et.nls <- proc.time()-start.time)

#10.71492 :  0.46100 0.83301 1.67800 0.50000 
#  It.   1, fac=           1, eval (no.,total): ( 1,  1): new dev = 10.6664
#10.66639 :  0.4662509 0.8398304 1.7211810 0.4959086 
#  It.   2, fac=           1, eval (no.,total): ( 1,  2): new dev = 10.6666
#  It.   2, fac=         0.5, eval (no.,total): ( 2,  3): new dev = 10.6664
#10.66639 :  0.4660565 0.8395190 1.7103179 0.4980595 
#  It.   3, fac=           1, eval (no.,total): ( 1,  4): new dev = 10.6667
#  It.   3, fac=         0.5, eval (no.,total): ( 2,  5): new dev = 10.6664
#  It.   3, fac=        0.25, eval (no.,total): ( 3,  6): new dev = 10.6664
#10.66636 :  0.4660900 0.8395955 1.7187953 0.4962464 
#  It.   4, fac=         0.5, eval (no.,total): ( 1,  7): new dev = 10.6663
#10.66635 :  0.4661271 0.8397103 1.7185536 0.4963137 
#  It.   5, fac=           1, eval (no.,total): ( 1,  8): new dev = 10.6669
#...
#  It.   5, fac= 0.000976562, eval (no.,total): (11, 18): new dev = 10.6663
#There were 30 warnings (use warnings() to see them)
# 1: Stepsize reduced below the minimum with parvec = 0.461, 0.833, 1.68, 0.5 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#...
#29: Stepsize reduced below the minimum with parvec = 0.466, 0.84, 1.72, 0.496 in trying to optimize 450 coefficients;  using suboptimal coefficients;  saved in '..CSTRfn.coef1.gen.5' in: CSTRfn(parvec = pv, datstruct = datstruct, fitstruct = fitstruct,   ...
#30: step factor 0.000488281 reduced below 'minFactor' of 0.000976562 in: nls(formula = ~CSTRres(kref = kref, EoverR = EoverR, a = a,   ...

nlsFit.a <- CSTRfn(parvec=coef(nlsFit), datstruct, fitStrH1234,
                     CSTRbasis, lambda)

save(file="CSTR-section12-nls-soln.Rdata",
     list=c("nlsFit", "nlsFit.a") )

#%  get the singular values of the jacobian
#s = svd(full(jacobian));

(s1234a <- svd(nlsFit.a$Dres)$d)
#%  log 10 of condition number

log10(s1234a[1]/s1234a[4])
# 3.017316 as before 

#%  get final solution values

#[res, Dres, fitstruct] = CSTRfn(parvec, datstruct, fitstruct, ...
#                                         CSTRbasis, lambda, 0);

#%  display the parameter estimates

#disp(['Initial   values: ', num2str(parvec0)])
#disp(['Estimated values: ', num2str(parvec)])
#disp(['True      values: ', num2str(parvectru)])
rbind("Initial values"=parvec0,
      "Estimated values"=coef(nlsFit),
      "True valules"=parvectru)

#%  get the final coefficient values for the solutions

str(nlsFit)
str(nlsFit.a)
#coef     = fitstruct.coef0;
coef. <- nlsFit.a$fitstruct$coef0

Ccoefest = coef.[1:nbasis]

Tcoefest = coef.[(nbasis+1):(2*nbasis)]

Cfdest = fd(Ccoefest, CSTRbasis);
Tfdest = fd(Tcoefest, CSTRbasis);

Cvecest = eval.fd(tobs, Cfdest);
Tvecest = eval.fd(tobs, Tfdest);

Cvectru = eval.fd(tobs, Cfdtru$fd);
Tvectru = eval.fd(tobs, Tfdtru$fd);

#%  display maximum absolute errors

100.*quantile(abs(Cvecest - Cvectru)/Cvectru)
#          0%          25%          50%          75%         100% 
#0.0004968317 0.0538577251 0.1150360258 0.1818828124 1.6100816866 
Cerrpct = 100.*median(abs(Cvecest - Cvectru)/Cvectru)

100.*quantile(abs(Tvecest - Tvectru)/Tvectru)
#          0%          25%          50%          75%         100% 
#0.0001489455 0.0137511402 0.0263061496 0.0450480997 0.3984394519
(Terrpct = 100.*median(abs(Tvecest - Tvectru)/Tvectru))


#%  plot the estimated and true solutions

#figure(6)
#subplot(2,1,1)
op <- par(mfrow=c(2,1))

#plot(tobs, Cvectru, 'r-', tobs, Cvecest, 'b-', tobs, Cobs, 'b.')
plot(tobs, Cvectru, type="l", col="red", xlab="", ylab="C(t)",
     main="Concentration (red = true, blue = estimate)",
     ylim=c(1.2, 1.8))
lines(tobs, Cvecest, col="blue")
points(tobs, Cobs, col="blue")
     
#ylabel('\fontsize{16} C(t)')
#title('\fontsize{16} Concentration (red = true, blue = estimated)')
#axis([0, Tlim, 1.2, 1.8])

#subplot(2,1,2)
#plot(tobs, Tvectru, 'r-', tobs, Tvecest, 'b-')
#ylabel('\fontsize{16} T(t)')
#title('\fontsize{16} Temperature')
#axis([0, Tlim, 330, 360])

plot(tobs, Tvectru, type="l", col="red", xlab="", ylab="C(t)",
     main="Temperature", ylim=c(330, 360))
lines(tobs, Tvecest, col="blue")
points(tobs, Tobs, col="blue")

par(op)

##
##%%  13.  Fit the hot solution with the cool parameter estimates
##
#%  replace true by estimated parameter values

#fitstruct.kref   = parvec(1);      
#fitstruct.EoverR = parvec(2);   
#fitstruct.a      = parvec(3);   
#fitstruct.b      = parvec(4);

fitStrHot.nls <- fitStrHot
parvecCoolEst <- coef(nlsFit)
fitStrHot.nls$kref <- parvecCoolEst[1]
fitStrHot.nls$EoverR <- parvecCoolEst[2]
fitStrHot.nls$a <- parvecCoolEst[3]
fitStrHot.nls$b <- parvecCoolEst[4]

#yinit = [Cinit_hot, Tinit_hot];
# yinitHot from 5.4 above 

#%  Matlab:  load hot input into fitstruct
#   R:  Retained from 5.4 above in fitStrHot 
#fitstruct.Tcin = Tcin_hot;
all.equal(fitStrHot$Tcin, hotStepInput[, "Tcin"])

#%  solve  differential equation with true parameter values

#h = waitbar(0,'Simulating Actual Model Output...');
#[t, yest_hot] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
#                            fitstruct, 'all_hot_step', Tlim);
#close(h)

hotStepSoln.nls <- lsoda(y=yinitHot, times=tspan, func=CSTR2, parms=
  list(fitstruct=fitStrHot.nls, condition='all.hot.step', Tlim=Tlim) )

#%  set up separate variables for concentration and temperature

#C_hotest = yest_hot(:,1);
#T_hotest = yest_hot(:,2);
str(hotStepSoln.nls)
C.hotest <- hotStepSoln.nls[, "Conc"]
T.hotest <- hotStepSoln.nls[, "Temp"]

#Cerrpct = 100.*median(abs(C_hotest - C_hot)./Cinit_hot)
#Terrpct = 100.*median(abs(T_hotest - T_hot)./Tinit_hot)

100.*quantile(abs(C.hotest - C_hot)/Cinit_hot)
#       0%       25%       50%       75%      100% 
# 0.000000  1.438755  1.802374  2.507454 11.481121
#Cerrpct = 100.*median(abs(Cvecest - Cvectru)/Cvectru)

100.*quantile(abs(T.hotest - T_hot)/Tinit_hot)
#        0%        25%        50%        75%       100% 
#0.00000000 0.03298232 0.05759536 0.08835659 0.53245573 

#% plot the estimated and true hot solutions

#figure(8)
#subplot(2,1,1)

op <- par(mfrow=c(2,1))
plot(tspan, C_hot, type="l", col="red", ylab="C(t)", xlab="",
     main="Concentration (red = true, blue = estimates)")
lines(tspan, C.hotest, col="blue")

plot(tspan, T_hot, type="l", col="red", ylab="T(t)", xlab="",
     main="Temperature")
lines(tspan, T.hotest, col="blue")

par(op)
# Matlab plot suggests it got a better fit for this.  

#lhdl = plot(t, C_hot, 'r-'); 
#set(lhdl, 'LineWidth', 1);
#lhdl = line(t, C_hotest);
#set(lhdl, 'LineWidth', 1, 'color', 'b');
#ylabel('\fontsize{16} C(t)')
#title('\fontsize{16} Concentration (red = true, blue = estimated)')
#axis([0, Tlim, 0.1, 0.48])
#subplot(2,1,2)
#lhdl = plot(t, T_hot, 'r-'); 
#set(lhdl, 'LineWidth', 1);
#lhdl = line(t, T_hotest);
#set(lhdl, 'LineWidth', 1, 'color', 'b');
#ylabel('\fontsize{16} T(t)')
#title('\fontsize{16} Temperature')
#axis([0, Tlim, 370, 420])

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

start.time <- proc.time()
nlsFitTemp <- nls(formula=~CSTRres(kref=kref, EoverR=EoverR, a=a, b=b,
   datstruct=datstr, fitstruct=fitstr, CSTRbasis=Cb, lambda=lam),
              data=list(datstr=datstruct, fitstr=fitStrHTemp,
                             Cb=CSTRbasis, lam=lambda), 
    start=as.list(parvec0),
    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE)
(et <- proc.time()-start.time)
#1105510 :  0.46100 0.83301 1.67800 0.50000 
#  It.   1, fac=           1, eval (no.,total): ( 1,  1):Warning message:
#singular gradient in: nls(formula = ~CSTRres(kref = kref, EoverR = EoverR, a = a, b = b,

#**** NO SOLUTION:  quit at starting values
#**** because 'singular gradient' estimated.  


















#%  get final solution values

#[res, Dres, fitstruct] = CSTRfn(parvec, datstruct, fitstruct, ...
#                                         CSTRbasis, lambda, 0);

#%  display the parameter estimates

#disp(['Initial   values: ', num2str(parvec0)])
#disp(['Estimated values: ', num2str(parvec)])
#disp(['True      values: ', num2str(parvectru)])

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

start.time <- proc.time()
nlsFitConc <- nls(formula=~CSTRres(kref=kref, EoverR=EoverR, a=a, b=b,
   datstruct=datstr, fitstruct=fitstr, CSTRbasis=Cb, lambda=lam),
              data=list(datstr=datstruct, fitstr=fitStrHConc,
                             Cb=CSTRbasis, lam=lambda), 
    start=as.list(parvec0),
    control=nls.control(printEval=TRUE, warnOnly=TRUE), trace=TRUE)
(et.Conc <- proc.time()-start.time)
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
nlsFitConc. <- CSTRfn(parvec=coef(nlsFitConc), datstruct,
                      fitStrHConc, CSTRbasis, lambda)

#%  display the parameter estimates

#disp(['Initial   values: ', num2str(parvec0)])
#disp(['Estimated values: ', num2str(parvec)])
#disp(['True      values: ', num2str(parvectru)])
rbind("Initial values"=parvec0,
      "Estimated values"=coef(nlsFitConc),
      "True valules"=parvectru)

#%  get the final coefficient values for the solutions

#coef     = fitstruct.coef0;
coefConc <- nlsFitConc.$fitstruct$coef0
#Ccoefest = coef(1:nbasis);
CcoefConc.est <- coefConc[1:nbasis]
#Tcoefest = coef(nbasis+1:2*nbasis);
TcoefConc.est <- coefConc[(nbasis+1):(2*nbasis)]

#Cfdest = fd(Ccoefest, CSTRbasis);
Cfd.ConcEst <- fd(CcoefConc.est, CSTRbasis)
#Tfdest = fd(Tcoefest, CSTRbasis);
Tfd.ConcEst <- fd(TcoefConc.est, CSTRbasis)

#Cvecest = eval_fd(tobs, Cfdest);
CvecConcEst <- eval.fd(tobs, Cfd.ConcEst)
#Tvecest = eval_fd(tobs, Tfdest);
TvecConcEst <- eval.fd(tobs, Tfd.ConcEst)

#Cvectru = eval_fd(tobs, Cfdtru)
#Tvectru = eval_fd(tobs, Tfdtru);
# Conputed above following 'nlsFit' computation 

#%  display maximum absolute errors

quantile((CvecConcEst-Cvectru)/Cvectru)
#           0%           25%           50%           75%          100% 
#-0.0109490925 -0.0014349335 -0.0006695048  0.0003116547  0.0029172212 
#Cerrpct = 100.*median(abs(Cvecest - Cvectru)./Cvectru)
CerrConcPct <- 100*median(abs((CvecConcEst-Cvectru)/Cvectru))
# 0.0934

quantile((TvecConcEst-Tvectru)/Tvectru)
#           0%           25%           50%           75%          100% 
#-4.597640e-04 -9.884341e-05  1.719107e-07  1.619007e-04  1.755867e-03 
#Terrpct = 100.*median(abs(Tvecest - Tvectru)./Tvectru)
TerrConcPct <- 100*median(abs((TvecConcEst-Tvectru)/Tvectru))
# 0.0124

#%  plot the estimated and true solutions

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

op <- par(mfrow=c(2,1))
plot(tspan, C_hot, type="l", col="red", ylab="C(t)", xlab="",
     main="Concentration (red = true, blue = estimates)")
lines(tspan, C.hotest, col="blue")

plot(tspan, T_hot, type="l", col="red", ylab="T(t)", xlab="",
     main="Temperature")
lines(tspan, T.hotest, col="blue")

par(op)

#%%  15.  Fit the hot solution with the cool parameter estimates

#%  replace true by estimated parameter values

#fitstruct.kref   = parvec(1);      
#fitstruct.EoverR = parvec(2);   
#fitstruct.a      = parvec(3);   
#fitstruct.b      = parvec(4);   

#yinit = [Cinit_hot, Tinit_hot];

#%  load hot input into fitstruct

#fitstruct.Tcin = Tcin_hot;

#%  solve  differential equation with true parameter values

#h = waitbar(0,'Simulating Actual Model Output...');
#[t, yest_hot] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
#                            fitstruct, 'all_hot_step', Tlim);
#close(h)

#%  set up separate variables for concentration and temperature

#C_hotest = yest_hot(:,1);
#T_hotest = yest_hot(:,2);

#% plot the estimated and true hot solutions

#figure(10)
#subplot(2,1,1)
#plot(t, C_hot, 'r-', t, C_hotest, 'b-')
#ylabel('\fontsize{16} C(t)')
#title('\fontsize{16} Concentration (red = true, blue = estimated)')
#axis([0, Tlim, 0, 0.5])
#subplot(2,1,2)
#plot(t, T_hot, 'r-', t, T_hotest, 'b-')
#ylabel('\fontsize{16} T(t)')
#title('\fontsize{16} Temperature')
#axis([0, Tlim, 370, 420])

#tst <- c(1,2,3,2)
#tst2 <- rep(tst, 4)
#tst3 <- rep(tst, c(1, 4, 8, 4))

#tst; tst2; tst3
#sd(tst)
#sd(tst2)
#sd(tst3)

