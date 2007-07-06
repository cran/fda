
%%  1.  Introduction

%  Estimate the four parameters defining two differential equations
%  that describe the behavior of the output concentration
%  and temperature for a nonisotherman continuously stirred 
%  tank reactor, abbreviated nonisothermal CSTR.  
%  The system of two nonlinear equations has five forcing or  
%  input functions.
%  These equations are taken from
%  Marlin, T. E. (2000) Process Control, 2nd Edition, McGraw Hill,
%  pages 899-902.

%  Last modified 16 April 2007 by Spencer Graves
%  Previously modified 26 October 2005

%addpath('c:/matlab7/fdaM')
addpath('d:/spencerg/statmtds/fda/Matlab/fdaM')

%%  2.  Set up the problem 

%  set up the system of equations with known constants
%  and values of the four parameters

%  store known system constants in a struct variable fitstruct

fitstruct.V    = 1.0;       %  volume in m^3
fitstruct.Cp   = 1.0;       %  concentration in cal/(g.K)
fitstruct.rho  = 1.0;       %  density in g/m^3
fitstruct.delH = -130.0;    %  cal/kmol
fitstruct.Cpc  = 1.0;       %  concentration in cal/(g.K)
fitstruct.rhoc = 1.0;       %  cal/kmol
fitstruct.Tref = 350;       % reference temperature

% store true values of known parameters 

EoverRtru = 0.83301;     %  E/R in units K/1e4
kreftru   = 0.4610 ;      %  reference value
atru      = 1.678;       %  a in units (cal/min)/K/1e6
btru      = 0.5;         %  dimensionless exponent

% enter these parameter values into fitstruct

fitstruct.kref   = kreftru;      
fitstruct.EoverR = EoverRtru;   %  kref = 0.4610
fitstruct.a      = atru;        %  a in units (cal/min)/K/1e6
fitstruct.b      = btru;        %  dimensionless exponent

%% 3.  Choose the input design

%  Outputs have oscillating responses to lowering 
%  the temperature of the coolant when the tank is
%  operating at a high temperature, around 365 deg K.
%  Outputs respond smoothly to changes when the tank
%  is operating at a cooler temperature, around 330 deg K. 

%  These two designs differ in terms of the temperature
%  at which the reactor is operated.  The "cool" design
%  using a coolant temperature of 330 deg K, and the
%  outputs respond smoothly to step changes in input.
%  The "hot" design uses a coolant temperature of 365 deg K,
%  and both concentration and temperature show sharp
%  oscillations when an input is changed, and especially
%  when coolant temperature is dropped, at which point
%  the reactor is nearly unstable.  

%  set up time values at which output values are to be
%  computed in the numerical solution to the equation.
 
Tlim  = 64;   % reaction observed over interval [0, Tlim]
delta = 1/12; % observe every five seconds
tspan = (0:delta:Tlim);  
nspan = length(tspan);

%%  4.  Display the two input conditions

%  set up input function data for both cool and hot
%  operations, and plot inputs

%[F,CA0,T0,Tcin_cool,Fc] = CSTR2in(tspan, 'all_cool_step');
[F,CA0,T0,Tcin_cool] = CSTR2in(tspan, 'all_cool_step');
[F,CA0,T0,Tcin_hot, Fc] = CSTR2in(tspan, 'all_hot_step');

%  plot inputs

figure(1)
subplot(5,1,1)
phdl = plot(tspan, F);
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{13} F(t)')
title('\fontsize{13} Input flow rate')
axis([0,Tlim,0.4,1.6])
subplot(5,1,2)
phdl = plot(tspan, CA0);
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{13} C_0(t)')
title('\fontsize{13} Input concentration')
axis([0,Tlim,1.7,2.3])
subplot(5,1,3)
phdl = plot(tspan, T0);
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{13} T_0(t)')
title('\fontsize{13} Input temperature')
axis([0,Tlim,300,350])
subplot(5,1,4)
phdl = plot(tspan, Tcin_cool, 'b', ...
            tspan, Tcin_hot,  'r');
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{13} T_{cin}(t)')
title('\fontsize{13} Coolant temperature (red = hot, blue = cool)')
axis([0,Tlim,325,375])
subplot(5,1,5)
phdl = plot(tspan, Fc);
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{13} F_c(t)')
title('\fontsize{13} Coolant flow rate')
axis([0,Tlim,8,22])

%%  5.  Solve the equations for both the hot and cool conditions

%  5.0.  set constants for ODE solver

odeoptions = odeset('RelTol',1e-7,'AbsTol',1e-7);

%  5.1.  cool condition solution

%  initial conditions

Cinit_cool = 1.5965;    %  initial concentration in kmol/m^3
Tinit_cool = 341.3754;  %  initial temperature in deg K
yinit = [Cinit_cool, Tinit_cool];

%  load cool input into fitstruct

fitStrCool = fitstruct; 
fitStrCool.Tcin = Tcin_cool;

%  5.2.  solve  differential equation with true parameter values

h = waitbar(0,'Simulating Cool Output...');
[t, y_cool] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
                            fitStrCool, 'all_cool_step', Tlim);
close(h)

%  5.3.  set up separate variables for concentration and temperature

C_cool = y_cool(:,1);
T_cool = y_cool(:,2);

%  5.4.  hot condition solution

%  initial conditions

Cinit_hot  = 0.2651;    %  initial concentration in kmol/m^3
Tinit_hot  = 394.0532;  %  initial temperature in deg K
yinit = [Cinit_hot, Tinit_hot];

%  load hot input into fitstruct

fitStrHot = fitstruct; 
fitStrHot.Tcin = Tcin_hot;

%  solve  differential equation with true parameter values

h = waitbar(0,'Simulating Hot Output...');
[t, y_hot] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
                            fitStrHot, 'all_hot_step', Tlim);
close(h)

%  set up separate variables for concentration and temperature

C_hot = y_hot(:,1);
T_hot = y_hot(:,2);

%  5.5.  plot deterministic model outputs

figure(2)
subplot(2,1,1)
lhdl = plot(t, C_cool, 'b-', t, C_hot, 'r-');
set(lhdl, 'LineWidth', 2);
hold on
lhdl = plot([44,44],[0,2], 'm--');
set(lhdl, 'LineWidth', 1);
lhdl = plot([48,48],[0,2], 'm--');
set(lhdl, 'LineWidth', 1);
hold off
title('Actual Model Output')       
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Output concentration (red = hot, blue = cool)')
axis([0,Tlim,0,2])
subplot(2,1,2)
lhdl = plot(t, T_cool, 'b-', t, T_hot, 'r-');
set(lhdl, 'LineWidth', 2);
hold on
lhdl = plot([44,44],[330,420], 'm--');
set(lhdl, 'LineWidth', 1);
lhdl = plot([48,48],[330,420], 'm--');
set(lhdl, 'LineWidth', 1);
hold off
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Output temperature')
axis([0,Tlim,330,420])


%%  6.  Set up the B-spline basis for representing solutions.

%  The splines must have discontinuous derivatives and 
%  continuous function values at points where there are 
%  step changes in an input.

%  this set up is for the cool operation.  The hot
%  operation requires twice as many knots.

brk    = (4:4:60);  %  times of step changes
knots  = 0:1/3:Tlim;  %  knots at every 20 secs.
knots  = sort([knots, brk, brk]);
nknots = length(knots);  %  223 knots

%  set up the basis for the expansion as an order 4 
%  B-spline basis with these knots 

norder = 4;
nbasis = length(knots) + norder - 2;
CSTRbasis0 = create_bspline_basis([0,Tlim],nbasis,norder,knots);

%  set up quadrature points and weights for estimating
%  the integrals in the roughness penalties

nquadint = 5;
[CSTRbasis, quadpts, quadwts] = quadset(nquadint, CSTRbasis0);

nquad = length(quadpts);  %  4032 quadrature points

%%  7.  Compute fits to error less data and display results

%  7.1.  compute the expansion for the errorless data 

fdParobj  = fdPar(CSTRbasis, 2, 1e-8);

Cfdtru = smooth_basis(t, C_cool, fdParobj);
Tfdtru = smooth_basis(t, T_cool, fdParobj);

% Transfer smooth_basis results to R 

save CSTRfdtru -v6 Cfdtru Tfdtru; 
% load CSTRfdtru ;

%  7.2.  plot the fit to the solution using the expansion

figure(3)
subplot(2,1,1)
plotfit_fd(C_cool, t, Cfdtru)
xlabel('')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plotfit_fd(T_cool, t, Tfdtru)
xlabel('')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

%  7.3.  compute weights to correct for differences in
%  scale for the two variables

Cvectru = eval_fd(t, Cfdtru);
Cwt     = var(Cvectru);

Tvectru = eval_fd(t, Tfdtru);
Twt     = var(Tvectru);

wt   = [Cwt, Twt];

%%   8.  Set up some data with noise

tobs = 1:1/3:Tlim;  %  observations at every 20 secs.
N = length(tobs);   %  190 observations

%  8.1.  fix standard errors 

stderrfac = 0.2;  %  standard error of measurements
stderrC   = stderrfac*sqrt(Cwt);
stderrT   = stderrfac*sqrt(Twt);

%  8.2.  compute errorless values

Cvec = eval_fd(tobs, Cfdtru);
Tvec = eval_fd(tobs, Tfdtru);

%  8.3.  add errors

Cobs    = Cvec + randn(N,1).*stderrC;
Tobs    = Tvec + randn(N,1).*stderrT;

% Save the simulations to be read by R 
% to make it easier to compare answers  
save CSTRsim -v6 Cobs Tobs; 
% load CSTRsim; 

%  8.4.  smooth the noisy data using smoothing splines

lambdaC = 1e0;
CfdPar  = fdPar(CSTRbasis, 2, lambdaC);
Cfdsmth = smooth_basis(tobs, Cobs, CfdPar);

lambdaT = 1e-1;
TfdPar  = fdPar(CSTRbasis, 2, lambdaT);
Tfdsmth = smooth_basis(tobs, Tobs, TfdPar);

% Save the Smooths to be read by R 
% to make it easier to compare answers  
save CSTRsmth -v6 CfdPar Cfdsmth TfdPar Tfdsmth; 
% load CSTRsmth ; 

%  8.5.  plot the data

figure(4)
subplot(2,1,1)
plot(tobs, Cobs, '.');
xlabel('')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plot(tobs, Tobs, '.');
xlabel('')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

%  8.6.  plot the data and the smooths

figure(5)
subplot(2,1,1)
plotfit_fd(Cobs, tobs, Cfdsmth);
xlabel('')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plotfit_fd(Tobs, tobs, Tfdsmth);
xlabel('')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

%%  9.  Load data into struct variable datstruct

%  this is information that only depends on the
%  data and sampling design

%  9.1.  weights for variables

datstruct.Cwt = Cwt;
datstruct.Twt = Twt;

%  9.2.  data

datstruct.y = [Cobs, Tobs];

%  9.3.  basis values at sampling points

basismat            = eval_basis(tobs, CSTRbasis);
Dbasismat           = eval_basis(tobs, CSTRbasis, 1);
datstruct.basismat  = basismat;
datstruct.Dbasismat = Dbasismat;

%  9.4.  forcing function values at quadrature points

[Fquad, CA0quad, T0quad, Tcinquad, Fcquad] = ...
                     CSTR2in(quadpts, 'all_cool_step');

datstruct.F    = Fquad;
datstruct.CA0  = CA0quad;
datstruct.T0   = T0quad;
datstruct.Tcin = Tcinquad;
datstruct.Fc   = Fcquad;

%  9.5.  basis values at quadrature points

quadbasismat            = eval_basis(quadpts, CSTRbasis);
Dquadbasismat           = eval_basis(quadpts, CSTRbasis, 1);
datstruct.quadpts       = quadpts;
datstruct.quadwts       = quadwts;
datstruct.quadbasismat  = quadbasismat;
datstruct.Dquadbasismat = Dquadbasismat;

%%  10.  Analyze the noisy data to estimate all the parameters.

%  This is a very difficult problem because all four
%  parameters affect the speed of the reaction in nearly
%  identical ways.  Convergence is slow.  

%% 11.  Preliminary explorations of sampling variances

%  Compute a solution using the true parameters as a 
%  starting point and zero iterations

%  11.1.  Specify that both output variables are measured

fitStrHot.fit = [1,1];

%  11.2.  Use coefficients from spline smoothing as initial
%  values for coefficient vectors

Ccoef = getcoef(Cfdsmth);
Tcoef = getcoef(Tfdsmth);
fitStrHot.coef0 = [Ccoef; Tcoef];

%  11.3.  Specify that all parameters will be estimated

fitStrHot.estimate = [1, 1, 1, 1];  
estind = find(fitStrHot.estimate == 1);

%  11.4.  set up the true parameter values

parvec0 = [kreftru, EoverRtru, atru, btru];

%  11.5.  specify smoothing parameter values

lambda = 1e1.*[lambdaC, lambdaT];

%  11.6.  calculate the estimates

[res, jacobian] = CSTRfn(parvec0, datstruct, fitStrHot, ...
                         CSTRbasis, lambda);
                     
% Save to compare with R 
save CSTR1 -v6 res jacobian; 
                     
%  get the singular values of the jacobian

s = svd(full(jacobian));

%  log 10 of condition number

log10(s(1)/s(4))
% 3.019 on 2007.05.30 
% previously:  
%  3.32 for multiplier of 1e0
%  3.11 for multiplier of 1e1
%  3.11 for multiplier of 1e2

%  Sampling variances will vary by over 6 orders of magnitude!

%  11.7.  Use only the first two parameters

fitStrH12 = fitStrHot ; 
fitStrH12.estimate = [1, 1, 0, 0];  
estind12 = find(fitStrH12.estimate == 1);

parvec12 = [kreftru, EoverRtru];

%    calculate the estimates

% second use of CSTRfn 
[res12, jacobian12] = CSTRfn(parvec12, datstruct, fitStrH12, ...
                         CSTRbasis, lambda);
% Save for comparison with R 
save CSTR12 -v6 res12 jacobian12; 

%    get the singular values of the jacobian

s12 = svd(full(jacobian12));

%  log 10 of condition number

log10(s12(1)/s12(2))
% 0.8261 on 2007.05.30 
% previous number 0.7393

%  11.8.  Use only the first and third parameters

fitStrH13 = fitStrHot ; 
fitStrH13.estimate = [1, 0, 1, 0];  
estind13 = find(fitStrH13.estimate == 1);

parvec13 = [kreftru, atru];

%    calculate the estimates

%# CSTRfn use 3 
[res13, jacobian13] = CSTRfn(parvec13, datstruct, fitStrH13, ...
                         CSTRbasis, lambda);
                                     
%    get the singular values of the jacobian

s13 = svd(full(jacobian13));

%  log 10 of condition number

log10(s13(1)/s13(2))
% 1.1789 on 2007.05.30 
% previously  0.0610

%  11.9.  Use only the third and fourth parameters

fitStrH34 = fitStrHot; 
fitStrH34.estimate = [0, 0, 1, 1];  
estind34 = find(fitStrH34.estimate == 1);

parvec34 = [atru, btru];

%  calculate the estimates
% CSTRfn use 4 
[res34, jacobian34] = CSTRfn(parvec34, datstruct, fitStrH34, ...
                         CSTRbasis, lambda);
                                     
%  get the singular values of the jacobian

s34 = svd(full(jacobian34));

%  log 10 of condition number

log10(s34(1)/s34(2))
% 2.1200 on 2007.05.30 
% previously 3.0344

%  11.10.  Use only the first, second and third parameters

fitStrH123 = fitStrHot ; 
fitStrH123.estimate = [1, 1, 1, 0];  
estind123 = find(fitStrH123.estimate == 1);

parvec123 = [kreftru, EoverRtru, atru];

%  calculate the estimates
% CSTRfn use 5 
[res123, jacobian123] = CSTRfn(parvec123, datstruct, fitStrH123, ...
                         CSTRbasis, lambda);
                                                          
%  get the singular values of the jacobian

s123 = svd(full(jacobian123));

%  log 10 of condition number

log10(s123(1)/s123(3))
% 1.2357 on 2007.06.03 by Spencer Graves
% previously 0.7495

save CSTR5 -v6 ;
% load CSTR5 ; 
%%  12.  Estimate all four parameters

%  set the optimization parameters for the outer optimization
%  only function values are used at this point

tolval = 1e-8;
optionsCSTRfn = optimset('LargeScale', 'on', 'Display', ...
                          'iter', 'MaxIter', 50, ...
                          'TolCon', tolval, 'TolFun', tolval, ...
                          'TolX',   tolval, 'TolPCG', tolval, ...
                          'Jacobian', 'on');
                      
%  use coefficients from spline smoothing as initial
%  values for coefficient vectors

Ccoef = getcoef(Cfdsmth);
Tcoef = getcoef(Tfdsmth);

%fitstruct.coef0 = [Ccoef; Tcoef];
fitStrH1234 = fitStrHot ; 
fitStrH1234.coef0 = [Ccoef; Tcoef];

%  Use only the first, second and third parameters

fitStrH1234.fit = [1,1];


%  set up some initial values for parameters

parvec1234 = [0.4, 0.8, 1.7, 0.5];

%  set up the true parameter values

parvectru = [kreftru, EoverRtru, atru, btru];

%  specify smoothing parameter values

lambda = 10.*[lambdaC, lambdaT];

%  calculate the estimates
% CSTRfn use 6 
tic;
[parvec1234, resnorm1234, residual1234, exitflag1234, output1234, ...
    lambdaout1234, jacobian1234] = ...
    lsqnonlin(@CSTRfn, parvec1234, [], [], optionsCSTRfn, ...
              datstruct, fitStrH1234, CSTRbasis, lambda);
et = toc

save cstr4nls  -v6 parvec1234 resnorm1234 residual1234 ... 
    exitflag1234 output1234 lambdaout1234 jacobian1234 ; 
%  get the singular values of the jacobian

% test CSTRfn again 
parvecTst = [.461, .83301, 1.678, .5]+0.001*pi.^(1./(1:4)) ; 

[resTst jacobianTst] = CSTRfn(...
    parvecTst, datstruct, fitStrH1234, CSTRbasis, lambda) ; 

svd(full(jacobianTst)) 

tic;
[parvec1234o, resnorm1234o, residual1234o, exitflag1234o, output1234o, ...
    lambdaout1234o, jacobian1234o] = ...
    lsqnonlin(@CSTRfn, parvec1234, [], [], optionsCSTRfn, ...
              datstruct, fitStrH1234, CSTRbasis, lambda);
et = toc

s1234 = svd(full(jacobian1234));

%  log 10 of condition number

%log10(s1234(1)/s1234(3))
log10(s1234(1)/s1234(4))
% 3.0146

%  get final solution values

[resAll, DresAll, fitstructAll] = CSTRfn(parvec1234, datstruct, ...
        fitStrH1234, CSTRbasis, lambda, 0);

%  display the parameter estimates

disp(['Initial   values: ', num2str(parvec1234)])
disp(['Estimated values: ', num2str(parvec)])
disp(['True      values: ', num2str(parvectru)])

%  get the final coefficient values for the solutions

coef     = fitStrH1234.coef0;
Ccoefest = coef(1:nbasis);
Tcoefest = coef(nbasis+1:2*nbasis);

Cfdest = fd(Ccoefest, CSTRbasis);
Tfdest = fd(Tcoefest, CSTRbasis);

Cvecest = eval_fd(tobs, Cfdest);
Tvecest = eval_fd(tobs, Tfdest);

Cvectru = eval_fd(tobs, Cfdtru);
Tvectru = eval_fd(tobs, Tfdtru);

%  display maximum absolute errors

Cerrpct = 100.*median(abs(Cvecest - Cvectru)./Cvectru)
Terrpct = 100.*median(abs(Tvecest - Tvectru)./Tvectru)

%  plot the estimated and true solutions

figure(6)
subplot(2,1,1)
plot(tobs, Cvectru, 'r-', tobs, Cvecest, 'b-', tobs, Cobs, 'b.')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plot(tobs, Tvectru, 'r-', tobs, Tvecest, 'b-')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

%%  13.  Fit the hot solution with the cool parameter estimates

%  replace true by estimated parameter values

%fitstruct.kref   = parvec(1);      
%fitstruct.EoverR = parvec(2);   
%fitstruct.a      = parvec(3);   
%fitstruct.b      = parvec(4);   

fitStrHC = fitStrHot; 
fitStrHC.kref   = parvec1234(1);      
fitStrHC.EoverR = parvec1234(2);   
fitStrHC.a      = parvec1234(3);   
fitStrHC.b      = parvec1234(4);   

yinit = [Cinit_hot, Tinit_hot];

%  load hot input into fitstruct

fitStrHC.Tcin = Tcin_hot;

%  solve  differential equation with true parameter values

h = waitbar(0,'Simulating Actual Model Output...');
[t, yest_hot] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
                            fitStrHC, 'all_hot_step', Tlim);
close(h)

%  set up separate variables for concentration and temperature

C_hotest = yest_hot(:,1);
T_hotest = yest_hot(:,2);

Cerrpct = 100.*median(abs(C_hotest - C_hot)./Cinit_hot)
Terrpct = 100.*median(abs(T_hotest - T_hot)./Tinit_hot)

% plot the estimated and true hot solutions

figure(8)
subplot(2,1,1)
lhdl = plot(t, C_hot, 'r-'); 
set(lhdl, 'LineWidth', 1);
lhdl = line(t, C_hotest);
set(lhdl, 'LineWidth', 1, 'color', 'b');
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 0.1, 0.48])
subplot(2,1,2)
lhdl = plot(t, T_hot, 'r-'); 
set(lhdl, 'LineWidth', 1);
lhdl = line(t, T_hotest);
set(lhdl, 'LineWidth', 1, 'color', 'b');
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 370, 420])

%%  14.  Estimate kref, EoverR and a with only 
%   one response variable observed.  

% 14.1.  With only temperature observed 

fitStrHTemp = fitStrHot
fitStrHTemp.fit = [0,1];

%  use coefficients from spline smoothing as initial
%  values for coefficient vectors

fitStrHTemp.coef0 = [Ccoef; Tcoef];

%  set up the true parameter values

parvectru = [kreftru, EoverRtru, atru, btru];
fitStrHTemp.b = btru;

%  calculate the estimates

tic;
[parvecT, resnormT, residualT, exitflagT] = ...
    lsqnonlin(@CSTRfn, parvec1234, [], [], optionsCSTRfn, ...
              datstruct, fitStrHTemp, CSTRbasis, lambda);
toc

%  get final solution values

[resT, DresT, fitStrFinalT] = CSTRfn(parvec, datstruct, ... 
                    fitStrHTemp, CSTRbasis, lambda, 0);

%  display the parameter estimates

disp(['Initial   values: ', num2str(parvec1234)])
disp(['Estimated values: ', num2str(parvecT)])
disp(['True      values: ', num2str(parvectru)])

%  get the final coefficient values for the solutions

coefT     = fitStrFinalT.coef0;
CcoefestT = coef(1:nbasis);
TcoefestT = coef(nbasis+1:2*nbasis);

CfdestT = fd(CcoefestT, CSTRbasis);
TfdestT = fd(TcoefestT, CSTRbasis);

CvecestT = eval_fd(tobs, CfdestT);
TvecestT = eval_fd(tobs, TfdestT);

CvectruT = eval_fd(tobs, CfdtruT);
TvectruT = eval_fd(tobs, TfdtruT);

%  display maximum absolute errors

CerrpctT = 100.*median(abs(Cvecest - Cvectru)./Cvectru)
TerrpctT = 100.*median(abs(Tvecest - Tvectru)./Tvectru)

%  plot the estimated and true solutions

figure('9T')
subplot(2,1,1)
plot(tobs, Cvectru, 'r-', tobs, CvecestT, 'b-')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plot(tobs, Tvectru, 'r-', tobs, TvecestT, 'b-')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

% 14.2.  With only concentration observed 

fitStrHConc = fitStrHot
fitStrHConc.fit = [1, 0];

%  use coefficients from spline smoothing as initial
%  values for coefficient vectors

fitStrHConc.coef0 = [Ccoef; Tcoef];

%  set up the true parameter values

parvectru = [kreftru, EoverRtru, atru, btru];
fitStrHConc.b = btru;

%  calculate the estimates

tic;
[parvecC, resnormC, residualC, exitflagC] = ...
    lsqnonlin(@CSTRfn, parvec1234, [], [], optionsCSTRfn, ...
              datstruct, fitStrHConc, CSTRbasis, lambda);
toc

%  get final solution values

[resC, DresC, fitStrFinalC] = CSTRfn(parvec, datstruct, ... 
                fitStrHConc, CSTRbasis, lambda, 0);

%  display the parameter estimates

disp(['Initial   values: ', num2str(parvec1234)])
disp(['Estimated values: ', num2str(parvecC)])
disp(['True      values: ', num2str(parvectru)])

%  get the final coefficient values for the solutions

coefC     = fitStrFinalC.coef0;
CcoefestC = coefC(1:nbasis);
TcoefestC = coefC(nbasis+1:2*nbasis);

CfdestC = fd(CcoefestC, CSTRbasis);
TfdestC = fd(TcoefestC, CSTRbasis);

CvecestC = eval_fd(tobs, CfdestC);
TvecestC = eval_fd(tobs, TfdestC);

CvectruC = eval_fd(tobs, CfdtruC);
TvectruC = eval_fd(tobs, TfdtruC);

%  display maximum absolute errors

CerrpctC = 100.*median(abs(Cvecest - Cvectru)./Cvectru)
TerrpctC = 100.*median(abs(Tvecest - Tvectru)./Tvectru)

%  plot the estimated and true solutions

figure('9C')
subplot(2,1,1)
plot(tobs, Cvectru, 'r-', tobs, CvecestC, 'b-')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plot(tobs, Tvectru, 'r-', tobs, TvecestC, 'b-')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

%%  15.  Fit the hot solution with the cool parameter estimates

%  replace true by estimated parameter values

%fitstruct.kref   = parvec(1);      
%fitstruct.EoverR = parvec(2);   
%fitstruct.a      = parvec(3);   
%fitstruct.b      = parvec(4);   

fitStrHTC = fitStrHTemp; 
fitStrHTC.kref   = parvec(1);      
fitStrHTC.EoverR = parvec(2);   
fitStrHTC.a      = parvec(3);   
fitStrHTC.b      = parvec(4);   

yinit = [Cinit_hot, Tinit_hot];

%  load hot input into fitstruct

fitStrHTC.Tcin = Tcin_hot;

%  solve  differential equation with true parameter values

h = waitbar(0,'Simulating Actual Model Output...');
[t, yest_hot] = ode45(@CSTR2, tspan, yinit, odeoptions, ...
                            fitStrHTC, 'all_hot_step', Tlim);
close(h)

%  set up separate variables for concentration and temperature

C_hotest = yest_hot(:,1);
T_hotest = yest_hot(:,2);

% plot the estimated and true hot solutions

figure(10)
subplot(2,1,1)
plot(t, C_hot, 'r-', t, C_hotest, 'b-')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 0, 0.5])
subplot(2,1,2)
plot(t, T_hot, 'r-', t, T_hotest, 'b-')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 370, 420])



