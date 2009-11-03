
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
R_HOME = getenv('R_HOME') ; % Where is R installed?  
Rpath = R_HOME(1:(end-4)) ; % Drop '\bin' 
fdaPath = fullfile(Rpath, 'library\fda') ; 
fdaM = fullfile(fdaPath, 'Matlab\fdaM') ; 
addpath(fdaM) ; 
CSTRpath = fullfile(fdaPath, 'scripts\CSTR') ; 
addpath(CSTRpath) ; 

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

save CSTR_fig1Matlab -v6; 

%%  5.  Solve the equations for both the hot and cool conditions

%  5.0.  set constants for ODE solver

odeoptions = odeset('RelTol',1e-7,'AbsTol',1e-7);

%  5.1.  cool condition solution

%  initial conditions

Cinit_cool = 1.5965;    %  initial concentration in kmol/m^3
Tinit_cool = 341.3754;  %  initial temperature in deg K
yinitCool = [Cinit_cool, Tinit_cool];

%  load cool input into fitstruct

fitStrCool = fitstruct; 
fitStrCool.Tcin = Tcin_cool;

%  5.2.  solve  differential equation with true parameter values

h = waitbar(0,'Simulating Cool Output...');
[t, y_cool] = ode45(@CSTR2, tspan, yinitCool, odeoptions, ...
                            fitStrCool, 'all_cool_step', Tlim);
close(h)

%  5.3.  set up separate variables for concentration and temperature

C_cool = y_cool(:,1);
T_cool = y_cool(:,2);

%  5.4.  hot condition solution

%  initial conditions

Cinit_hot  = 0.2651;    %  initial concentration in kmol/m^3
Tinit_hot  = 394.0532;  %  initial temperature in deg K
yinitHot = [Cinit_hot, Tinit_hot];

%  load hot input into fitstruct

fitStrHot = fitstruct; 
fitStrHot.Tcin = Tcin_hot;

%  solve  differential equation with true parameter values

h = waitbar(0,'Simulating Hot Output...');
[t, y_hot] = ode45(@CSTR2, tspan, yinitHot, odeoptions, ...
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

save CSTR_fig2Matlab -v6; 

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

save CSTR_fig3Matlab -v6; 

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

save CSTR_fig4Matlab -v6; 


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

save CSTR_fig5Matlab -v6; 

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
% manual check of jacobian(:, 4) 
res4=CSTRfn(parvec0+[0 0 0 .01],datstruct,fitStrHot,CSTRbasis,lambda);
jacob4 = (res4-res)/0.01 ; 
figure(6) ; 
plot(jacob4, jacobian(:, 4)) ; 
rat = jacobian(:, 4) ./ jacob4 ; 
median(rat) ;
% Mostly off by a factor of 0.0124 
mean(abs(rat - 0.0124) > 0.001) ; 
% but with 17% of observations 'outliers', 
% deviating from this by more than 0.001 

% manual check of jacobian(:, 1) 
res4k=CSTRfn(parvec0+[.01 0 0 0],datstruct,fitStrHot,CSTRbasis,lambda);
jacob4k = (res4k-res)/0.01 ; 
figure(7) 
plot(jacob4k, jacobian(:, 1) ) 

res4r=CSTRfn(parvec0+[0 .01 0 0],datstruct,fitStrHot,CSTRbasis,lambda);
jacob4r = (res4r-res)/0.01 ; 
figure(8) 
plot(jacob4r, jacobian(:, 2) ) 

res4a=CSTRfn(parvec0+[0 0 .01 0],datstruct,fitStrHot,CSTRbasis,lambda);
jacob4a = (res4a-res)/0.01 ; 
figure(9) 
plot(jacob4a, jacobian(:, 3) ) 

% Other parameters OK.  

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

save CSTR5 -v6; 
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
[parvec1234o, resnorm1234o, residual1234o, exitflag1234o, ... 
        output1234o, lambdaout1234o, jacobian1234o] = ...
    lsqnonlin(@CSTRfn, parvec1234, [], [], optionsCSTRfn, ...
              datstruct, fitStrH1234, CSTRbasis, lambda);
et = toc ;

save cstr4nls  -v6 parvec1234o resnorm1234o residual1234o ... 
    exitflag1234o output1234o lambdaout1234o jacobian1234o ; 
% load cstr4nls ; 
%  get the singular values of the jacobian

s1234 = svd(full(jacobian1234o));

%  log 10 of condition number

%log10(s1234(1)/s1234(3))
log10(s1234(1)/s1234(4))
% 3.0146

%  get final solution values

[resAll, DresAll, fitstructAll] = CSTRfn(parvec1234o, datstruct, ...
        fitStrH1234, CSTRbasis, lambda, 0);

%  display the parameter estimates

disp(['Initial   values: ', num2str(parvec1234)])
disp(['Estimated values: ', num2str(parvec1234o)])
disp(['True      values: ', num2str(parvectru)])

%  get the final coefficient values for the solutions

coef     = fitstructAll.coef0;
Ccoefest1234o = coef(1:nbasis);
Tcoefest1234o = coef(nbasis+1:2*nbasis);

Cfdest1234o = fd(Ccoefest1234o, CSTRbasis);
Tfdest1234o = fd(Tcoefest1234o, CSTRbasis);

Cvecest1234o = eval_fd(tobs, Cfdest1234o);
Tvecest1234o = eval_fd(tobs, Tfdest1234o);

Cvectru = eval_fd(tobs, Cfdtru);
Tvectru = eval_fd(tobs, Tfdtru);

%  display maximum absolute errors

Cerrpct1234o = 100.*median(abs(Cvecest1234o - Cvectru)./Cvectru) ;
fprintf('Cerrpct1234o = %f\n',Cerrpct1234o)
Terrpct1234o = 100.*median(abs(Tvecest1234o - Tvectru)./Tvectru) ;
fprintf('Terrpct1234o = %f\n',Terrpct1234o)

%  plot the estimated and true solutions

figure(6)
subplot(2,1,1)
plot(tobs, Cvectru, 'r-', tobs, Cvecest1234o, 'b-', tobs, Cobs, 'b.')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plot(tobs, Tvectru, 'r-', tobs, Tvecest1234o, 'b-')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

save CSTR_fig6Matlab -v6; 

% (0) To compare computation with R, starting from here:  
% load CSTR_fig6Matlab ; 
% (1) Open 'CSTR_demo.R' in R 
% (2) Create a connection in R for Matlab 
% (matlab <- Matlab())
% (3) In Matlab, addpath to 'MatlabServer'
%     at system.file("externals", package="R.matlab")
addpath('C:\Program Files\R\R-2.5.1\library\R.matlab\externals') ; 
% (4) Lock up Matlab until realsed by close(matlab) in R 
Matlab ; 
% (5) In R, open the connection 
% (isOpenMatlab <- open(matlab))
% (6) In R, transfer '.CSTRres0.trace' to Matlab
% setVariable(matlab, R_CSTRres0trace = .CSTRres0.trace)
% (7) In R, close(matlab) to allow interactive work in Matlab 
nRows = size(R_CSTRres0trace, 1) ; 
% 106 
% (8) In Matlab, add a column to R_CSTRres0trace:  
RmatlabCSTRres0trace = [R_CSTRres0trace NaN*ones(nRows, 1)] ; 
% (9) Fill the new column with CSTRfn(R_CSTRres0trace(i, 1:4), ...) 
for i = 1:nRows 
    parvec = R_CSTRres0trace(i, 1:4) ;      
    resi = CSTRfn(parvec, datstruct, fitStrH1234, CSTRbasis, lambda) ; 
    RmatlabCSTRres0trace(i, 6) = sum(resi(:).^2) ; 
end % for i = 1:nRows 

%(10) In Matlab 
MatlabServer % to transfer control back to R.
%(10) In R, reopen the connection 
% (isOpenMatlab <- open(matlab))
%(11) Get the results
% CSTRres0.R.Matlab.trace <-  getVariable(matlab, 'R_CSTRres0trace')

%%  13.  Fit the hot solution with the cool parameter estimates

%  replace true by estimated parameter values

%fitstruct.kref   = parvec(1);      
%fitstruct.EoverR = parvec(2);   
%fitstruct.a      = parvec(3);   
%fitstruct.b      = parvec(4);   

fitStrHC = fitStrHot; 
fitStrHC.kref   = parvec1234o(1);      
fitStrHC.EoverR = parvec1234o(2);   
fitStrHC.a      = parvec1234o(3);   
fitStrHC.b      = parvec1234o(4);   

% from 5.4 above
%yinitHot = [Cinit_hot, Tinit_hot];

%  load hot input into fitstruct
plot(Tcin_hot);
fitStrHC.Tcin = Tcin_hot;

%  solve  differential equation with true parameter values

h = waitbar(0,'Simulating Actual Model Output...');
[t, yestHC] = ode45(@CSTR2, tspan, yinitHot, odeoptions, ...
                            fitStrHC, 'all_hot_step', Tlim);
close(h)

%  set up separate variables for concentration and temperature

C_HC = yestHC(:,1);
T_HC = yestHC(:,2);

CerrpctHC = 100.*median(abs(C_HC - C_hot)./Cinit_hot) ;
fprintf('CerrpctHC = %f\n',CerrpctHC)
TerrpctHC = 100.*median(abs(T_HC - T_hot)./Tinit_hot) ; 
fprintf('TerrpctHC = %f\n',TerrpctHC)

% plot the estimated and true hot solutions

figure(8)
subplot(2,1,1)
lhdl = plot(t, C_hot, 'r-'); 
set(lhdl, 'LineWidth', 1);
lhdl = line(t, C_HC);
set(lhdl, 'LineWidth', 1, 'color', 'b');
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 0.1, 0.48])
subplot(2,1,2)
lhdl = plot(t, T_hot, 'r-'); 
set(lhdl, 'LineWidth', 1);
lhdl = line(t, T_HC);
set(lhdl, 'LineWidth', 1, 'color', 'b');
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 370, 420])

save CSTR_fig8Matlab -v6; 
% load CSTR_fig8Matlab ; 
% save CSTR_fig8Mat -v6 t Tlim C_hot C_HC T_hot T_HC ; 

%%  14.  Estimate kref, EoverR and a with only 
%   one response variable observed.  

% 14.1.  With only temperature observed 

fitStrHTemp = fitStrHot ; 
fitStrHTemp.fit = [0,1];

%  use coefficients from spline smoothing as initial
%  values for coefficient vectors

fitStrHTemp.coef0 = [Ccoef; Tcoef];

%  set up the true parameter values

parvectru = [kreftru, EoverRtru, atru, btru];
fitStrHTemp.b = btru;

%  calculate the estimates

tic;
[parvecHTemp, resnormHTemp, residualHTemp, exitflagHTemp] = ...
    lsqnonlin(@CSTRfn, parvec1234, [], [], optionsCSTRfn, ...
              datstruct, fitStrHTemp, CSTRbasis, lambda);
toc 

optionsCSTRfn.MaxIter = 100 ; 
tic;
[parvecHTemp, resnormHTemp, residualHTemp, exitflagHTemp] = ...
    lsqnonlin(@CSTRfn, parvec1234, [], [], optionsCSTRfn, ...
              datstruct, fitStrHTemp, CSTRbasis, lambda);
toc 

%  get final solution values

[resHTemp, DresHTemp, fitStrFinalHTemp] = CSTRfn(parvecHTemp, ...
        datstruct, fitStrHTemp, CSTRbasis, lambda, 0);

%  display the parameter estimates

disp(['Initial   values: ', num2str(parvec1234)])
disp(['Estimated values: ', num2str(parvecHTemp)])
disp(['True      values: ', num2str(parvectru)])

%  get the final coefficient values for the solutions

coefHTemp     = fitStrFinalHTemp.coef0;
CcoefestHTemp = coefHTemp(1:nbasis);
TcoefestHTemp = coefHTemp(nbasis+1:2*nbasis);

CfdestHTemp = fd(CcoefestHTemp, CSTRbasis);
TfdestHTemp = fd(TcoefestHTemp, CSTRbasis);

CvecestHTemp = eval_fd(tobs, CfdestHTemp);
TvecestHTemp = eval_fd(tobs, TfdestHTemp);

Cvectru = eval_fd(tobs, Cfdtru);
Tvectru = eval_fd(tobs, Tfdtru);

%  display maximum absolute errors

CerrpctHTemp = 100.*median(abs(CvecestHTemp - Cvectru)./Cvectru) ;

TerrpctHTemp = 100.*median(abs(TvecestHTemp - Tvectru)./Tvectru) ;
fprintf('CerrpctHTemp = %f\n',CerrpctHTemp) ; 
fprintf('TerrpctHTemp = %f\n',TerrpctHTemp) ; 

%  plot the estimated and true solutions

figure(9)
subplot(2,1,1)
plot(tobs, Cvectru, 'r-', tobs, CvecestHTemp, 'b-')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plot(tobs, Tvectru, 'r-', tobs, TvecestHTemp, 'b-')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

save CSTR_fig9Matlab -v6; 
% load CSTR_fig9Matlab -v6 ; 
save CSTR_fig9Mat tobs Tlim Cvectru CvecestHTemp Tvectru TvecestHTemp ; 
save CSTR_fig9Mat0 tobs Tlim ; 
save CSTR_fig9Mat1 Cvectru Tvectru ; 
save CSTR_fig9Mat2 CvecestHTemp TvecestHTemp ; 

csvwrite('fit9Mat.csv', [Tobs Cvectru CvecestHTemp ... 
    Tvectru TvecestHTemp]) ; 


% 14.2.  With only concentration observed 

fitStrHConc = fitStrHot ; 
fitStrHConc.fit = [1, 0];

%  use coefficients from spline smoothing as initial
%  values for coefficient vectors

fitStrHConc.coef0 = [Ccoef; Tcoef];

%  set up the true parameter values

parvectru = [kreftru, EoverRtru, atru, btru];
fitStrHConc.b = btru;

%  calculate the estimates

tic;
[parvecHConc, resnormHConc, residualHConc, exitflagHConc] = ...
    lsqnonlin(@CSTRfn, parvec1234, [], [], optionsCSTRfn, ...
              datstruct, fitStrHConc, CSTRbasis, lambda);
toc

%  get final solution values

[resHConc, DresHConc, fitStrFinalHConc] = CSTRfn(parvecHConc, datstruct, ... 
                fitStrHConc, CSTRbasis, lambda, 0);

%  display the parameter estimates

disp(['Initial   values: ', num2str(parvec1234)])
disp(['Estimated values: ', num2str(parvecHConc)])
disp(['True      values: ', num2str(parvectru)])

%  get the final coefficient values for the solutions

coefHConc     = fitStrFinalHConc.coef0;
CcoefestHConc = coefHConc(1:nbasis);
TcoefestHConc = coefHConc(nbasis+1:2*nbasis);

CfdestHConc = fd(CcoefestHConc, CSTRbasis);
TfdestHConc = fd(TcoefestHConc, CSTRbasis);

CvecestHConc = eval_fd(tobs, CfdestHConc);
TvecestHConc = eval_fd(tobs, TfdestHConc);

Cvectru = eval_fd(tobs, Cfdtru);
Tvectru = eval_fd(tobs, Tfdtru);

%  display maximum absolute errors

CerrpctHConc = 100.*median(abs(CvecestHConc - Cvectru)./Cvectru) ;
TerrpctHConc = 100.*median(abs(TvecestHConc - Tvectru)./Tvectru) ;
fprintf('CerrpctHConc = %f\n',CerrpctHConc) ; 
fprintf('TerrpctHConc = %f\n',TerrpctHConc) ; 

%  plot the estimated and true solutions

figure(10) % was figure(9): second copy for obs. Conc only 
subplot(2,1,1)
plot(tobs, Cvectru, 'r-', tobs, CvecestHConc, 'b-')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 1.2, 1.8])
subplot(2,1,2)
plot(tobs, Tvectru, 'r-', tobs, TvecestHConc, 'b-')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 330, 360])

save CSTR_fig10Matlab -v6; 
% load CSTR_fig10Matlab ; 

%%  15.  Fit the hot solution with the cool parameter estimates

%  replace true by estimated parameter values

%fitstruct.kref   = parvec(1);      
%fitstruct.EoverR = parvec(2);   
%fitstruct.a      = parvec(3);   
%fitstruct.b      = parvec(4);   

fitStrHTC = fitStrHTemp; 
fitStrHTC.kref   = parvecHConc(1);      
fitStrHTC.EoverR = parvecHConc(2);   
fitStrHTC.a      = parvecHConc(3);   
fitStrHTC.b      = parvecHConc(4);   

% from 5.4 above
%yinitHot = [Cinit_hot, Tinit_hot];

%  load hot input into fitstruct

fitStrHTC.Tcin = Tcin_hot;

%  solve  differential equation with true parameter values

h = waitbar(0,'Simulating Actual Model Output...');
[t, yestHTC] = ode45(@CSTR2, tspan, yinitHot, odeoptions, ...
                            fitStrHTC, 'all_hot_step', Tlim);
close(h)

%  set up separate variables for concentration and temperature

C.HTC= yestHTC(:,1);
T.HTC = yestHTC(:,2);

% plot the estimated and true hot solutions

figure(11) % was figure(10) 
subplot(2,1,1)
plot(t, C_hot, 'r-', t, C.HTC, 'b-')
ylabel('\fontsize{16} C(t)')
title('\fontsize{16} Concentration (red = true, blue = estimated)')
axis([0, Tlim, 0, 0.5])
subplot(2,1,2)
plot(t, T_hot, 'r-', t, T.HTC, 'b-')
ylabel('\fontsize{16} T(t)')
title('\fontsize{16} Temperature')
axis([0, Tlim, 370, 420])

save CSTR_fig11Matlab -v6; 

