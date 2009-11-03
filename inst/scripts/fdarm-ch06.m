%
% Ramsay, Hooker & Graves (2009)
% Functional Data Analysis with R and Matlab (Springer)

%  Remarks and disclaimers

%  These R commands are either those in this book, or designed to 
%  otherwise illustrate how R can be used in the analysis of functional
%  data.  
%  We do not claim to reproduce the results in the book exactly by these 
%  commands for various reasons, including:
%    -- the analyses used to produce the book may not have been
%       entirely correct, possibly due to coding and accuracy issues
%       in the functions themselves 
%    -- we may have changed our minds about how these analyses should be 
%       done since, and we want to suggest better ways
%    -- the R language changes with each release of the base system, and
%       certainly the functional data analysis functions change as well
%    -- we might choose to offer new analyses from time to time by 
%       augmenting those in the book
%    -- many illustrations in the book were produced using Matlab, which
%       inevitably can imply slightly different results and graphical
%       displays
%    -- we may have changed our minds about variable names.  For example,
%       we now prefer "yearRng" to "dayrange" for the weather data.
%    -- three of us wrote the book, and the person preparing these scripts
%       might not be the person who wrote the text
%  Moreover, we expect to augment and modify these command scripts from time
%  to time as we get new data illustrating new things, add functionality
%  to the package, or just for fun.

%
% ch. 6.  Descriptions of Functional Data
%

%  Set up some strings for constructing paths to folders.
%  These strings should be modified so as to provided access
%  to the specified folders on your computer.

%  Path to the folder containing the Matlab functional data analysis
%  software

fdaMPath = 'c:/Program Files/MATLAB/R2009a/fdaM';

addpath(fdaMPath)

%  Path to the folder containing the examples

examplesPath = [fdaMPath,'/examples'];

addpath(examplesPath)

%
% Section 6.1 Some Functional Descriptive Statistics
%

%  set up a functional data object for log precipitation

%   ----------  Statistics for the log precipitation data  ---------------

%  path to the daily weather data example folder

weatherPath = [examplesPath,'/weather'];

addpath(weatherPath)

load daily

tempav  = daily.tempav;
precav  = daily.precav;
yearRng = daily.rng;
time    = daily.time;
place   = daily.place;
dayperiod = daily.period;

%  define the harmonic acceleration operator

Lcoef   = [0,(2*pi/dayperiod)^2,0];    
harmaccelLfd = vec2Lfd(Lcoef, yearRng); 

%  organize data to have winter in the center of the plot

dayOfYearShifted = [182:365, 1:181];

%  change 0's to 0.05 mm in precipitation data

prectmp = precav;
for j=1:35
    index = prectmp(:,j)==0;
    prectmp(index,j) = 0.05;
end

%  set up functional data object for log precipitation

logprecav = log10(prectmp(dayOfYearShifted,:));

%  set up a saturated basis: as many basis functions as observations

daybasis  = create_fourier_basis(yearRng, 365);

%  smooth data with lambda that minimizes GCV

lambda      = 1e6;
fdParobj    = fdPar(daybasis, harmaccelLfd, lambda);
logprec_fd = smooth_basis(time, logprecav, fdParobj);
fdnames{1} = 'Day (July 1 to June 30)';
fdnames{2,1} = 'Weather Station';
fdnames{2,2} = place;
fdnames{3} = 'Log 10 Precipitation (mm)';
logprec_fd  = putnames(logprec_fd, fdnames);

%  elementary pointwise mean and standard deviation

meanlogprec   = mean(logprec_fd);
stddevlogprec = std_fd(logprec_fd);

% Section 6.1.1 The Bivariate Covariance Function v(s; t)

logprecvar_bifd = var_fd(logprec_fd);

weektime        = linspace(0,365,53);
logprecvar_mat  = eval_bifd(weektime, weektime, logprecvar_bifd);

% Figure 6.1

surf(weektime, weektime, logprecvar_mat)
xlabel('\fontsize{13} Day (July 1 to June 30)')
ylabel('\fontsize{13} Day (July 1 to June 30)')
zlabel('\fontsize{13} variance(log10 precip)')

contour(weektime, weektime, logprecvar_mat);
xlabel('\fontsize{13} Day (July 1 to June 30)')
ylabel('\fontsize{13} Day (July 1 to June 30)')

% Figure 6.2

%  Matlab doesn't indicate contour level values in the plot,
%  but instead uses an optional colorbar to indicate height.

contour(weektime, weektime, logprecvar_mat);
xlabel('\fontsize{13} Day (July 1 to June 30)')
ylabel('\fontsize{13} Day (July 1 to June 30)')
colorbar

%
% Section 6.2 The Residual Variance-Covariance Matrix Se
%
%  (no computations in this section)

%
% Section 6.3 Functional Probes rho[xi]
%

% see section 6.5 below

%
% Section 6.4 Phase-plane Plots of Periodic Effects
%

%  -----------  Phase-plane plots for the nondurable goods data  ---------

%  path to the nondurable goods index example folder

goodsindexPath = [examplesPath,'/goodsindex'];

addpath(goodsindexPath)

%  load the goodsindex data from the .mat file goodsindex.mat

load goodsindex

ndur        = goodsindex.ndur;
yeartime    = goodsindex.yeartime;
nondurables = goodsindex.nondurables;
lognondur   = goodsindex.lognondur;

%  set up a basis for smoothing log nondurable goods index

goodsRng = [1919,2000];
nbasis = 979;
norder = 8;
goodsbasis = create_bspline_basis(goodsRng, nbasis, norder);

%  set up fdPar object

LfdobjNonDur= int2Lfd(4);
lambda = 1e-11;
goodsLfdPar = fdPar(goodsbasis, LfdobjNonDur, lambda);

logNondur_fd = smooth_basis(yeartime, lognondur, goodsLfdPar);

% Fig. 6.3 The log nondurable goods index for 1964 to 1967

t64_67 = linspace(1964, 1967, 601);
logNondur_vec = eval_fd(t64_67, logNondur_fd);

sel64_67 = find(1964 <= yeartime & yeartime <= 1967);
ymin = min(lognondur(sel64_67));
ymax = max(lognondur(sel64_67));
phdl=plot(yeartime(sel64_67), lognondur(sel64_67), 'bo', ...
          t64_67, logNondur_vec, 'b-', ...
     [1965,1965], [ymin,ymax], 'b--', [1966,1966], [ymin,ymax], 'b--');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Year')
ylabel('\fontsize{13} Log_{10} Nondurable Goods Index')
axis([1964,1967,1.615,1.725])

% Section 6.4.1.  Phase-plane Plots Show Energy Transfer

% Figure 6.4.  Phase-plane plot for a simple harmonic function

tvec     = linspace(0,1,101);
sinvec   = sin(2*pi*tvec);
cosvec   = cos(2*pi*tvec);
Dsinvec  = 2*pi*cosvec;
D2sinvec = -(2*pi).^2*sinvec;

phdl=plot(Dsinvec, D2sinvec, 'b-', ...
          [min(Dsinvec),max(Dsinvec)],   [0,0], 'b:', ...
          [0,0], [min(D2sinvec),max(D2sinvec)], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Velocity')
ylabel('\fontsize{13} Acceleration')
text(-4,  45, '\fontsize{12} Max. potential energy')
text(-4, -45, '\fontsize{12} Max. potential energy')
text(-9.5,  5, '\fontsize{12} Max.')
text(-9.5,  0, '\fontsize{12} kinetic')
text(-9.5, -5, '\fontsize{12} energy')
text( 6.5,  5, '\fontsize{12} Max.')
text( 6.5,  0, '\fontsize{12} kinetic')
text( 6.5, -5, '\fontsize{12} energy')
axis([-10,10,-50,50])

% Section 6.4.2 The Nondurable Goods Cycles

%  -----------  Phase-plane plots for the nondurable goods data  ---------

% Figure 6.5

startyr   = 64;   %   starting year
yearindex = startyr + 1900 - 1919;  % indices of year

%  Set up the derivatives for phase-plane plotting.

index    = (1:(13)) + yearindex*12;
durrange = [min(yeartime(index)), max(yeartime(index))];
durfine  = linspace(durrange(1),durrange(2),401)';
durcrse  = linspace(durrange(1),durrange(2),13)';

D1f = eval_fd(durfine, logNondur_fd, 1);
D2f = eval_fd(durfine, logNondur_fd, 2);
D1c = eval_fd(durcrse, logNondur_fd, 1);
D2c = eval_fd(durcrse, logNondur_fd, 2);

%  The phase-plane plotting.

monthlabs = ['jan';'F  ';'mar';'A  ';'M  ';'Jun'; ...
             'Jly';'Aug';'S  ';'Oct';'N  ';'D  '];

plot([0,0], [-12,12], 'b:', [-.75,.75], [0,0], 'b:')
x0 = durrange(1);
hold on
indj = find(durfine >= x0 & durfine <= x0+1);
phdl  = plot(D1f(indj), D2f(indj), 'b-');
set(phdl, 'LineWidth', 2);
indexc =  1:12; 
for k=1:12
    plot(D1c(indexc(k)), D2c(indexc(k)), 'o')
    text(D1c(indexc(k))+0.05, D2c(indexc(k))+0.5, monthlabs(k,:))
end
hold off
axis([-.5, .5, -12, 12])
xlabel('\fontsize{13} Velocity') 
ylabel('\fontsize{13} Acceleration')

% sec. 6.4.3.  Phase-Plane Plotting the Growth of Girls

%  -------------  Phase-plane diagrams for the growth data  ---------------

%  path to the Berkeley growth data example folder

growthPath = [examplesPath,'/growth'];

addpath(growthPath)

%  load the growth data struct object from the .mat file growth.mat

load growth

%  the dimensions of the data
ncasem = growth.ncasem;
ncasef = growth.ncasef;
nage   = growth.nage;
%  the heights of the 54 girls
heightmat = growth.hgtfmat;
%  the 31 ages of measurement
age = growth.age;

%  define the range of the ages and set up a fine mesh of ages

ageRng  = [1,18];

norder = 6;
breaks = age;
nbasis = length(breaks) + norder - 2;
gr_basis = create_bspline_basis(ageRng, nbasis, norder, breaks);

children = 1:10;
coef     = zeros(nbasis,10);
Wfd0     = fd(coef, gr_basis);
gr_Lfd    = 3;
gr_lambda = 10^(-1.5);
gr_fdPar1_5 = fdPar(Wfd0, gr_Lfd, gr_lambda);

[Wfd, beta] = smooth_monotone(age, heightmat(:,children), gr_fdPar1_5);

agefine = linspace(1,18,101);
i11_7  = find(abs(agefine-11.7) == min(abs(agefine-11.7)));

betamat = repmat(beta(2,:),101,1);

velffine = betamat.*eval_monfd(Wfd, agefine, 1);
accffine = betamat.*eval_monfd(Wfd, agefine, 2);

phdl=plot(velffine(:,1),      accffine(:,1),      'b-', ...
          velffine(i11_7, 1), accffine(i11_7, 1), 'bo', ...
          [0,12], [0,0], 'b:');
set(phdl, 'LineWidth', 2);
hold on
for i = 2:10
    phdl=plot(velffine(:,i),      accffine(:,i),      'b-', ...
        velffine(i11_7, i), accffine(i11_7, i), 'bo');
    set(phdl, 'LineWidth', 2);
end
hold off
xlabel('\fontsize{13} Velocity (cm/yr)')
ylabel('\fontsize{13} Acceleration (cm/yr^2)')
axis([0,12,-6,3])

%
% Section 6.5 Confidence Intervals for Curves and their Derivatives
%

% sec. 6.5.1.  Two Linear mappings Defining a Probe Value

%  -----------------  Probe values for Canadian weather data  ------------

%  define a probe to emphasize mid-winter

dayvec  = linspace(0,365,101);
xivec   = exp(20*cos(2*pi*(dayvec-197)/365));
xibasis = create_bspline_basis(yearRng,13);
xifd    = smooth_basis(dayvec, xivec, xibasis);

plot(xifd)

%  define bases for temperature and precipitation

tempbasis = create_fourier_basis(yearRng, 65);
precbasis = create_fourier_basis(yearRng,365);

%  probe values for basis functions with respect to xifd

tempLmat = inprod(tempbasis, xifd);
precLmat = inprod(precbasis, xifd);

% sec. 6.5.3.  Confidence Limits for Prince Rupert's Log Precipitation

%  smooth data with lambda that minimizes GCV getting
%  all of the output up to matrix y2cMap

lambda      = 1e6;
fdParobj    = fdPar(daybasis, harmaccelLfd, lambda);

[logprec_fd, df, gcv, SSE, penmat, y2cMap] = ...
              smooth_basis(time, logprecav, fdParobj);
fdnames{1}   = 'Day (July 1 to June 30)';
fdnames{2,1} = 'Weather Station';
fdnames{2,2} = place;
fdnames{3}   = 'Log 10 Precipitation (mm)';
logprec_fd   = putnames(logprec_fd, fdnames);

%  compute the residual matrix and variance vector

logprecmat = eval_fd(time, logprec_fd);
logprecres = logprecav - logprecmat;
logprecvar = sum(logprecres.^2, 2)/(35-1);

%  smooth log variance vector

lambda      = 1e8;
resfdParobj = fdPar(daybasis, harmaccelLfd, lambda);

logvar_fd = smooth_basis(time, log(logprecvar), resfdParobj);

%  evaluate the exponentiated log variance vector and
%  set up diagonal error variance matrix SigmaE

varvec      = exp(eval_fd(time, logvar_fd));
SigmaE      = diag(varvec);

%  compute variance covariance matrix for fit

c2rMap        = eval_basis(time, daybasis);
Sigmayhat     = c2rMap * y2cMap * SigmaE * y2cMap' * c2rMap';

%  extract standard error function for yhat

logprec_stderr= sqrt(diag(Sigmayhat));

%  plot Figure 6.6

logprecvec29 = eval_fd(time,logprec_fd(29));

phdl = plot(time, logprecvec29,                 'b-', ...
            time, logprec29 + 2*logprec_stderr, 'b--', ...
            time, logprec29 - 2*logprec_stderr, 'b--', ...
            time, logprecav(:,29),              'bo');
set(phdl, 'LineWidth', 2);
xlabel('\fontsize{13} Day (July 1 to June 30)')
ylabel('\fontsize{13} Log 10 Precipitation (mm)')
title('\fontsize{13} Prince Rupert')
axis([0,365,0.2,1.3])

%
% Section 6.6 Some Things to Try
%
% (exercises for the reader)
