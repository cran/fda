%
%
% Ramsay, Hooker & Graves (2009)
% Functional Data Analysis with R and Matlab (Springer)
%

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
% ch. 5.  Smoothing: Computing Curves from Noisy Data
%

%  Set up some strings for constructing paths to folders.
%  These strings should be modified so as to provided access
%  to the specified folders on your computer.

%  Path to the folder containing the Matlab functional data analysis
%  software

fdaMPath = 'c:/Program Files/MATLAB/R2009a/fdaM';

addpath(fdaMPath)

%  -------------------  Smoothing the growth data  ------------------------

%  Path to the folder containing the examples

examplesPath = [fdaMPath,'/examples'];

addpath(examplesPath)

%
% Section 5.1.  Regression Splines: Smoothing by Regression Analysis
%

%  path to the Berkeley growth data example folder

growthPath = [examplesPath,'/growth'];

addpath(growthPath)

%  load the growth data struct object from the .mat file growth.mat

load growth

%  the dimensions of the data
ncasem = growth.ncasem;
nage   = growth.nage;
%  the heights of the 54 girls
heightmat = growth.hgtfmat;
%  the 31 ages of measurement
age = growth.age;

%  define the range of the ages and set up a fine mesh of ages

ageRng  = [1,18];
agefine = linspace(1, 18, 501);

%  set up basis with 12 basis functions

nbasis     = 12;
norder     = 6;
heightbasis12 = create_bspline_basis(ageRng, nbasis, norder);

%  fit the data by least squares

basismat     = eval_basis(age, heightbasis12);
heightcoef   = basismat\heightmat;
heighthatmat = basismat*heightcoef;

%  fit the data using function smooth_basis, which does the same thing.

[height_fd, height_df, height_gcv] = ...
              smooth_basis(age, heightmat, heightbasis12);
          
%  compute

y2cMap = (heightbasis12'*heightbasis12)\heightbasis12';

%
% Section 5.2.  Data Smoothing with Roughness Penalties
%

% section 5.2.2 The Roughness Penalty Matrix R

%  ---------------  Smoothing the Canadian weather data  ------------------

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

%  define a basis for daily temperature data

nbasis    = 65;
tempbasis = create_fourier_basis(yearRng,nbasis);

%  define the harmonic acceleration operator

%  When the coefficients of the linear differential operator
%  are constant, the Lfd object can be set up more simply
%  by using function vec2Lfd as follows:

%  The first argument is a vector of coefficients for the
%  operator, and the second argument is the range over which
%  the operator is defined.

harmaccelLfd = vec2Lfd([0,(2*pi/dayperiod)^2,0], yearRng);

%  compute the penalty matrix R

Rmat = eval_penalty(tempbasis, harmaccelLfd);

% section 5.2.4 Defining Smoothing by Functional Parameter Objects

%  -------  Smoothing the growth data with a roughness penalty  -----------

%  set up a basis for the growth data
%  with knots at ages of height measurement

norder      = 6;
nbasis      = length(age) + norder - 2;
heightbasis = create_bspline_basis(ageRng, nbasis, norder, age);

%  define a functional parameter object for smoothing

heightLfd    = 4;
heightlambda = 0.01;
heightfdPar  = fdPar(heightbasis, heightLfd, heightlambda);

%  smooth the data

height_fd = smooth_basis(age, heightmat, heightfdPar);

% section 5.2.5 Choosing Smoothing Parameter lambda by minimizing GCV

loglam         = -6:0.25:0;
nloglam        = length(loglam);
Gcvsave        = zeros(nloglam,1);
Dfsave         = Gcvsave;
for i = 1:nloglam
    heightLfd    = 4;
    heightLambda = 10^loglam(i);
    heightdPari  = fdPar(heightbasis, heightLfd, heightLambda);
    [hgtfdi, hgtdfi, hgtgcvi] = smooth_basis(age, heightmat, heightdPari);
    Gcvsave(i) = sum(hgtgcvi);
    Dfsave(i)  = hgtdfi;
end

% Figure 5.1.

phdl=plot(loglam, Gcvsave, 'o-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log_{10}\lambda')
ylabel('\fontsize{13} GCV(\lambda)')

%
% 5.3.  Case Study: The Log Precipitation Data
%

%  organize data to have winter in the center of the plot

dayOfYearShifted = [182:365, 1:181];

%  change 0's to 0.05 mm in precipitation data

prectmp = precav;
for j=1:35
    index = find(prectmp(:,j)==0);
    prectmp(index,j) = 0.05;
end

%  set up functional data object for log precipitation

logprecav = log10(prectmp(dayOfYearShifted,:));

%  set up a saturated basis: as many basis functions as observations

daybasis  = create_fourier_basis(yearRng, 365);

%  see above for setting the harmonic acceleration operator

%  step through values of log(lambda)

loglam         = 4:0.25:9;
nloglam        = length(loglam);
Gcvsave        = zeros(nloglam,1);
Dfsave         = Gcvsave;
for ilam = 1:nloglam
    disp(['log10 lambda = ',num2str(loglam(ilam))])
    lambda        = 10^loglam(ilam);
    fdParobj      = fdPar(daybasis, harmaccelLfd, lambda);
    [logprec_fd, logprec_df, logprec_gcv] = ...
                    smooth_basis(time, logprecav, fdParobj);
    Dfsave(ilam)  = logprec_df;
    Gcvsave(ilam) = sum(logprec_gcv);
end

% Figure 5.2.

phdl=plot(loglam, Gcvsave, 'bo-'); 
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log_{10}\lambda')
ylabel('\fontsize{13} GCV(\lambda)')

%  smooth data with minimizing value of lambda

lambda      = 1e6;
fdParobj    = fdPar(daybasis, harmaccelLfd, lambda);
logprec_fd = smooth_basis(time, logprecav, fdParobj);
fdnames{1} = 'Day (July 1 to June 30)';
fdnames{2,1} = 'Weather Station';
fdnames{2,2} = place;
fdnames{3} = 'Log 10 Precipitation (mm)';
logprec_fd  = putnames(logprec_fd, fdnames);

%  plot the functional data object

plot(logprec_fd)

% plotfit_fd:  Pauses between plots
% *** --->>> input required press any key to advance to the next plot

plotfit_fd(logprecav, time, logprec_fd)

%
% Section 5.4 Positive, Monotone, Density
%             and Other Constrained Functions
%

%   ----------------  Positive smoothing of precipitation  ----------------

lambda      = 1e3;
WfdParobj   = fdPar(daybasis, harmaccelLfd, lambda);
VanPrec     = precav(dayOfYearShifted, 26);
VanWfd      = smooth_pos(time, VanPrec, WfdParobj);

VanPrecFit = exp(eval_fd(time, VanWfd));

phdl = plot(time, VanPrec, 'o', time, VanPrecFit, 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day (July 1 to June 30)')
ylabel('\fontsize{13} Millimeters')
title('\fontsize{13} Vancouver Precipitation');

% 5.4.2.  Monotone smoothing

%  ------------  5.4.2.1  Monotone smoothing  of the tibia data  ----------

%  path to the growth of the infant's tibia example folder

growthPath = [examplesPath,'/growth'];

%  load the infant's data from file infantGrowth.mat

load infantGrowth

%  set up the data for analysis

day = infantGrowth.day;
tib = infantGrowth.length;
n   = length(tib);

%  a basis for monotone smoothing

nbasis = 42;
Wbasis   = create_bspline_basis([1,n], nbasis);

%  the fdPar object for smoothing

infantLfd    = 2;
infantlambda = 1e-4;
WfdPar       = fdPar(Wbasis, infantLfd, infantlambda);

%  smooth the data

[Wfd, beta] = smooth_monotone(day, tib, WfdPar);

%  compute fit and derivatives of fit

dayfine  = linspace(1,n,151);
tibhat   = beta(1)+beta(2)*eval_monfd(dayfine ,Wfd);
Dtibhat  =         beta(2)*eval_monfd(dayfine, Wfd, 1);
D2tibhat =         beta(2)*eval_monfd(dayfine, Wfd, 2);

%  plot height

subplot(1,2,1)
phdl=plot(day, tib, 'bo', dayfine, tibhat, 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Tibia Length (mm)')
axis('square')

%  plot velocity

subplot(1,2,2)
phdl=plot(dayfine, Dtibhat, 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Tibia Velocity (mm/day)')
axis('square')

%  plot acceleration

subplot(1,1,1)
phdl=plot(dayfine, D2tibhat, 'b-', [1,n], [0,0], 'b--');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Tibia Acceleration (mm/day/day)')

% ---------  5.4.2.2  Monotone smoothing the Berkeley female data  --------

% 
%  Compute the monotone smoothing of the Berkeley female growth data.
%

%  the growth data are set up above

%  the data

height = heightmat;
ncasef = size(height,2);

%  an order 6 bspline basis with knots at ages of measurement

norder = 6;
nbasis = nage + norder - 2;
wbasis = create_bspline_basis(ageRng, nbasis, norder, age);

%  define the roughness penalty for function W

Lfdobj    = 3;          %  penalize curvature of acceleration
lambda    = 10^(-0.5);  %  smoothing parameter
cvecf     = zeros(nbasis, ncasef);
Wfd0      = fd(cvecf, wbasis);
growfdPar = fdPar(Wfd0, Lfdobj, lambda);

%  monotone smoothing

[Wfd, betaf, heighthatfd] = smooth_monotone(age, height, growfdPar);

%  Set up functional data objects for the acceleration curves 
%  and their mean.  Suffix UN means 'unregistered'.

accelfdUN     = deriv_fd(heighthatfd, 2);
accelmeanfdUN = mean(accelfdUN);

%  plot unregistered curves

subplot(1,1,1)
plot(accelfdUN)
xlabel('\fontsize{13} Age') 
ylabel('\fontsize{13} Acceleration (cm/yr/yr)')
axis([1,18,-6,4])

% 5.4.3.  Probability density function for Regina precipitaton

%  -----------  Density function for Regina Precipitation  ----------------

%  path to the daily weather data example folder

weatherPath = [examplesPath,'/weather'];

addpath(weatherPath)

%  load data from file ReginaPrecip.mat

load ReginaPrecip

%  plot the empirical quantile function

NR = length(ReginaPrecip);
plot(1:NR, sort(ReginaPrecip), 'b-')
xlabel('\fontsize{13} Rank of rainfall')
ylabel('\fontsize{13} Ordered daily rainfall (mm)' )

%  precipitation has an extremely long tail.
%  here we only use precipitations between 2 and 45 mm/day.

sel2_45 = find(2 <= ReginaPrecip & ReginaPrecip <= 45);
RegPrec = sort(ReginaPrecip(sel2_45));
N       = length(RegPrec);
precRng = [2,45];

%  set up spline basis for log precipitation density
%  with knots logarithmicaly spaced

breaks     = 10.^linspace(log10(2), log10(45), 11);
breaks(1)  = precRng(1);
breaks(11) = precRng(2);
Wnbasis    = length(breaks) + 2;
Wbasis     = create_bspline_basis(precRng, Wnbasis, 4, breaks);

%  set up the functional parameter object

precLfd = 2;
Wlambda = 1e-1;
WfdPar  = fdPar(Wbasis, precLfd, Wlambda);

%  estimate the density

[Wfd, C] = density_fd(RegPrec, WfdPar);

%  plot the density with the knot locations

Zfine = linspace(2,45,201);
Wfine = eval_fd(Zfine, Wfd);
Pfine = exp(Wfine)/C;
Plim  = max(Pfine);

phdl=plot(Zfine, Pfine, 'b-'); 
set(phdl, 'LineWidth', 2)
hold on
for i=1:9
    plot([breaks(i+1),breaks(i+1)], [0,Plim], 'b--')
end
hold off
xlabel('\fontsize{13} Precipitation (mm)')
ylabel('\fontsize{13} Probability Density')

%
% Section 5.5 Assessing the Fit to the Log Precipitation Data
%

%  logprec_fd was set up above

logprecmat = eval_fd(time, logprec_fd);
logprecres = logprecav - logprecmat;

%  variance across stations

logprecvar1 = sum(logprecres.^2, 2)/35;

%  variance across time

logprecvar2 = sum(logprecres.^2, 1)/(365-12);

% Figure 5.7

%  indices for 'Winnipeg', 'Regina', 'Churchill', 'Montreal', 'St. Johns'
rt    = [1,12,17,19,20];
%  indices for 'Yellowknife', 'Resolute', 'Vancouver', 'Iqaluit',
%              'Pr. George', 'Pr. Rupert'
lft   = [26,28,29,32,33,35];
%  index for 'Edmonton'
below = 23;
%  index for 'Halifax'
top   = 2;

%  plot with labels on selected points.  This is easier in R!

phdl=plot((1:35), sqrt(logprecvar2), 'bo');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Station Number')
ylabel('\fontsize{13} Standard Deviation across Day')
for i=1:5
    text(rt(i)+0.5, sqrt(logprecvar2(rt(i))), ...
         ['\fontsize{12} ',place(rt(i),:)])
end
for i=1:5
    text(lft(i)-7,  sqrt(logprecvar2(lft(i))), ...  
        ['\fontsize{12} ',place(lft(i),:)])
end
i=1;
    text(below(i)-3.5, sqrt(logprecvar2(below(i)))-0.005,  ...
        ['\fontsize{12} ',place(below(i),:)])
i=1;
    text(top(i)-2,   sqrt(logprecvar2(top(i)))+0.005,  ... 
        ['\fontsize{12} ',place(top(i),:)])
axis([0,35,0.11,0.26])

% Figure 5.8

logstddev_fd = smooth_basis(time, log(logprecvar1)/2, fdParobj);
logprecvar1fit = exp(eval_fd(time, logstddev_fd));

phdl=plot(time, sqrt(logprecvar1), 'bo', time, logprecvar1fit, 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day')
ylabel('\fontsize{13} Standard seviation across stations')
axis([0,365,0.1,0.35])

%
% Section 5.6 Details for the fdPar Class and smooth.basis Function
%

help(fdPar)
help(smooth_basis)

%
% Section 5.8 Some Things to Try
%
% (exercises for the reader)

%
% Section 5.7 More to Read
%
