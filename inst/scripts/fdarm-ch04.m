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
% ch. 4  How to Build Functional Data Objects
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
% Section 4.1 Adding Coefficients to Bases to Define Functions
%

%  4.1.1 Coefficient Vectors, Matrices and Arrays

%  path to the daily weather data example folder

weatherPath = [examplesPath,'/weather'];

addpath(weatherPath)

load daily

tempav  = daily.tempav;
precav  = daily.precav;
yearRng = daily.rng;
time    = daily.time;
place   = daily.place;

nbasis    = 65;
tempbasis = create_fourier_basis(yearRng, nbasis);
coefmat   = zeros(nbasis, 35);
tempfd    = fd(coefmat, tempbasis);

% 4.1.2 Labels for Functional Data Objects

%  the growth data

growth_fdnames{1} = 'Age (years)';
growth_fdnames{2} = 'Child';
growth_fdnames{3} = 'Height (cm)';

%  the daily temperature data with names for stations

station{1} = 'Weather Station';
station{2} = place;

temp_fdnames{1} = 'Day'; 
temp_fdnames{2} = station;
temp_fdnames{3} = 'Mean temperature (deg C)';

%
% 4.2 Methods for Functional Data Objects
%

%  Two order 2 splines over unit interval

unitRng = [0,1];
nbasis = 2;
norder = 2;
bspl2  = create_bspline_basis(unitRng, nbasis, norder);
plot(bspl2)

%  a pair of straight lines

tstFn1 = fd([-1; 2], bspl2);
tstFn2 = fd([ 1; 3], bspl2);

subplot(3,1,1)
plot(tstFn1)
xlabel('')
ylabel('Line 1')
subplot(3,1,2)
plot(tstFn2)
xlabel('')
ylabel('Line 2')

%  sum of lines

fdsumobj = tstFn1 + tstFn2;
subplot(3,1,3)
plot(fdsumobj)
xlabel('')
ylabel('Line 1 + Line 2')

%  difference between lines

fddifobj = tstFn2 - tstFn1;
subplot(3,1,3)
plot(fddifobj)
xlabel('')
ylabel('Line 2 - Line 1')

%  product of two straight lines

fdprdobj = tstFn1 .* tstFn2;
subplot(3,1,3)
plot(fdprdobj)
xlabel('')
ylabel('Line 1 * Line 2')

%  square a straight line

fdsqrobj = tstFn1.^2;
subplot(3,1,3)
plot(fdsqrobj)
xlabel('')
ylabel('Line 1 ^2')

%  square root of a line with negative values:  illegal

a    = 0.5;
fdrootobj = tstFn0.^a;
subplot(3,1,3)
plot(fdrootobj)

%  square root of a square:  this illustrates the hazards of
%  fractional powers when values are near zero.  The right answer is
%  two straight line segments with a discontinuity in the first
%  derivative.  It would be better to use order two splines and
%  put a knot at the point of discontinuity, but the power method
%  doesn't know how to do this.

fdrootobj = fdsqrobj.^a;
subplot(3,1,3)
plot(fdrootobj)

%  square root of a quadratic without values near zero:  no problem

fdrootobj = (fdsqrobj + 1).^a;
subplot(3,1,3)
plot(fdrootobj)

%  reciprocal of a function with zero values:  illegal operation

a    = (-1);
fdinvobj = tstFn0.^a;
subplot(3,1,3)
plot(fdinvobj)

%  reciprocal of a function with near zero values:  a foolish thing
%  to do and the power function fails miserably

fdinvobj = fdsqrobj.^a;
subplot(3,1,3)
plot(fdinvobj)

%  reciprocal of a positive function with no values near zero

fdinvobj = (fdsqrobj+1).^a;
subplot(3,1,3)
plot(fdinvobj)

%  near reciprocal of a positive function with no values near zero

a = -0.99;
fdpowobj = (fdsqrobj+1).^a;
subplot(3,1,3)
plot(fdpowobj)

% compute mean temperature in two ways and plot the difference

Tempfd     = smooth_basis(time, tempav, tempbasis);
meanTempfd = mean(Tempfd);
sumTempfd  = sum(Tempfd);
subplot(1,1,1)
plot((meanTempfd-sumTempfd.*(1/35)));

%  plot the temperature for Resolute and add the Canadian mean 

plot(Tempfd(35))
line(meanTempfd)
axis([0,365,-35,20])

%  evaluate the derivative of mean temperature and plot

DmeanTempVec = eval_fd(day5, meanTempfd, 1);
plot(day5, DmeanTempVec)

%  evaluate and plot the harmonic acceleration of mean temperature

Lbasis  = create_constant_basis(dayrange);  %  create a constant basis
Lcoef   = [0,(2*pi/dayperiod)^2,0];    %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

LmeanTempVec = eval_fd(day5, meanTempfd, harmaccelLfd);

plot(day5, LmeanTempVec, 'b-', yearRng, [0,0], 'b--')
xlabel('\fontsize{13} Day')
ylabel('\fontsize{13} Harmonic Acceleration')
axis([0,365,-0.05,0.05])

% Figure 4.1.

%  plot temperatures with winter in the center of the plot

dayOfYearShifted = [182:365, 1:181];

temp_fd = smooth_basis(time, tempav(dayOfYearShifted,:), tempbasis);

plot(temp_fd) 
xlabel('\fontsize{13} Day (July 1 to June 30)')
ylabel('\fontsize{13} Mean temperature (deg. C)')

% Section 4.2.1 Illustration: Sinusoidal Coefficients

% Figure 4.2

tenRng   = [0,10];
basis13  = create_bspline_basis(tenRng, 13);
tvec     = linspace(0,1,13)';
sinecoef = sin(2*pi*tvec);
sinefd   = fd(sinecoef, basis13);
tfine    = linspace(0,10,101);
sinevec  = eval_fd(tfine, sinefd);
phdl=plot(tfine, sinevec, 'b-', tvec*10, sinecoef, 'o');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} f(t)')

%
% Section 4.3 Smoothing using Regression Analysis
%

% Section 4.3.1 Plotting the January Thaw

% Figure 4.3

% The data are in file MontrealTemp.mat

load MontrealTemp

thawdata = MontrealTemp(16:47,:);
thawtime = ((16:47)+0.5);
plot(thawtime, mean(thawdata,2), 'o-')
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Temperature (deg C)')

thawRng = [16,48];
nbasis = 7;
thawbasis    = create_bspline_basis(thawRng, nbasis);
thawbasismat = eval_basis(thawtime, thawbasis);

% Figure 4.4

thawcoef = thawbasismat\thawdata;
thawfdnames{1} = 'Day';
thawfdnames{2} = 'Year';
thawfdnames{3} = 'Temperature (deg C)';
thawfd = fd(thawcoef, thawbasis, thawfdnames);
thawfine = linspace(16,48,101);
thawmat = eval_fd(thawfine, thawfd);
phdl = plot(thawfine, thawmat, 'b-', thawRng, [0,0], 'b--');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Temperature (deg C)')
axis([thawRng,-35,10])

% Figure 4.5

plotfit_fd(thawdata(:,1), thawtime, thawfd(1))

%
% Section 4.4 The Linear Differential Operator or Lfd Class
%

omega = 2*pi/365;
const_basis = create_constant_basis(yearRng);

Lcoef   = [0,omega^2,0];    
wfd     = fd(Lcoef,const_basis);  
wfdcell = fd2cell(wfd);         
harmaccelLfd = Lfd(3, wfdcell);  

accelLfd = int2Lfd(2);

disp(class(accelLfd))
disp(class(harmaccelLfd))

Tempfd = smooth_basis(time, tempav, tempbasis);

D2tempfd = deriv_fd(Tempfd, 2);
Ltempfd  = deriv_fd(Tempfd, harmaccelLfd);

%
% Section 4.5 Bivariate Functional Data Objects:
%             Functions of Two Arguments
%

Bspl2 = create_bspline_basis(unitRng,2,1);
Bspl3 = create_bspline_basis(unitRng,3,2);

corrmat  = reshape((1:6)/6, 2, 3);
bBspl2_3 = bifd(corrmat, Bspl2, Bspl3);

%
% Section4.6 The structure of the fd and Lfd Classes
%

help fd
help Lfd

%
% Section 4.7 Some Things to Try
%
% (exercises for the reader)
