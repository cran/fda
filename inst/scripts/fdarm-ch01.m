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

%  We want to call you attention particularly to a new package in the
%  R project at www.R-project.org called 'fds' containing
%  functional data that has appeared since the publication of 
%  book.  It contains most of the data that we analyze here, as well as
%  many new data sets.

%
% Chapter 1.  Introduction
%

%  Set up some strings for constructing paths to folders.
%  These strings should be modified so as to provided access
%  to the specified folders on your computer.

%  Path to the folder containing the Matlab functional data analysis
%  software

fdaMPath = 'c:/Program Files/MATLAB/R2009a/fdaM';

addpath(fdaMPath)

%  Path to the folder containing the examples

%  Here is the full path name:

examplesPath = 'c:/Program Files/MATLAB/R2009a/fdaM/examples';

%  Here is a path name constructed from fdaMPath:

examplesPath = [fdaMPath,'/examples'];

addpath(examplesPath)

%
% Section 1.1.  What are functional data?
%

% Plot the Berkeley growth data

%  path to the Berkeley growth data example folder

growthPath = [examplesPath,'/growth'];

addpath(growthPath)

% --------------  Plot the Berkeley growth data  --------------------------

%  load the growth data struct object from the .mat file growth.mat

load growth

%  the dimensions of the data

ncasem = growth.ncasem;
ncasef = growth.ncasef;
nage   = growth.nage;

%  the heights of the 39 boys

hgtmmat = growth.hgtmmat;

%  the heights of the 54 girls

hgtfmat = growth.hgtfmat;

%  the 31 ages of measurement

age = growth.age;

%  define the range of the ages and set up a fine mesh of ages

ageRng  = [1,18];
agefine = linspace(ageRng(1), ageRng(2), 501);

% Monotone smooth (see section 5.4.2)
% B-spline of order 6 = quintic polynomials
% so the acceleration will be cubic

gr_norder = 6;
gr_breaks = age;
gr_nbasis = length(gr_breaks) + gr_norder - 2;
gr_basis  = create_bspline_basis(ageRng, gr_nbasis, gr_norder, gr_breaks);

% Consider only the first 10 girls

children = 1:10;

% matrix of starting values for coefficients

cvecf    = zeros(gr_nbasis, ncasef);

% Create an initial functional data object

gr_fd0   = fd(cvecf, gr_basis);

% Create an initial functional parameter object
% Lfdobj = 3 to penalize the rate of change of acceleration

gr_Lfd    = 3;

% Figure 1.15 was created with lambda = 10^(-1.5);
% we use that also for Figure 1.1

gr_lambda = 10^(-1.5);

%  Define the functional parameter object

gr_fdPar  = fdPar(gr_fd0, gr_Lfd, gr_lambda);

%  Monotonically smooth the female data  

[gr_Wfd, gr_beta, gr_yhatfd] = ...
          smooth_monotone(age, hgtfmat(:,children), gr_fdPar);

%  Evaluate height at a fine mesh of values

hgtf_vec = eval_fd(agefine, gr_yhatfd);

%  Plot Figure 1.1

plot(age, hgtfmat(:,children), 'bo', agefine, hgtf_vec, 'b-')
xlabel('\fontsize{13}Age (years)') 
ylabel('\fontsize{13}Height (cm)')
axis([1,18,60,190])

% Figure 1.2

accfvec = eval_fd(agefine, gr_yhatfd, 2);
accfvecmean = mean(accfvec,2);
phdl = plot(agefine, accfvec, 'b-', [1,18], [0,0], 'b:');
set(phdl, 'LineWidth', 1)
lhdl = line(agefine, accfvecmean);
set(lhdl, 'LineWidth', 2, 'LineStyle', '--', 'color', 'b')
xlabel('\fontsize{13}Age (years)') 
ylabel('\fontsize{13}Acceleration (cm/yr^2)')
axis([1,18,-4,2])

%  ----------------  Plot nondurable manufacturing goods index  -----------

%  path to the nondurable goods index example folder

goodsindexPath = [examplesPath,'/goodsindex'];

addpath(goodsindexPath)

%  load the goodsindex data from the .mat file goodsindex.mat

load goodsindex

ndur        = goodsindex.ndur;
yeartime    = goodsindex.yeartime;
nondurables = goodsindex.nondurables;
lognondur   = goodsindex.lognondur;

% Figure 1.3

phdl = plot(yeartime, nondurables);
set(phdl, 'LineWidth', 1)
xlabel('\fontsize{13} Year')
ylabel('\fontsize{13} Nondurable Goods Index')
axis([1919,2000,0,120])

%  compute linear trend for log index

xmat = ones(ndur,2);
xmat(:,2) = yeartime';
lognondurhat = xmat*(xmat\lognondur);

% Figure 1.4

phdl = plot(yeartime, log10(nondurables), 'b-', ...
            yeartime, lognondurhat, 'b--');
set(phdl, 'LineWidth', 1)
xlabel('\fontsize{13} Year')
ylabel('\fontsize{13} Log10 Nondurable Goods Index)')
axis([1919,2000,0.7,2.2])

%  ------------------------  Plot refinery data  --------------------------

%  path to the refinery data example folder

refineryPath = [examplesPath,'/refinery'];

addpath(refineryPath)

%  input the data from file refinery.mat

load refinery

%  Input and output values have baselines 0 here

Time   = refinery.Time;    %  observation time
Reflux = refinery.Reflux -  20.4013;  %  reflux flow
Tray47 = refinery.Tray47 - 215.2583;  %  tray 47 level

% Figure 1.5

subplot(2,1,1)
plot(Time, Tray47, '.')
xlabel('') 
ylabel('\fontsize{13} Tray 47 level')
subplot(2,1,2)
plot(Time, Reflux, '.') 
xlabel('\fontsize{13} Time (min)')
ylabel('\fontsize{13} Reflux flow')

%
% Section 1.2.  Multivariate functional data
%

%  --------------------------  Plot gait data  ----------------------------

%  path to the gait data example folder

gaitPath = [examplesPath,'/gait'];

addpath(gaitPath)

%  load the data from gait.mat

load gait

hip  = gait.hip;
knee = gait.knee;
gaitarray = gait.gaitarray;

Time = [0,linspace(0.025, 0.975, 20),1]';

Gait = zeros(22, 39, 2);
Gait(2:21,:,1) = hip;
Gait(2:21,:,2) = knee;

% Interpolate a common value for 0 & 1

gait01       = 0.5.*(gaitarray(20,:,:) + gaitarray(1,:,:));
Gait( 1,:,:) = gait01;
Gait(22,:,:) = gait01;

before = 17:22;
Before = [-0.225:0.05:-0.025, 0];

after = 1:6;
After = [1, 1.025:0.05:1.225];

% Figure 1.6

subplot(2,1,1)
plot(Time, Gait(:,:,1), 'b-', ...
     Before, Gait(before,:,1), 'b:', After, Gait(after,:,1), 'b:') 
xlabel('')
ylabel('\fontsize{13}Hip Angle (degrees)')
axis([-0.25, 1.25,-10,80]) 
subplot(2,1,2)
plot(Time, Gait(:,:,2), 'b-', ...
     Before, Gait(before,:,2), 'b:', After, Gait(after,:,2), 'b:') 
xlabel('\fontsize{13} Time (portion of gait cycle)')
ylabel('\fontsize{13} Knee Angle (degrees)')
axis([-0.25, 1.25,-10,80]) 

GaitMean = mean(Gait, 2);

xtemp = [Gait(:,:,1), GaitMean(:, 1)];
ytemp = [Gait(:,:,2), GaitMean(:, 2)];

% Figure 1.7

subplot(1,1,1)
plot(Gait(:,1,1),   Gait(:,1,2),   'bo-', ...
     GaitMean(:,1), GaitMean(:,2), 'bo--')
xlabel('\fontsize{13} Hip angle (degrees)') 
ylabel('\fontsize{13} Knee angle (degrees)')
axis([min(min(xtemp)),max(max(xtemp)),min(min(ytemp)),max(max(ytemp))])
i4 = 5:4:21;
LETTERS = ['B';'C';'D';'E';'A'];
text(Gait(i4,1,1), Gait(i4,1,2),   LETTERS)
text(GaitMean(i4,1),    GaitMean(i4,2),      LETTERS)

%  ---------------------  Plot "fda" handwriting data  --------------------

%  path to the handwriting of "fda" data example folder

handwritPath = [examplesPath,'/handwrit'];

addpath(handwritPath)

%  load the data from file fda.mat

load fda

handwrit = fda.fdaarray;
fdatime  = fda.fdatime;
fdarange = fda.fdarange;

% Figure 1.8

plot(100*handwrit(:,:,1), 100*handwrit(:,:,2), '-')
xlabel('') 
ylabel('')

% Figure 1.9

%  -------------  Plot "statistics" in Chinese handwriting data  ----------

%  path to the handwriting of "statistical science" in Chinese
%  example folder

ChinaScriptPath = [examplesPath,'/ChinaScript'];

addpath(ChinaScriptPath)

%  load the data from file ChinaScript.mat

load ChinaScript

StatSciChinese = ChinaScript.penpos;

mark = 1:12:601;

i  = 1;
StatSci1 = StatSciChinese(:,i,:);
% Where does the pen leave the paper?
thresh = quantile(StatSci1(:, 3), .8);

sel1 = find(StatSci1(:,3) >= thresh);
StatSci1(sel1,1:2) = NaN;
phdl=plot(StatSci1(:,1), StatSci1(:,2), 'b-', ...
          StatSci1(mark, 1), StatSci1(mark, 2), 'o');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} X') 
ylabel('\fontsize{13} Y')

%
% Section 1.3.  Functional models for nonfunctional data
%

% Figure 1.10:

% Based on a functional logistic regression.
% Unfortunately, the current 'fda' package does NOT include
% code to estimate a functional logistic regression model.

%
% Section 1.4.  Some functional data analyses
%

%  -----------------------  Plot Canadian weather data  -------------------

%  Plot temperature curves for four Canadian weather stations

%  path to the daily weather data example folder

weatherPath = [examplesPath,'/weather'];

addpath(weatherPath)

load daily

tempav = daily.tempav;
precav = daily.precav;
rng    = daily.rng;
time   = daily.time;
place  = daily.place;

fig1_11Stns = [12, 23, 29, 35];
fig1_11Temp = tempav(:, fig1_11Stns);

% Comment from the end of section 1.5.1:
% "the temperature data in Figure 1.11 were fit using
%  smoothing splines".

day5 = time(1:5:365);
Temp_fourier    = create_fourier_basis(rng, 13);
fig1_11Temp_fd  = smooth_basis(time, fig1_11Temp, Temp_fourier);
fig1_11Temp_mat = eval_fd(day5, fig1_11Temp_fd); 

% Figure 1.11

plot(day5, fig1_11Temp_mat)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13}Mean Temperature (deg C)')
axis([0,365,-40,30])
legend(place(12,:),place(23,:),place(29,:),place(35,:),'Location','South')

%  plot the harmonic accelerations of temperature curves

%  set up the harmonic acceleration linear differential operator

Lbasis  = create_constant_basis(rng);  %  create a constant basis
Lcoef   = [0,(2*pi/365)^2,0];    %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

%  evaluate the harmonic acceleration curves

Lfd_Temp_mat = eval_fd(day5, fig1_11Temp_fd, harmaccelLfd);

% Figure 1.12

plot(day5, Lfd_Temp_mat, '-', [0,365], [0,0], 'b--')
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} L-Temperature')
axis([0,365,[-10,8]*1e-4])
legend(place(12,:),place(23,:),place(29,:),place(35,:),'Location','South')

%
% Section 1.5.  The first steps in a functional data analysis
%

%  Plot the precipitation for Prince Rupert

P_RupertPrecip = precav(:,29);

TempBasis = create_fourier_basis(rng, 13);
TempfdPar = fdPar(TempBasis, harmaccelLfd, 10^7);
P_Rupert_Prec_fd  = smooth_basis(time, P_RupertPrecip, TempfdPar);
P_Rupert_Prec_mat = eval_fd(day5, P_Rupert_Prec_fd);


% Figure 1.13

plot(time, P_RupertPrecip, 'b.') 
lhdl = line(day5, P_Rupert_Prec_mat);
set(lhdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Precipitation (mm)')
axis([0,365,0,18])

%  -----------------------  Plot pinch force data  ------------------------

%  path to the pinch force data example folder

pinchPath = [examplesPath,'/pinch'];

addpath(pinchPath)

%  load the data from file pinch.mat

load pinch

time  = pinch.time;
force = pinch.force;

% Figure 1.14

plot(time, force, 'b-')
xlabel('\fontsize{13}seconds')
ylabel('\fontsize{13}Force (N)')

%  ----------  Plot phase-plane diagrams for the growth data  -------------

midspurt_index = find(abs(agefine-11.7) == min(abs(agefine-11.7)));

% Use the fit from Figure 1.1

hgtf_vel = eval_fd(gr_yhatfd, agefine, 1);
hgtf_acc = eval_fd(gr_yhatfd, agefine, 2);

% Figure 1.15

phdl=plot(hgtf_vel(:,1), hgtf_acc(:,1), 'b-', ...
          [0,12], [0,0], 'b:', ...
          hgtf_vel(midspurt_index,1), ...
          hgtf_acc(midspurt_index,1), 'b.');
set(phdl, 'LineWidth', 2)
hold on
for i = 2:10
     phdl=plot(hgtf_vel(:,i), hgtf_acc(:,i), 'b-', ...
               hgtf_vel(midspurt_index,i), ...
               hgtf_acc(midspurt_index,i), 'bo');
     set(phdl, 'LineWidth', 2)
end
hold off
xlabel('\fontsize{13} Velocity (cm/yr)')
ylabel('\fontsize{13} Acceleration (cm/yr^2)')
axis([0,12,-5,2])

%
% Section 1.6.  Exploring variability in functional data
%
% (no data analysis in this section)

%
% Section 1.7.  Functional linear models
%
% (no data analysis in this section)

%
% Section 1.8.  Using derivatives in functional data analysis
%
% (no data analysis in this section)

%
% Section 1.9.  Concluding remarks
%
% (no data analysis in this section)

%
% Section 1.10.  Some Things to Try
%
