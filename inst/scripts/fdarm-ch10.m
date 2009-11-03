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
% ch.  10.  Linear Models for Functional Responses
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
% Section 10.1  Functional Responses and an Analysis of Variance Model
%

%  ----------------- Climate zone effects for temperature -----------------

%  Section 10.1.1 Climate Region Effects on Temperature

%  path to the daily weather data example folder

weatherPath = [examplesPath,'/weather'];

addpath(weatherPath)

load daily

%  set up the data for the analysis

tempav       = daily.tempav;
precav       = daily.precav;
yearRng      = daily.rng;
time         = daily.time;
place        = daily.place;
dayperiod    = daily.period;
region_names = daily.region_names;
regions      = daily.regions;

% tempfd from chapter 9

%  define the harmonic acceleration operator

harmaccelLfd = vec2Lfd([0,(2*pi/dayperiod)^2,0], dayperiod);

%  organize data to have winter = the center of the plot

dayOfYearShifted = [182:365, 1:181];
tempavShifted    = tempav(dayOfYearShifted,:);

tempbasis65 = create_fourier_basis(yearRng,65);
lambda      = 1e6;
tempfdPar65 = fdPar(tempbasis65, harmaccelLfd, lambda);
tempfd      = smooth_basis(time, tempavShifted, tempfdPar65);
           
station{1} = 'Weather Station';
station{2} = place;
fdnames{1} = 'Day (July 1 to June 30)';
fdnames{2} = station;
fdnames{3} = 'Temperature (deg C)';
tempfd     = putnames(tempfd, fdnames);

%  augment tempfd by adding a 36th observation with temp(t) = 0

coef     = getcoef(tempfd);
coef36   = [coef,zeros(65,1)];
temp36fd = fd(coef36,tempbasis65,fdnames);

%  set up the covariates

p = 5;
regionCell = cell(p,1);
regionCell{1} = ones(36,1);
for j=2:p
    xj = zeros(36,1);
    xj(regions == j-1) = 1;
    xj(36) = 1;
    regionCell{j} = xj;
end

%  set up the regression coefficient list

betabasis = create_fourier_basis(yearRng, 11);
betafdPar = fdPar(betabasis);
betaCell  = cell(p,1);
for j = 1:p 
    betaCell{j} = betafdPar;
end

%  carry out the functional analysis of variance

fRegressStr = fRegress(temp36fd, regionCell, betaCell);

%  extract the estimated regression coefficients and y-values

betaestCell = fRegressStr.betahat;
regionFit   = fRegressStr.yhat;

% Figure 10.1

for j = 1:p 
    subplot(2,3,j)
    plot(getfd(betaestCell{j}))
    xlabel('')
    ylabel(['\fontsize{13} ',region_names(j,:)])
end
subplot(2,3,6)
plot(regionFit)
xlabel('') 
ylabel('\fontsize{13} Prediction')

% ------  10.1.2 Trends in Sea Bird Populations on Kodiak Island  ---------

% 10.1.2 Trends = Sea Bird Populations on Kodiak Island

seabirdPath = [examplesPath,'/seabirds'];

addpath(seabirdPath)

load seabirds

%  select only the data for sites Uyak and Uganik, which have data
%  from 1986 to 2005, except for 1998

sites        = seabirds.Bay;
birdlabels   = seabirds.birdlabels;
birdcountmat = seabirds.birdcountmat;

%  year indices for which there are data.

selYear = [1:12, 14:20];  
yearObs = selYear + 1985;
yearRng = [1986,2005];

%  keep only rows for which there are data

birdcountmat = birdcountmat(selYear,:);

%  Compute mean counts taken over both sites and transects

meanCounts = (birdcountmat(:, 1:13) + ...
              birdcountmat(:,14:26))/2;

%  Compute log (base 10) mean counts

logCounts = log10(meanCounts);

% Figure 10.2

shellfishindex = [1,2,5,6,12,13];
fishindex      = [3,4,7,8,9,10,11];
logcountlim    = [min(min(logCounts)), max(max(logCounts))];

meanShellfish = mean(meanCounts(:, shellfishindex), 2);

subplot(2,1,1)
plot(yearObs, logCounts(:,shellfishindex), 'bo-', ...
     yearObs, log10(meanShellfish),        'b--', ...
     yearRng, [0,0], 'b:')
xlabel('\fontsize{13} year') 
ylabel('\fontsize{13} Log_{10} count')
title('\fontsize{13} Shellfish Diet')
axis([yearRng, logcountlim])

meanfish = mean(meanCounts(:,fishindex), 2);

subplot(2,1,2)
plot(yearObs, logCounts(:,fishindex), 'bo-', ...
     yearObs, log10(meanfish),        'b--', ...
     yearRng, [0,0], 'b:')
xlabel('\fontsize{13} year') 
ylabel('\fontsize{13} Log_{10} count')
title('\fontsize{13} fish Diet')
axis([yearRng, logcountlim])

%  Compute mean counts taken over transects only within sites
%  so we have 2 observations for each bird species each year.
%  Two of these counts are zero, and are replaced by 1/(2*n)

meanCounts2 = birdcountmat;

for j = 1:13
    for k = 1:2 
        jcol = j + (k-1)*13;
        n = sum(meanCounts2(:,jcol));
        zeroind = find(meanCounts2(:,jcol) == 0);
        if ~isempty(zeroind) 
            meanCounts2(zeroind,jcol) = 1/(2*n);
        end
    end
end

logCounts2 = log10(meanCounts2);

%  Represent log mean counts exactly with a polygonal basis

birdbasis = create_polygonal_basis([1,20],selYear);
[birdfd2, df, gcv, SSE, penmat, y2cMap] = ...
               smooth_basis(selYear, logCounts2, birdbasis);

%  -----------------------------------------------------------------
%  After some preliminary analyses we determined that there was no
%  contribution from either site or food*site interaction.
%  Now we use a reduced model with only a feed effect,
%  but we add bird effects, which were seen = the plot to be
%  strong.  Birds are nested within feed groups, and either their
%  effects must sum to zero within each group, or we must designate
%  a bird = each group as a baseline, and provide dummy variables
%  for the remainder.  We opt for the latter strategy.
%  -----------------------------------------------------------------

%  The design matrix contains an intercept dummy variable, a
%  feed dummy variable, and dummy variables for birds, excluding
%  the second bird = each group, which turns out to be the each
%  group's most abundant species, and which is designated as the
%  baseline bird for that group.

Zmat0 = zeros(26,15);

%  Intercept or baseline effect

Intercept = ones(26,1);

%  Crustacean/Mollusc feeding effect:  a contrast between the two groups

foodindex = [1,2,5,6,12,13];
fooddummy = zeros(26,1);
fooddummy([foodindex,foodindex+13]) = 1;

%  Bird effect, one for each species

birddummy = diag(ones(13,1));
birdvarbl = [birddummy;birddummy];

%  fill the columns of the design matrix

Zmat0(:,1)    = Intercept;
Zmat0(:,2)    = fooddummy;
Zmat0(:,3:15) = birdvarbl;

%  Two extra dummy observations are added to the functional data
%  object for log counts, and two additional rows are added to
%  the design matrix to force the bird effects within each diet
%  group to equal 0.

birdcoef3 = [getcoef(birdfd2), zeros(19,2)];
birdfd3   = fd(birdcoef3, birdbasis);

Zmat = [Zmat0; zeros(2,15)];
Zmat(28,fishindex+2) = 1;
Zmat(27,shellfishindex+2)  = 1;

p = 15;
xfdcell = cell(p,1);
betacell = xfdcell;
for j = 1:p
    xfdcell{j} = Zmat(:,j);
end

%  set up the functional parameter object for the regression fns.
%  use cubic b-spline basis for intercept and food coefficients

betabasis1 = create_bspline_basis([1,20],21,4,selYear);
lambda     = 10;
betafdPar1 = fdPar(betabasis1,2,lambda);
betacell{1} = betafdPar1;
betacell{2} = betafdPar1;
betabasis2 = create_constant_basis([1,20]);
betafdPar2 = fdPar(betabasis2);
for j = 3:15 
    betacell{j} = betafdPar2;
end

birdRegressStr = fRegress(birdfd3, xfdcell, betacell);

betaestcell = birdRegressStr{4};

% Figure 10.3 is produced = Section 10.2.2 below
% after estimating the smoothing parameter = Section 10.1.3
%
% Here we plot the regression parameters
% without the confidence intervals.

subplot(2,1,1)
plot(getfd(betaestcell{1}))
subplot(2,1,2)
plot(getfd(betaestcell{2}))

%
% Section 10.1.3 Choosing Smoothing Parameters
%

%  Choose the level of smoothing by minimizing cross-validated
%  error sums of squares.

loglam = -2:0.25:4;
SSE_CV = zeros(length(loglam),1);
betafdPari = betafdPar1;
CV_obs = 1:26;
wt     = ones(28,1);
for i = 1:length(loglam)
    disp(loglam(i))
    betafdPari = putlambda(betafdPari, 10^loglam(i));
    betacelli = betacell;
    for j = 1:2 
        betacelli{j} = betafdPari;
    end
    CVi = fRegress_CV(birdfd3, xfdcell, betacelli, wt, CV_obs);
    SSE_CV(i) = CVi;
end

%  Figure 10.4

phdl = plot(loglam, SSE_CV, 'bo-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log smoothing parameter')
ylabel('\fontsize{13} cross validated sum of squares')

%  Cross-validation is minimized at something like lambda = sqrt(10),
%  although the discontinous nature of the CV function is disquieting.

betafdPar1 = putlambda(betafdPar1, 10^0.5);
for j = 1:2 
    betacell{j} = betafdPar1;
end

%  carry out the functional regression analysis

fRegressStr = fRegress(birdfd3, xfdcell, betacell);

%  plot regression functions

betanames{1} = 'Intercept';
betanames{2} = 'Food Effect';

birdBetaestcell = fRegressStr{4};

for j = 1:2 
    subplot(2,1,j)
    betaestParfdj = birdBetaestcell{j};
    betaestfdj    = getfd(betaestParfdj);
    betaestvecj   = eval_fd(selYear, betaestfdj);
	plot(selYear, betaestvecj, 'b-')
    xlabel('\fontsize{13} Year')
    ylabel(['\fontsize{13} ',betanames{j}])
end

%  plot predicted functions

yhatfdobj = fRegressStr{5};
subplot(1,1,1)
plotfit_fd(logCounts2, selYear, yhatfdobj(1:26))

%
% Section 10.2 Functional Responses with Functional Predictors:
%              The Concurrent Model
%

%  Section 10.2.2 Confidence Intervals for Regression Functions

birdYhatmat = eval_fd(selYear, yhatfdobj(1:26));
rmatb       = logCounts2 - birdYhatmat;
SigmaEb     = cov(rmatb');

birdbetastderrCell = fRegress_stderr(fRegressStr, y2cMap, SigmaEb);

plotbeta(birdBetaestcell(1:2), birdbetastderrCell(1:2))

% Section 10.2.3 Knee Angle Predicted from Hip Angle

%  path to the gait data example folder

gaitPath = [examplesPath,'/gait'];

addpath(gaitPath)

%  load the data from gait.mat

load gait

hip  = gait.hip;
knee = gait.knee;
gaitarray = gait.gaitarray;

gaittime  = 0.5:1:19.5;
gaitrange = [0,20];
gaitfine = linspace(0,20,101)';

harmaccelLfd20 = vec2Lfd([0, (2*pi/20).^2, 0], gaitrange);

nbasis = 21;
gaitbasis = create_fourier_basis(gaitrange, nbasis);

gaitLoglam = -4:0.25:0;
nglam = length(gaitLoglam);

% First select smoothing for the raw data

gaitSmoothStats      = zeros(nglam, 3);
gaitSmoothStats(:,1) = gaitLoglam;

%  loop through smoothing parameters

for ilam = 1:nglam 
    lambdai = 10^gaitLoglam(ilam);
    gaitfdPari = fdPar(gaitbasis, harmaccelLfd20, lambdai); 
    [gaitfdi, dfi, gcvi] = smooth_basis(gaittime, gaitarray, gaitfdPari);
    gaitSmoothStats(ilam, 2) = dfi;
    gaitSmoothStats(ilam, 3) = sum(sum(gcvi));
end

%  display and plot GCV criterion and degrees of freedom

disp(gaitSmoothStats)

subplot(1,1,1)
plot(gaitSmoothStats(:,1), gaitSmoothStats(:,3), 'bo-')
xlabel('\fontsize{13} log_{10} lambda')
ylabel('\fontsize{13} GCV')

%  two panel display of df and GCV

subplot(2,1,1)
plot(gaitSmoothStats(:,1), gaitSmoothStats(:,3), 'bo-')
xlabel('\fontsize{13} log_{10} lambda')
ylabel('\fontsize{13} GCV')
subplot(2,1,2)
plot(gaitSmoothStats(:,1), gaitSmoothStats(:,2), 'bo-')
xlabel('\fontsize{13} log_{10} lambda')
ylabel('\fontsize{13} degrees of freedom')

%    GCV is minimized with lambda = 10^(-1.5).

gaitfdPar = fdPar(gaitbasis, harmaccelLfd20, 10^(-1.5)); 

[gaitfd, df, gcv, SSE, penmat, y2cMap] = ...
                     smooth_basis(gaittime, gaitarray, gaitfdPar);

gaitfdnames{1}   = 'Normalized time'; 
gaitfdnames{2}   = 'Child'; 
gaitfdnames{3}   = cell(2,1);
gaitfdnames{3,1} = 'Angle';
gaitfdnames{3,2} = ['Hip ', 'Knee'];

gaitfd = putnames(gaitfd, gaitfdnames);

hipfd  = gaitfd(:,1);
kneefd = gaitfd(:,2);

% Figure 10.5

%  compute mean knee function and its first two derivatives

kneefdMean    = mean(kneefd);

D1kneevecMean = eval_fd(gaitfine, kneefdMean, 1);
D2kneevecMean = eval_fd(gaitfine, kneefdMean, 2);

%  plot the mean knee function and its first two derivatives in a 
%  3-panel plot

%  plot mean knee function

subplot(3,1,1)
kneeMeanvec = eval_fd(gaitfine, kneefdMean);
plot(gaitfine, kneeMeanvec, 'b-', ...
     [7.5, 7.5], [0, 80], 'b--', [14.7, 14.7], [0, 80], 'b--') 
xlabel('')
title('\fontsize{13} Mean Knee Angle')

%  plot first derivative of mean knee function

subplot(3,1,2)
D1kneeRng = [min(D1kneevecMean), max(D1kneevecMean)];
plot(gaitfine, D1kneevecMean, 'b-', [0,20], [0,0], 'b:', ...
     [7.5, 7.5], D1kneeRng, 'b--', [14.7, 14.7], D1kneeRng, 'b--') 
xlabel('')
title('\fontsize{13} Mean Knee Angle Velocity')

%  plot second derivative of mean knee function

subplot(3,1,3)
D2kneeRng = [min(D2kneevecMean), max(D2kneevecMean)];
plot(gaitfine, D2kneevecMean, 'b-', [0,20], [0,0], 'b:', ...
     [7.5, 7.5], D2kneeRng, 'b--', [14.7, 14.7], D2kneeRng, 'b--') 
xlabel('')
title('\fontsize{13} Mean Knee Angle Acceleration')

% Figure 10.6

%  plot a phase/plane diagram for the mean knee function

%  values of the first and second derivative at observation times

D1vec = eval_fd(gaittime, kneefdMean, 1);
D2vec = eval_fd(gaittime, kneefdMean, 2);

%  labels for these points

plotlab = cell(20,1);
for i=1:20,  plotlab{i} = [' ',num2str(i)];  end

%  plot the phase/plane diagram using function phaseplanePlot

subplot(1,1,1)
[phdl, thdl] = phaseplanePlot(kneefdMean, 1, 2, plotlab, gaittime);
set(phdl, 'LineWidth', 2)
set(thdl, 'fontsize', 12)
xlabel('\fontsize{13} Knee Velocity') 
ylabel('\fontsize{13} Knee Acceleration')
axis([-25,25,-15,15])

%  Set up a  functional linear regression to predict knee angle
%  from hip angle

%  It is essential to remove previous cell arrays with the same names

%  define new covariate cell array with the first covariate being
%  the constant function and the second being hip angle

hipxfdCell{1} = ones(39,1); 
hipxfdCell{2} = hipfd;

%  define regression function cell array kneebetaCell

kneebetafdPar   = fdPar(gaitbasis, harmaccelLfd20);
kneebetaCell    = cell(2,1);
kneebetaCell{1} = kneebetafdPar;
kneebetaCell{2} = kneebetafdPar;

%  the regression analysis

kneefRegressStr = fRegress(kneefd, hipxfdCell, kneebetaCell);

% Figure 10.7

% cell array for estimated regression coefficients

kneebetaestCell = kneefRegressStr.betahat;

%  evaluate regression functions 

kneeIntercept = eval_fd(gaitfine, getfd(kneebetaestCell{1}));
kneeHipCoef   = eval_fd(gaitfine, getfd(kneebetaestCell{2}));

% Plot intercept along with the mean knee function

subplot(2,1,1)
phdl=plot(gaitfine, kneeIntercept, 'b-', gaitfine, kneeMeanvec, 'b--', ...
          [0,20], [0,0], 'b:', ...
          [7.5,7.5], [0,80], 'b--', [14.7,14.7], [0,80], 'b--');
set(phdl, 'LineWidth', 2)
xlabel('') 
ylabel('')
title('\fontsize{16} Intercept and Mean Knee Angle')

% Plot hip coefficient

ylim = [0, max(kneeHipCoef)];

subplot(2,1,2)
phdl=plot(gaitfine, kneeHipCoef, 'b-', [0,20], [0,0], 'b:', ...
          [7.5,7.5], ylim, 'b--', [14.7,14.7], ylim, 'b--');
set(phdl, 'LineWidth', 2)
xlabel('') 
ylabel('')
title('\fontsize{16} Hip Coefficient')

% Squared multiple correlation

%  fitted knee functions

kneehatfd  = kneefRegressStr.yhat;

%  matrix of residual values computed at observation times

kneehatmat = eval_fd(gaitfine, kneehatfd);
kneemat    = eval_fd(gaitfine, kneefd);
resmat1    = kneemat - kneehatmat;

%  matrix of residuals for fitting the data with the mean knee angle

kneemeanvec = eval_fd(gaitfine, mean(kneefd));
resmat0     = kneemat - kneemeanvec * ones(1,39);

SSE0  = sum(resmat0.^2, 2);
SSE1  = sum(resmat1.^2, 2);
RSQRD = (SSE0-SSE1)./SSE0;

% Plot Hip Coefficient & Squared Multiple Correlation

subplot(1,1,1)
phdl=plot(gaitfine, kneeHipCoef, 'b-', gaitfine, RSQRD, 'b--', ...
          [7.5,7.5], ylim, 'b--', [14.7,14.7], ylim, 'b--');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} ') 
ylabel('\fontsize{13} ') 
title('\fontsize{16} Hip Coefficient and Squared Multiple Correlation')
legend('\fontsize{12} Hip Coef.', '\fontsize{12} R^2')

% Figure 10.8

%  covariance matrix for residuals

resmat = eval_fd(gaittime, kneefd) - eval_fd(gaittime, kneehatfd);
SigmaE = cov(resmat');

%  compute cell array of standard error functions

kneebetastderrCell = fRegress_stderr(kneefRegressStr, y2cMap, SigmaE);

%  plot regression functions with 95% confidence intervals

plotbeta(kneebetaestCell, kneebetastderrCell)

% Figure 10.9

%  predict knee acceleration from hip acceleration

clear hipxfdcell  %  important to clear earlier version of cell array

hipxfdCell2    = cell(2,1);
hipxfdCell2{1} = ones(39,1); 
hipxfdCell2{2} = deriv(hipfd, 2);

knee_accelfd = deriv(kneefd, 2);

kneefRegressStr2 = fRegress(knee_accelfd, hipxfdCell2, kneebetaCell);

% cell array for estimated regression coefficients

kneebetaestCell2 = kneefRegressStr2.betahat;

%  evaluate regression functions 

kneeIntercept2 = eval_fd(gaitfine, getfd(kneebetaestCell2{1}));
kneeHipCoef2   = eval_fd(gaitfine, getfd(kneebetaestCell2{2}));

% Squared multiple correlation

%  fitted knee functions

knee_accelhatfd  = kneefRegressStr2.yhat;

%  matrix of residual values computed at observation times

knee_accelmat = eval_fd(gaitfine, knee_accelfd);

resmat1 = knee_accelmat - eval_fd(gaitfine, knee_accelhatfd);

%  matrix of residuals for fitting the data with the mean knee angle

knee_accelmeanvec = eval_fd(gaitfine, mean(knee_accelfd), 2);
resmat0  = knee_accelmat - knee_accelmeanvec * ones(1,39);

SSE0  = sum(resmat0.^2, 2);
SSE1  = sum(resmat1.^2, 2);
RSQRD = (SSE0-SSE1)./SSE0;


ylim=[0, max(kneeHipCoef2)];
phdl = plot(gaitfine, kneeHipCoef2, 'b-', ...
            [7.5,7.5], ylim, 'b--', [14.7,14.7], ylim, 'b--');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} ') 
ylabel('\fontsize{13} ') 
title('\fontsize{13} Hip acceleration')

% Squared multiple correlation

phdl = plot(gaitfine, RSQRD, 'b-', ...
            [7.5,7.5], [0,1], 'b--', [14.7,14.7], [0,1], 'b--');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} ') 
ylabel('\fontsize{13} ') 
title('\fontsize{13} R^2')

%
% Section 10.3 Beyond the Concurrent Model
%
%  (no computations = this section)

%
% Section 10.4 A Functional Linear Model for Swedish Mortality
%

% SwedeLogHazardPath = [examplesPath,'/SwedeLogHazard'];
% 
% addpath(SwedeLogHazardPath)
% 
% load SwedeLogHazard
% load SwedeLogHazard1914

%  The Swedish

%  Figure 10.10

%  indices of years 1751, 1810 and 1860

threeyears = [1751, 1810, 1860] - 1750;

tempLogHazard = [SwedeLogHazard(:,threeyears), SwedeLogHazard1914];

SwedeTime = (0:80)';

phdl = plot(SwedeTime, tempLogHazard, '-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} age')
ylabel('\fontsize{13} log hazard')
legend('\fontsize{12} 1751', '\fontsize{12} 1810', ...
       '\fontsize{12} 1860', '\fontsize{12} 1914', ...
       'location', 'north')

%  plot the data for 1751 to 1884

plot(SwedeLogHazard)

%  set up a saturated order six basis for smoothing the data

SwedeRng = [0,80];

nbasis = 85;
norder =  6;
SwedeBasis = create_bspline_basis(SwedeRng,nbasis,norder);

%  set the smoothing parameters

lambda = 1e-7;
SwedefdPar = fdPar(SwedeBasis, 4, lambda);

%  smooth the data

Swedefd  = smooth_basis(SwedeTime,SwedeLogHazard,SwedefdPar);

plotfit_fd(SwedeLogHazard, SwedeTime, Swedefd)

% set up coefficient functions

nbasis = 23;
SwedeBetaBasis = create_bspline_basis(SwedeRng, nbasis);

SwedeBeta0fdPar = fdPar(SwedeBetaBasis, 2, 1e-5);

SwedeBeta1fd    = bifd(zeros(23), SwedeBetaBasis, SwedeBetaBasis);
SwedeBeta1fdPar = bifdPar(SwedeBeta1fd, 2, 2, 1e3, 1e3);

p = 2;
SwedeBetaCell = cell(p, 1);
SwedeBetaCell{1} = SwedeBeta0fdPar;
SwedeBetaCell{2} = SwedeBeta1fdPar;

%  Define the dependent and independent variables

NextYear = Swedefd(2:144);
LastYear = Swedefd(1:143);

Swede_linmodStr = linmod(NextYear, LastYear, SwedeBetaCell);

Swede_beta1fd = eval_bifd(SwedeTime, SwedeTime, ...
                          Swede_linmodStr.beta);

% Figure 10.11

surfc(SwedeTime, SwedeTime, Swede_beta1fd)
xlabel('\fontsize{13} age')
ylabel('\fontsize{13} age')
zlabel('\fontsize{13} \beta(s,t)')

%
% Section 10.5 Permutation Tests of Functional Hypotheses
%
%  Section 10.5.1 Functional t-Tests

%  ---------------  Comparing male and female growth data  ----------------


%  path to the Berkeley growth data example folder

growthPath = [examplesPath,'/growth'];

addpath(growthPath)

%  load the growth data struct object from the .mat file growth.mat

load growth

% Figure 10.12

phdl = plot(age, hgtmmat(:, 1:10), 'b--', ...
            age, hgtfmat(:, 1:10), 'b-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} age')
ylabel('\fontsize{13} height (cm)')
axis([1,18, 60,200])

nbasis = 35;
norder =  6;
growthbasis = create_bspline_basis([1,18], nbasis, norder, age);
growfdPar = fdPar(growthbasis, 3, 10^(-0.5));

hgtffd = smooth_basis(age,hgtfmat,growfdPar);
hgtmfd = smooth_basis(age,hgtmmat,growfdPar);

tresStr = tperm_fd(hgtffd, hgtmfd);

% Figure 10.13

% Section 10.5.2 Functional F-Tests

%  ---------------  Testing for no effect of climate zone  ----------------

% temp36fd, regionCell, betaCell from Section 10.1.1 above

F_res = Fperm_fd(temp36fd, regionCell, betaCell);

% Figure 10.14

qval = 0.95;
phdl = plot(argvals, F_res.Fvals, 'b-', ...
            argvals, F_res.qvals_pts, 'b--', ...
            [min(argvals), max(argvals)], [qval,qval], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} day')
ylabel('\fontsize{13} F-statistic')
legend('\fontsize{13} Observed Statistic', ...
       ['\fontsize{13} pointwise ', 1-q, ' critical value'], ...
       ['\fontsize{13} maximum ',   1-q, ' critical value'])

%
% 10.6 Details for R Functions fRegress, fRegress_CV and fRegress_stderr
%

help fRegress
help fRegress_CV
help fRegress_stderr

%
% 10.7 Details for Function plotbeta
%

help plotbeta

%
% 10.8 Details for Function linmod
%
help linmod

%
% 10.9 Details for Functions Fperm_fd and tperm_fd
%
help(Fperm_fd)
help(tperm_fd)

%
% Section 10.10 Some Things to Try
%
% (exercises for the reader)

%
% Section 10.11  More to Read
%
