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
% ch. 9.  Functional Linear Models for Scalar Response
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
% Section 9.1 Functional Linear regression with a Scalar response
%

%  (no computations in this section)

%
% Section 9.2 A Scalar Response Model for Log Annual Precipitation
%

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

harmaccelLfd = vec2Lfd([0,(2*pi/dayperiod)^2,0],yearRng);    

%  organize data to have winter in the center of the plot

dayOfYearShifted = [182:365, 1:181];

%  change 0's to 0.05 mm in precipitation data

prectmp = precav;
for j=1:35
    index = find(prectmp(:,j)==0);
    prectmp(index,j) = 0.05;
end

%  set up functional data object for log precipitation

logprecav = log10(prectmp);

%  set up a saturated basis: as many basis functions as observations

daybasis  = create_fourier_basis(yearRng, 365);

%  smooth data with lambda that minimizes GCV

lambda     = 1e6;
fdParobj   = fdPar(daybasis, harmaccelLfd, lambda);
logprec_fd = smooth_basis(time, logprecav, fdParobj);
station{1} = 'Weather Station';
station{2} = place;
fdnames{1} = 'Day (July 1 to June 30)';
fdnames{2} = station;
fdnames{3} = 'Log 10 Precipitation (mm)';
logprec_fd  = putnames(logprec_fd, fdnames);

annualprec = log10(sum(precav,1))';

tempbasis65 = create_fourier_basis(yearRng,65);
[tempfd65, df, gcv, SSE, penmat, y2cMap] = ...
                       smooth_basis(time, tempav, tempbasis65);

%
% Section 9.3 Setting Up the Functional Linear Model
%
%  (no computations in this section)

%
% Section 9.4 Three Estimates of the Regression Coefficient
%             Predicting Annual Precipitation
%

tempcell{1} = fd(ones(1,35),create_constant_basis(yearRng));
tempcell{2} = tempfd65;

% 9.4.1 Low Dimensional Regression Coefficient Function beta

conbasis     = create_constant_basis(yearRng);
betabasis5   = create_fourier_basis(yearRng,5);
betacell1{1} = conbasis;
betacell1{2} = betabasis5;

fRegressStr1 = fRegress(annualprec,tempcell,betacell1);

betaestcell1  = fRegressStr1.betahat;
tempbetafd1   = getfd(betaestcell1{2});

% Figure 9.1

plot(tempbetafd1) 
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Beta for temperature')

%  display intercept

intercept = getcoef(getfd(betaestcell1{1}));
disp(['Constant intercept = ',num2str(intercept)])

% 0.0095 as in the book

annualprechat1 = fRegressStr1.yhat;
annualprecres1 = annualprec - annualprechat1;
SSE1_1 = sum(annualprecres1.^2);
SSE0   = sum((annualprec - mean(annualprec)).^2);

RSQ1   = (SSE0-SSE1_1)/SSE0;
disp(['Squared multiple correlation = ',num2str(RSQ1)])

Fratio1 = ((SSE0-SSE1_1)/5)/(SSE1_1/29);
disp(['F-ratio = ',num2str(Fratio1)])

% 22.6 as in the book

% 9.4.2 Coefficient beta Estimate Using a Roughness Penalty

% refit with 35 terms rather than 5 in the fourier basis

betabasis35 = create_fourier_basis(yearRng, 35);
lambda      = 10^12.5;
betafdPar   = fdPar(betabasis35, harmaccelLfd, lambda);

betacell2    = betacell1;
betacell2{2} = betafdPar;

annPrecTempStr = fRegress(annualprec, tempcell, betacell2);
betaestcell2   = annPrecTempStr.betahat;
annualprechat2 = annPrecTempStr.yhat;

disp(['Degrees of freedom for fit = ',num2str(annPrecTempStr.df)])

SSE1_2 = sum((annualprec-annualprechat2).^2);

RSQ2 = (SSE0 - SSE1_2)/SSE0;
disp(['Squared multiple correlation = ',num2str(RSQ2)])

% 0.75 as in the book

Fratio2 = ((SSE0-SSE1_2)/3.7)/(SSE1_2/30.3);
disp(['F-ratio = ',num2str(Fratio2)])

% 25.1 as in the book

% Figure 9.2

plot(annualprechat2, annualprec,     'o',  ...
     annualprechat2, annualprechat2, 'b-')
xlabel('\fontsize{13} annualprechat') 
ylabel('\fontsize{13} annualprec')

% Figure 9.3
% ... see section 9.4.4 below ...

plot(getfd(betaestcell2{2}))

% Compare with the constant fit:

betacell     = betacell1;
betacell{2}  = fdPar(conbasis);
fRegressStr = fRegress(annualprec, tempcell, betacell);

betaestcell   = fRegressStr.betahat;
annualprechat = fRegressStr.yhat;

SSE1 = sum((annualprec-annualprechat).^2);

RSQ = (SSE0 - SSE1)/SSE0;
disp(['Squared multiple correlation = ',num2str(RSQ)])

% 0.49 as in the book

Fratio = ((SSE0-SSE1)/1)/(SSE1/33);
disp(['F-ratio = ',num2str(Fratio)])

% 31.3 as in the book

% 9.4.3 Choosing Smoothing Parameters

loglam = 5:0.5:15;
nlam   = length(loglam);
SSE_CV = zeros(nlam,1);
for ilam = 1:nlam
  disp(['log lambda = ', num2str(loglam(ilam))])
  lambda     = 10^(loglam(ilam));
  betacelli  = betacell2;
  betafdPar2 = betacelli{2};
  betafdPar2 = putlambda(betafdPar2,lambda);
  betacelli{2} = betafdPar2;
  SSE_CVi      = fRegress_CV(annualprec, tempcell, betacelli);
  SSE_CV(ilam) = SSE_CVi;
end

%  Figure 9.4

phdl=plot(loglam, SSE_CV, 'bo-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log_{10} smoothing parameter lambda')
ylabel('\fontsize{13} Cross-validation score')

% 9.4.4 Confidence Intervals

N       = length(annualprec);
resid   = annualprec - annualprechat2;
SigmaE  = sum(resid.^2)/(35-annPrecTemp{10});
SigmaE  = SigmaE.*diag(ones(N,1));

%  obtain the standard error functions using function
%  fRegress_stderr

betastderrCell = fRegress_stderr(annPrecTemp, y2cMap, SigmaE);
betastderrfd   = betastderrCell{2};

betafdPar      = betaestcell2{2};
betafd         = getfd(betafdPar);

% Figure 9.3

betavec       = eval_fd(time, betafd);
betastderrvec = eval_fd(time, betastderrfd);

phdl=plot(time, betavec, 'b-', ...
          time, betavec+2*betastderrvec, 'b--', ...
          time, betavec-2*betastderrvec, 'b--', ...
          yearRng, [0,0], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Temperature Reg. Coeff.')
axis([yearRng, -6e-4, 1.2e-03])

%  Optionally, function fRegress can also compute the
%  standard error functions, provided the last two arguments
%  are supplied

temp = fRegress(annualprec, tempcell, betacell2, ones(N,1), ...
                y2cMap, SigmaE);
betastderrCell = temp{13};
betastderrfd   = betastderrCell{2};

% Section 9.4.5 Scalar Response Models by Functional Principal Components

daybasis365   = create_fourier_basis(yearRng, 365);
lambda        = 1e6;
tempfdPar365  = fdPar(daybasis365, harmaccelLfd, lambda);
tempfd365     = smooth_basis(time, tempav, tempfdPar365);

lambda     = 1e0;
tempfdPar  = fdPar(daybasis365, harmaccelLfd, lambda);
temppcastr = pca_fd(tempfd365, 4, tempfdPar);
harmonics  = temppcastr.harmfd;
scores     = temppcastr.harmscr;
Xmat       = [ones(35,1),scores];
pcacoefs   = Xmat\annualprec;
annprechat = Xmat*pcacoefs;

harmmat = eval_fd(time, harmonics);

annprecres = annualprec - annprechat;

annprecvar = sum(annprecres.^2)./(35-5);

betavec    = pcacoefs(2).*harmmat(:,1) + ...
             pcacoefs(3).*harmmat(:,2) + ...
             pcacoefs(4).*harmmat(:,3);
coefvar    = annprecvar.*diag(inv(Xmat'*Xmat));
betavarvec = coefvar(2).*harmmat(:,1).^2 + ...
             coefvar(3).*harmmat(:,2).^2 + ...
             coefvar(4).*harmmat(:,3).^2;

% Figure 9.5

phdl=plot(time, betavec, 'b-',  ...
          time, betavec+2*sqrt(betavarvec), 'b--', ...
          time, betavec-2*sqrt(betavarvec), 'b--', ...
          yearRng, [0,0], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Regression Coef.')
axis([yearRng,-6e-4,1.2e-03])

%
% Section 9.5 Statistical Tests
%

F_resstr = Fperm_fd(annualprec, tempcell, betacell);

disp(F_resstr.Fobs)

disp(F_resstr.qval)

%
% Section 9.6 Some Things to Try
%
% (exercises for the reader)

%
% Section 9.7  More to Read
%
