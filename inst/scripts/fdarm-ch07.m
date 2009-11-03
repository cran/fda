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
% ch. Chapter 7  Exploring Variation: Functional Principal
%                and Canonical Components analysis
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
% Section 7.1 An Overview of Functional PCA
%
%  (no computations in this section)

%
% Section 7.2 PCA with Function pca.fd
%

% Section 7.2.1 PCA of the Log Precipitation Data

%  set up a functional data object for log precipitation

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

Lcoef = [0,(2*pi/dayperiod)^2,0];    
harmaccelLfd = vec2Lfd(Lcoef, yearRng); 

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

%  smooth data with lambda that minimizes GCV

lambda      = 1e6;
fdParobj    = fdPar(daybasis, harmaccelLfd, lambda);
logprec_fd = smooth_basis(time, logprecav, fdParobj);
station{1} = 'Weather Station';
station{2} = place;
fdnames{1} = 'Day (July 1 to June 30)';
fdnames{2} = station;
fdnames{3} = 'Log 10 Precipitation (mm)';
logprec_fd  = putnames(logprec_fd, fdnames);

%  do principal component analysis with 2 components

nharm = 2;
logprec_pcastr = pca_fd(logprec_fd, nharm);

disp(logprec_pcastr.values(1:4))

% Figure 7.1

plot_pca_fd(logprec_pcastr)

%  The expansion supplied by the function is too large,
%  and here we supply a smaller value, 0.5

plot_pca(logprec_pcastr, 1, 0, 0.5)

% Figure 7.2

logprec_rotpcastr = varmx_pca_fd(logprec_pcastr);
plot_pca_fd(logprec_rotpcastr, 1, 0, 0.5)

% Figure 7.3

pcascores = logprec_rotpcastr.harmscr;

subplot(1,1,1)
phdl=plot(pcascores(:,1), pcascores(:,2), 'bo');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Rotated Harmonic I ')
ylabel('\fontsize{13} Rotated Harmonic II')

%  labels for selected score points must be plotted interactively

% Section 7.2.2 PCA of Log Precipitation Residuals
% logprecres = residuals from
% the smooths of the log precipitation curves in Chapter 5.

logprecmat = eval_fd(time, logprec_fd);
logprecres = logprecav - logprecmat;

% Figure 7.4

logprecres_fd = smooth_basis(time, logprecres, fdParobj);

plot(logprecres_fd)
xlabel('\fontsize{13} Day') 
ylabel('\fontsize{13} Residual (log 10 mm)')
axis([0,365,-0.07, 0.07])

% Figure 7.5

logprec_str1 = pca_fd(logprecres_fd, 1);
plot_pca_fd(logprec_str1, 1, 0, 0.01)

%
% Section 7.3 More Functional PCA Features
%

%  (no computations in this section)

%
% Section 7.4 PCA of joint X-Y Variation in Handwriting
%

%  Define time values and order 6 spline basis with 105 basis functions
%  This places a knot at every 23rd observation point, and is found to
%  correspond closely to spline smoothing results.

%  path to the handwriting of 'fda' data example folder

handwritPath = [examplesPath,'/handwrit'];

addpath(handwritPath)

%  load the data from file fda.mat

load fda

handwrit = fda.fdaarray;
fdatime  = fda.fdatime;
fdarange = fda.fdarange;

nbasis = 105;
norder = 6;
fdabasis = create_bspline_basis(fdarange, nbasis, norder);

%  set up the functional data structure

fda_fd = smooth_basis(fdatime, handwrit, fdabasis);

clear fdnames
fdnames{1} = 'Milliseconds';
fdnames{2} = 'Replications';
varnames = cell(2,1);
varnames{1} = 'coordinates';
varnames{2} = ['X';'Y'];
fdnames{3} = varnames;
fda_fd = putnames(fda_fd, fdnames);

%  plot the data

plot(fda_fd)

%  a principal components analysis

nharm = 3;
fdapcastr = pca_fd(fda_fd, nharm);

plot_pca_fd(fdapcastr, 1, 0, 0.2);

fdarotpcastr = varmx_pca_fd(fdapcastr);

plot_pca_fd(fdarotpcastr, 1, 0, 0.2)

fdaeig = fdapcastr.values;
neig = 12;
x = ones(neig-nharm,2);
x(:,2) = (nharm+1):neig;
y = log10(fdaeig((nharm+1):neig));
c = x\y;

% Figure 7.6

subplot(1,1,1)
plot(1:neig, log10(fdaeig(1:neig)), 'bo-', ...
     1:neig, c(1)+ c(2)*(1:neig), 'b--')
xlabel('\fontsize{13} Eigenvalue Number')
ylabel('\fontsize{13} Log10 Eigenvalue')

% Figure 7.7 varimax rotation 

%  set up mean function

fdameanfd  = mean(fda_fd);
fdameanmat = eval_fd(fdatime, fdameanfd);

%  evaluate the harmonics

harmfd  = fdarotpcastr.harmfd;
harmmat = eval_fd(fdatime, harmfd);

fdapointtime = linspace(0,2300,201);
fdameanpoint = eval_fd(fdapointtime, fdameanfd);
harmpointmat = eval_fd(fdapointtime, harmfd);

fac = 0.1;
harmplusmat = zeros(201,3,2);
harmminsmat = zeros(201,3,2);
for j = 1:3
    harmplusmat(:,j,:) = fdameanpoint(:,1,:) + fac*harmpointmat(:,j,:);
    harmminsmat(:,j,:) = fdameanpoint(:,1,:) - fac*harmpointmat(:,j,:);
end

j=3;
plot(fdameanmat(:,1,1)-0.035,  fdameanmat(:,1,2), 'b-', ...
     harmplusmat(:,j,1)-0.035, harmplusmat(:,j,2), 'b:', ...
     harmminsmat(:,j,1)-0.035, harmminsmat(:,j,2), 'b:')
xlabel('\fontsize{13} ') 
ylabel('\fontsize{13} ')
axis([-0.075,0.075,-0.04,0.04])

j=2;
hold on
plot(fdameanmat(:,1,1)+0.035,  fdameanmat(:,1,2), 'b-', ...
     harmplusmat(:,j,1)+0.035, harmplusmat(:,j,2), 'b:', ...
     harmminsmat(:,j,1)+0.035, harmminsmat(:,j,2), 'b:')
hold off

%
% Section 7.5 Exploring Functional Covariation
%             with Canonical Correlation Analysis
%

%  set up temp_fd

lambda   = 1e2; 
fdParobj = fdPar(daybasis, harmaccelLfd, lambda);
temp_fd  = smooth_basis(time, tempav, fdParobj);
fdnames{1} = 'Day (July 2 to June 30)';
fdnames{2} = 'Weather Station';
fdnames{3} = 'Mean temperature (deg. C)';
temp_fd = putnames(temp_fd, fdnames);

ccaLfd = 2;
ccalambda = 5e6;
ccafdPar = fdPar(daybasis, ccaLfd, ccalambda);

nharm = 3;
ccastr  = cca_fd(temp_fd, logprec_fd, nharm, ccafdPar, ccafdPar);

ccawt_temp    = ccastr.wtfdx;
ccawt_logprec = ccastr.wtfdy;
corrs         = ccastr.corrs;

disp(corrs(1:3))

ccawtmat_temp    = eval_fd(time, ccawt_temp);
ccawtmat_logprec = eval_fd(time, ccawt_logprec);

%  Figure 7.8

phdl=plot(time, ccawtmat_temp(:,1),    'b-', ...
          time, ccawtmat_logprec(:,1), 'b--', ...
          yearRng, [0,0], 'b:');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Day (July 1 to June 30)')
ylabel('\fontsize{13} Canonical Weight Functions')
axis([yearRng, -0.12, 0.08])
legend('Temp.', 'Log Prec.', 'Location', 'SouthWest')

%  Figure 7.9

ccascr_temp    = ccastr.varx;
ccascr_logprec = ccastr.vary;

placeindex = [35,30,31,19,33,25,24,17,16,8,14,12,15,10,27,6,1,29];

phdl=plot(ccascr_temp(:,1), ccascr_logprec(:,1), '*');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} Temperature Canonical Weight') 
ylabel('\fontsize{13} Log Precipitation Canonical Weight')
for i = 1:length(placeindex)
    text(ccascr_temp(placeindex(i),1)+10, ...
         ccascr_logprec(placeindex(i),1), ...
         place(placeindex(i),:))
end
axis([-40,80,-5,8])

%
% Section 7.6 Details for the pca_fd and cca_fd Functions
%
help pca_fd
help cca_fd

%
% Section 7.7 Some Things to Try
%
% (exercises for the reader)

%
% Section 7.8 More to Read
%
