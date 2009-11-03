%  -----------------------------------------------------------------------
%                 Registered Handwriting Data
%  -----------------------------------------------------------------------

% These data are the X-Y coordinates of 20 replications of writing
% the script "fda".  The subject was Jim Ramsay.  Each replication
% is represented by 1401 coordinate values.  The scripts have been 
% extensively pre-processed.  They have been adjusted to a common
% length that corresponds to 2.3 seconds or 2300 milliseconds, and
% they have already been registered so that important features in
% each script are aligned.
% 
% This analysis is designed to illustrate techniques for working
% with functional data having rather high frequency variation and
% represented by thousands of data points per record.  Comments
% along the way explain the choices of analysis that were made.
% 
% The final result of the analysis is a third order linear 
% differential equation for each coordinate forced by a 
% constant and by time.  The equations are able to reconstruct
% the scripts to a fairly high level of accuracy, and are also
% able to accommodate a substantial amount of the variation in
% the observed scripts across replications.  by contrast, a 
% second order equation was found to be completely inadequate.
% 
% An interesting suprise in the results is the role placed by
% a 120 millisecond cycle such that sharp features such as cusps
% correspond closely to this period.  This 110-120 msec cycle
% seems is usually seen in human movement data involving rapid
% movements, such as speech, juggling and so on.

%  Last modified 26 July 2006

% Add paths to functional data functions and to the data

addpath ('c:\Program Files\matlab\fdaM')
addpath ('c:\Program Files\matlab\fdaM\examples\handwrit')

%  Input the data.  These 20 records have already been
%  normalized to a common time interval of 2300 milliseconds
%  and have beeen also registered so that prominent features
%  occur at the same times across replications.
%  Time will be measured in milliseconds and space in meters.
%  The data will require a small amount of smoothing, since
%  an error of 0.5 mm is characteristic of the OPTOTRAK 3D
%  measurement system used to collect the data.

fid      = fopen('fdareg.dat','rt');
fdaarray = reshape(fscanf(fid,'%f'), [20,2,1401]);
fdaarray = permute(fdaarray,[3,1,2]);
fdaarray = fdaarray/1000;   %  convert spatial unit to meters

%  Set up time values and range.
%  It is best to choose milliseconds as a time scale
%  in order to make the ratio of the time
%  unit to the inter-knot interval not too
%  far from one.  Otherwise, smoothing parameter values
%  may be extremely small or extremely large.

fdatime  = linspace(0, 2300, 1401)';
fdarange = [0, 2300];

%  The basis functions will be B-splines, with a spline 
%  placed at each knot.  One may question whether so many
%  basis functions are required, but this decision is found to 
%  be  essential for stable derivative estimation up to the 
%  third order at and near the boundaries.

%  Order 7 was used to get a smooth third derivative, which
%  requires penalizing the size of the 5th derivative, which
%  in turn requires an order of at least 7.
%  This implies norder + no. of interior knots = 1399 + 7 = 1406 
%  basis functions.  

nbasis   = 1406;
norder   =    7;
fdabasis = create_bspline_basis(fdarange, nbasis, norder);

%  The smoothing parameter value 1e8 was chosen to obtain a
%  fitting error of about 0.5 mm, the known error level in
%  the OPTOTRACK equipment.

fdafd  = fd(zeros(nbasis,20,2), fdabasis);
lambda = 1e9;
fdaPar = fdPar(fdafd, 5, lambda);

%  set up the functional data structure 

[fdafd, df, gcv] = smooth_basis(fdatime, fdaarray, fdaPar);
%  Add suitable names for the dimensions of the data.
fdafd_fdnames{1} = 'Milliseconds';
fdafd_fdnames{2} = 'Replications';
fdafd_fdnames{3} = 'Metres';
fdafd = putnames(fdafd, fdafd_fdnames);

%  display degrees of freedom and total GCV criterion

disp(['degrees of freedom = ',num2str(df)])  %  about 115
totalgcv = sum(sum(gcv));
disp(['total GCV = ',num2str(totalgcv)])  %  about 115
RMSgcv = sqrt(totalgcv)*1000; % about 0.3 mm
disp(['RMS GCV = ',num2str(RMSgcv)])  %  about 115

%  plot the fit to the data

plotfit_fd(fdaarray, fdatime, fdafd);

%  plot all curves

plot(fdafd)

%  compute values of curves and the values of the curve

fdameanfd  = mean(fdafd);
fdamat     = eval_fd(fdatime, fdafd);
fdameanmat = squeeze(eval_fd(fdatime, fdameanfd));

%  Set up motor control clock cycle times at every 
%  119 milliseconds. 

cyclelength = 117;
cycle  = 0:cyclelength:2300;
ncycle = length(cycle);

%  evaluate curves at cycle times

fdamatcycle     = eval_fd(cycle, fdafd);
fdameanmatcycle = squeeze(eval_fd(cycle, fdameanfd));

%  Indices of cycle times corresponding to important features:
%  -- the cusp in "f", 
%  -- the the cusp in "d", 
%  -- the first cusp in "a", 
%  -- the rest after the first cusp in "a", and 
%  -- the second cusp in "a".
%  It is remarkable that these features correspond so closely
%  with clock cycle times!

featureindex   = [3, 5, 7, 10, 13, 16, 19];
fdafeature     = fdamatcycle(featureindex,:,:);
fdameanfeature = fdameanmatcycle(featureindex,:);

%  Plot mean, including both sampling points and fit
%  Points at cycle times are plotted as blue circles, and
%  points at feature times are plotted as red circles.

subplot(1,1,1);
lhdl = plot(fdameanmat(:,1),      fdameanmat(:,2),      'b-', ...
            fdameanmatcycle(:,1), fdameanmatcycle(:,2), 'bo', ...
            fdameanfeature(:,1),  fdameanfeature(:,2),  'ro');
set(lhdl, 'LineWidth', 2);
xlabel('\fontsize{16} Metres')
ylabel('\fontsize{16} Metres')
axis([-.040, .040,-.040, .040]);
title('\fontsize{16} Mean script');

%  Plot individual curves, including both sampling points and fit
%  also plot the mean curve in the background.
%  Note how closely individual curve features are tied to the
%  feature cycle times.

subplot(1,1,1);
for i = 1:20
  lhdl = plot(fdamat(:,i,1),        fdamat(:,i,2),        'b-', ...
              fdamatcycle(:,i,1),   fdamatcycle(:,i,2),   'bo', ...
              fdafeature(:,i,1),    fdafeature(:,i,2),    'go', ...
              fdameanmat(:,1),      fdameanmat(:,2),      'r--', ...
              fdameanmatcycle(:,1), fdameanmatcycle(:,2), 'ro', ...
              fdameanfeature(:,1),  fdameanfeature(:,2),  'mo');
  set(lhdl, 'LineWidth', 2);
  xlabel('\fontsize{16} Metres')
  ylabel('\fontsize{16} Metres')
  axis([-.040, .040,-.040, .040]);
  title(['Record ', num2str(i)]);
  pause;
end

%  Evaluate the three derivatives and their means

D1fdamat = eval_fd(fdatime, fdafd, 1);
D2fdamat = eval_fd(fdatime, fdafd, 2);
D3fdamat = eval_fd(fdatime, fdafd, 3);

D1fdameanmat = eval_fd(fdatime, fdameanfd, 1);
D2fdameanmat = eval_fd(fdatime, fdameanfd, 2);
D3fdameanmat = eval_fd(fdatime, fdameanfd, 3);

%  Plot the individual acceleration records.
%  In these plots, acceleration is displayed as
%  metres per second per second.

%  Cycle and feature times are plotted as vertical
%  dashed lines, un-featured cycle times as red 
%  dotted lines, and cycle times of features as 
%  heavier magenta solid lines.

subplot(1,1,1);
for i=1:20
    lhdl = plot(fdatime, 1e6*D2fdamat(:,i,1), '-', ...
        fdatime, 1e6*D2fdamat(:,i,2), '-', ...
        [0,2300],[0,0], 'r:');
    set(lhdl, 'Linewidth', 2)
    hold on
    plotrange = [-12,12];
    for k=1:length(cycle)
        lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
        set(lhdl, 'LineWidth', 1')
    end
    for j=1:length(featureindex)
        k = featureindex(j);
        lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
        set(lhdl, 'LineWidth', 2')
    end
    hold off
    xlabel('\fontsize{16} Milliseconds')
    ylabel('\fontsize{16} Meters/msec/msec')
    title(['\fontsize{16} Curve ',num2str(i)])
    legend('X', 'Y')
    axis([0, 2300, -12, 12]);
    pause
end

%  Compute and plot the acceleration magnitudes,  
%  also called the tangential accelerations.  

D2mag  = sqrt(D2fdamat(:,:,1).^2 + D2fdamat(:,:,2).^2);
D2magmean = mean(D2mag,2);

subplot(1,1,1);
plot(fdatime, 1e6*D2mag)
hold on
plotrange = [0,12];
for k=1:length(cycle)
    lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
    set(lhdl, 'LineWidth', 1')
end
for j=1:length(featureindex)
    k = featureindex(j);
    lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
    set(lhdl, 'LineWidth', 2')
end
hold off
xlabel('\fontsize{16} Milliseconds')
ylabel('\fontsize{16} Metres/sec/sec')
title('\fontsize{16} Acceleration Magnitude');
axis([0,2300,0,12]);

%  Plot the mean acceleration magnitude as well as
%  those for each curve.
%  Note the two rest cycles, one in "d" and one in "a"

subplot(1,1,1);
lhdl = plot(fdatime, 1e6*D2magmean);
set(lhdl, 'LineWidth', 2')
hold on
plotrange = [0,8];
for k=1:length(cycle)
    lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
    set(lhdl, 'LineWidth', 1')
end
for j=1:length(featureindex)
    k = featureindex(j);
    lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
    set(lhdl, 'LineWidth', 2')
end
hold off
xlabel('\fontsize{16} Milliseconds')
ylabel('\fontsize{16} Metres/sec/sec')
title('\fontsize{16} Mean acceleration Magnitude');
axis([0,2300,0,8]);

%  Plot each individual acceleration magnitude, along
%  with the mean magnitude as a green dashed line

subplot(1,1,1)
for i=1:20
    lhdl = plot(fdatime, 1e6*D2mag(:,i), '-', ...
        fdatime, 1e6*D2magmean, '--');
    set(lhdl, 'LineWidth', 2)
    hold on
    plotrange = [0,12];
    for k=1:length(cycle)
        lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
        set(lhdl, 'LineWidth', 1')
    end
    for j=1:length(featureindex)
        k = featureindex(j);
        lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
        set(lhdl, 'LineWidth', 2')
    end
    hold off
    xlabel('\fontsize{16} Milliseconds')
    ylabel('\fontsize{16} Metres/sec/sec')
    title(['\fontsize{16} Script ',num2str(i)])
    axis([0,2300,0,12])
    pause
end

%  ------------------------------------------------------------
%                 Principal Differential Analysis 
%  A third order equation forced by a constant and time is
%          estimated for X and Y coordinates separately.
%  Forcing with constant and time is required to allow for 
%  an arbitrary origin and the left-right motion in X.
%  ------------------------------------------------------------

difeorder  = 3;  %  order of equation

%  set up the two forcing functions

%  constant forcing
constbasis = create_constant_basis(fdarange);
constfd    = fd(ones(1,20), constbasis);
ufdcell{1} = constfd;               
% time forcing
linbasis   = create_monomial_basis(fdarange, 2);
lincoef    = zeros(2,20);
lincoef(2,:) = 1;
ufdcell{2} = fd(lincoef, linbasis); 

%  set up the corresponding weight functions

constfd    = fd(1, constbasis);
constfdPar = fdPar(constfd);
awtcell{1} = constfdPar;
awtcell{2} = constfdPar;

%  Define two basis systems for the derivative weight
%  functions.  One for a background analysis used as a
%  baseline and using constant weight functions, and
%  another using 125 basis functions that will be used
%  to generate the equaations.

%  Set the number of basis functions to match the 
%  119 msec clock time estimated above, so that the
%  gap between knots is 11.9 msec.

delta   = 11.9;
breaks  = [0:delta:2300, 2300];
nbreaks = length(breaks);
nwbasis = nbreaks + 2;
wbasis  = create_bspline_basis(fdarange, nwbasis);

%  ------------------------------------------------------------
%                   Analysis for coordinate X
%  ------------------------------------------------------------

%  Define the variable with a reduced basis to save
%  computation time

xbasis = create_bspline_basis(fdarange, 245);
xfdcell{1} = data2fd(squeeze(fdamat(:,:,1)), fdatime, xbasis);

%  Set up the derivative weight functions

bfd     = fd(zeros(1,1), constbasis);
bfdPar  = fdPar(bfd, 1, 0);
bwtcell{1,1} = bfdPar;
bwtcell{1,2} = bfdPar;
bwtcell{1,3} = bfdPar;

%  carry out principal differential analysis

[bestwtcell, aestwtcell, resfdcell] = ...
    pda_fd(xfdcell, bwtcell, awtcell, ufdcell, difeorder);

%  evaluate forcing functions

resfd  = resfdcell{1};
resmat = eval_fd(fdatime, resfd);

MSY = mean(mean(resmat.^2));

%  Set the number of basis functions to 125,

bfd     = fd(zeros(nwbasis,1), wbasis);
bfdPar  = fdPar(bfd, 1, 0);
bwtcell{1,1} = bfdPar;
bwtcell{1,2} = bfdPar;
bwtcell{1,3} = bfdPar;

%  carry out principal differential analysis

[bestwtcell, aestwtcell, resfdcell] = ...
    pda_fd(xfdcell, bwtcell, awtcell, ufdcell, difeorder);

%  evaluate forcing functions

resfd  = resfdcell{1};
resmat = eval_fd(fdatime, resfd);

%  compute a squared multiple correlation measure of fit
%  MSY = mean(mean(resmat.^2)); %  Used only with constant basis
%  Uncomment this line when the constant basis is used, and
%  comment it out otherwise.

MSE = mean(mean(resmat.^2));
RSQ = (MSY-MSE)/MSY;

disp(['R-squared = ',num2str(RSQ)])

%  Plot the weight functions

for j=1:3
    subplot(3,1,j)
    plot(getfd(bestwtcell{1,j}));
    ylabel(['Weight function ',num2str(j-1)]);
end

%  Plot the first derivative weight function.

b2fdX = getfd(bestwtcell{1,2});
b2vecX = eval_fd(fdatime, b2fdX);
b2meanX = mean(b2vecX);

subplot(1,1,1)
lhdl = plot(fdatime, b2vecX, '-', ...
            [0, 2300], [b2meanX, b2meanX], 'g--');
set(lhdl, 'LineWidth', 2);
hold on
plotrange = [0,6e-3];
for k=1:length(cycle)
    lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
    set(lhdl, 'LineWidth', 1')
end
for j=1:length(featureindex)
    k = featureindex(j);
    lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
    set(lhdl, 'LineWidth', 2')
end
hold off
xlabel('\fontsize{16} Time (milliseconds)')
ylabel('\fontsize{16} \beta_1(t)')
axis([0, 2300, 0, 6e-3])

%  plot the weight coefficient for the second derivative

b3fdX = getfd(bestwtcell{1,3});
b3vecX = eval_fd(fdatime, b3fdX);
b3meanX = mean(b3vecX);

subplot(1,1,1)
lhdl = plot(fdatime, b3vecX, '-', ...
            [0, 2300], [b3meanX, b3meanX], 'g--');
set(lhdl, 'LineWidth', 2);
hold on
plotrange = [-.05,.05];
for k=1:length(cycle)
    lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
    set(lhdl, 'LineWidth', 1')
end
for j=1:length(featureindex)
    k = featureindex(j);
    lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
    set(lhdl, 'LineWidth', 2')
end
hold off
xlabel('\fontsize{16} Time (milliseconds)')
ylabel('\fontsize{16} \beta_2(t)')
axis([0, 2300, plotrange])

%  display coefficients for forcing weight functions

disp(getcoef(getfd(aestwtcell{1})))
disp(getcoef(getfd(aestwtcell{2})))
 
%  plot all forcing functions

subplot(1,1,1)
plot(fdatime, 1e9*resmat, '-', fdatime, 1e9*D3fdameanmat(:,1), 'r:')
xlabel('\fontsize{16} Milliseconds')
ylabel('\fontsize{16} Meters/sec/sec/sec')
axis([0,2300,-200,200])

%  plot the mean forcing function along with third deriv.

resmeanfd  = mean(resfd);
resmeanvec = eval_fd(fdatime, resmeanfd);

subplot(1,1,1)
plot(fdatime, 1e9*resmeanvec, 'b-', ...
     fdatime, 1e9*D3fdameanmat(:,1), 'r:')
xlabel('\fontsize{16} Milliseconds')
ylabel('\fontsize{16} Meters/sec/sec/sec')
axis([0,2300,-200,200])

% Define a functional data object for the 
%  three derivative weight functions

wcoef1 = getcoef(getfd(bestwtcell{1}));
wcoef2 = getcoef(getfd(bestwtcell{2}));
wcoef3 = getcoef(getfd(bestwtcell{3}));

wcoef  = [wcoef1, wcoef2, wcoef3];
wfd    = fd(wcoef,wbasis);

%  solve the equation

xstart = eye(3);
xstart(1,1) =   fdameanmat(1,1);
xstart(2,2) = D1fdameanmat(1,1);
xstart(3,3) = D2fdameanmat(1,1);

odeoptions = odeset('RelTol', 1e-6);

[tp1, xp1] = ode45(@derivs, fdatime, xstart(:,1), odeoptions, wfd);
[tp2, xp2] = ode45(@derivs, fdatime, xstart(:,2), odeoptions, wfd);
[tp3, xp3] = ode45(@derivs, fdatime, xstart(:,3), odeoptions, wfd);

%  plot the three solutions

umatx   = [xp1(:,1),xp2(:,1),xp3(:,1)];
Dumatx  = [xp1(:,2),xp2(:,2),xp3(:,2)];
D2umatx = [xp1(:,3),xp2(:,3),xp3(:,3)];

subplot(3,1,1)
plot(fdatime, umatx, '-', [0,2300], [0,0], 'r:');
title('Function')

subplot(3,1,2)
plot(fdatime, Dumatx, [0,2300], [0,0], 'r:')
title('First Derivative');

subplot(3,1,3)
plot(fdatime, D2umatx, [0,2300], [0,0], 'r:')
title('Second Derivative');

%  plot fit to each curve ... hit any key after each plot

subplot(1,1,1)
index  = 1:20;
fdamat = eval_fd(fdatime, fdafd);
zmat   = [ones(1401,1),fdatime-1150,umatx];
for i = index
   xhat = zmat * (zmat\fdamat(:,i,1));
   lhdl = plot(fdatime, xhat, '-', ...
               fdatime, fdamat(:,i,1), '--');
   set(lhdl, 'LineWidth', 2)
   title(['X curve ',num2str(i)])
   axis([0, 2300, -0.04, 0.04])
   legend('\fontsize{13} Fit', '\fontsize{13} Observed', ...
          'Location', 'NorthWest')
   pause;
end

%  ------------------------------------------------------------
%  analysis for Y: a third order equation forced by a 
%  constant and time
%  ------------------------------------------------------------

%  Define the variable

yfdcell{1} = fdafd(:,2);

%  Define the variable with a reduced basis to save
%  computation time

ybasis = xbasis;
yfdcell{1} = data2fd(squeeze(fdamat(:,:,2)), fdatime, ybasis);

%  A constant basis is used for a background level of error

bfd     = fd(zeros(1,1), constbasis);
bfdPar  = fdPar(bfd, 1, 0);
bwtcell{1,1} = bfdPar;
bwtcell{1,2} = bfdPar;
bwtcell{1,3} = bfdPar;

%  carry out principal differential analysis

[bestwtcell, aestwtcell, resfdcell] = ...
    pda_fd(yfdcell, bwtcell, awtcell, ufdcell, difeorder);

%  evaluate forcing functions

resfd  = resfdcell{1};
resmat = eval_fd(fdatime, resfd);

MSY = mean(mean(resmat.^2));

%  Use 125 basis functions.

bfd     = fd(zeros(nwbasis,1), wbasis);
bfdPar  = fdPar(bfd, 1, 0);
bwtcell{1,1} = bfdPar;
bwtcell{1,2} = bfdPar;
bwtcell{1,3} = bfdPar;

%  carry out principal differential analysis

[bestwtcell, aestwtcell, resfdcell] = ...
    pda_fd(yfdcell, bwtcell, awtcell, ufdcell, difeorder);

%  evaluate forcing functions

resfd  = resfdcell{1};
resmat = eval_fd(fdatime, resfd);

%  compute a squared multiple correlation measure of fit
%  MSY = mean(mean(resmat.^2)); %  Used only with constant basis
%  Uncomment this line when the constant basis is used, and
%  comment it out otherwise.

MSE = mean(mean(resmat.^2));
RSQ = (MSY-MSE)/MSY;

disp(['R-squared = ',num2str(RSQ)])

%  Plot the weight functions

for j=1:3
    subplot(3,1,j)
    plot(getfd(bestwtcell{1,j}));
    ylabel(['Weight function ',num2str(j-1)]);
end

%  Plot the second derivative weight, defining the period
%  of a harmonic oscillator.

b2fdY = getfd(bestwtcell{1,2});
b2vecY = eval_fd(fdatime, b2fdY);
b2meanY = mean(b2vecY);

subplot(1,1,1)
lhdl = plot(fdatime, b2vecY, '-', ...
            [0, 2300], [b2meanY, b2meanY], 'g--');
set(lhdl, 'LineWidth', 2);
hold on
plotrange = [0,6e-3];
for k=1:length(cycle)
    lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
    set(lhdl, 'LineWidth', 1')
end
for j=1:length(featureindex)
    k = featureindex(j);
    lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
    set(lhdl, 'LineWidth', 2')
end
hold off
axis([0, 2300, 0, 6e-3])

%  Plot the third derivative weight

b3fdY = getfd(bestwtcell{1,3});
b3vecY = eval_fd(fdatime, b3fdY);
b3meanY = mean(b3vecY);

subplot(1,1,1)
lhdl = plot(fdatime, b3vecY, '-', ...
            [0, 2300], [b3meanY, b3meanY], 'g--');
set(lhdl, 'LineWidth', 2);
hold on
plotrange = [-0.07,0.07];
for k=1:length(cycle)
    lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
    set(lhdl, 'LineWidth', 1')
end
for j=1:length(featureindex)
    k = featureindex(j);
    lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
    set(lhdl, 'LineWidth', 2')
end
hold off
axis([0, 2300, -0.07, 0.07])

%  display coefficients for forcing weight functions

disp(getcoef(getfd(aestwtcell{1})))
disp(getcoef(getfd(aestwtcell{2})))
 
%  plot forcing functions

subplot(1,1,1)
plot(fdatime, 1e9*resmat, '-', ...
     fdatime, 1e9*D3fdameanmat(:,2), 'r:')
xlabel('\fontsize{16} Milliseconds')
ylabel('\fontsize{16} Meters/sec/sec/sec')
axis([0,2300,-200,200])

%  plot the mean forcing function along with second deriv.

resmeanfd  = mean(resfd);
resmeanvec = eval_fd(fdatime, resmeanfd);

plot(fdatime, 1e9*resmeanvec, 'b-', ...
     fdatime, 1e9*D3fdameanmat(:,2), 'r:')
xlabel('\fontsize{16} Milliseconds')
ylabel('\fontsize{16} Meters/sec/sec/sec')
axis([0,2300,-200,200])

% Define a functional data object for the 
%  three derivative weight functions

wcoef1 = getcoef(getfd(bestwtcell{1}));
wcoef2 = getcoef(getfd(bestwtcell{2}));
wcoef3 = getcoef(getfd(bestwtcell{3}));

wcoef  = [wcoef1, wcoef2, wcoef3];
wfd    = fd(wcoef,wbasis);

%  Set up a linear differential operator.
%  This isn't used in these analyses.

fdaLfd = Lfd(difeorder, fd2cell(wfd));

%  solve equation

ystart = eye(3);
ystart(1,1) =   fdameanmat(1,2);
ystart(2,2) = D1fdameanmat(1,2);
ystart(3,3) = D2fdameanmat(1,2);

odeoptions = odeset('RelTol', 1e-6);

[tp1, yp1] = ode45(@derivs, fdatime, ystart(:,1), odeoptions, wfd);
[tp2, yp2] = ode45(@derivs, fdatime, ystart(:,2), odeoptions, wfd);
[tp3, yp3] = ode45(@derivs, fdatime, ystart(:,3), odeoptions, wfd);

%  plot the three solutions

umaty   = [yp1(:,1),yp2(:,1),yp3(:,1)];
Dumaty  = [yp1(:,2),yp2(:,2),yp3(:,2)];
D2umaty = [yp1(:,3),yp2(:,3),yp3(:,3)];

subplot(3,1,1)
plot(fdatime, umaty, '-', [0,2300], [0,0], 'r:');
title('Function')

subplot(3,1,2)
plot(fdatime, Dumaty, [0,2300], [0,0], 'r:')
title('First Derivative');

subplot(3,1,3)
plot(fdatime, D2umaty, [0,2300], [0,0], 'r:')
title('Second Derivative');

%  plot fit to each curve ... hit any key after each plot

subplot(1,1,1)
index  = 1:20;
fdamat = eval_fd(fdatime, fdafd);
zmat   = [ones(1401,1),fdatime-1150,umaty];
for i = index
   yhat = zmat * (zmat\fdamat(:,i,2));
   lhdl = plot(fdatime, yhat, '-', ...
               fdatime, fdamat(:,i,2), '--');
   set(lhdl, 'LineWidth', 2)
   title(['X curve ',num2str(i)])
   axis([0, 2300, -0.04, 0.04])
   legend('\fontsize{13} Fit', '\fontsize{13} Observed', ...
          'Location', 'SouthEast')
   pause;
end

%  plot fit to each script ... hit any key after each plot

subplot(1,1,1)
index  = 1:20;
fdamat = eval_fd(fdatime, fdafd);
zmatx   = [ones(1401,1),fdatime-1150,umatx];
zmaty   = [ones(1401,1),fdatime-1150,umaty];
for i = index
   xhat = zmatx * (zmatx\fdamat(:,i,1));
   yhat = zmaty * (zmaty\fdamat(:,i,2));
   lhdl = plot(xhat, yhat, '-', ...
               fdamat(:,i,1), fdamat(:,i,2), '--');
   set(lhdl, 'LineWidth', 2)
   xlabel('\fontsize{16} X')
   ylabel('\fontsize{16} Y')
   title(['\fontsize{16} Script ',num2str(i)])
   axis([-0.05,0.05,-0.05,0.05])
   legend('\fontsize{13} Fit', '\fontsize{13} Observed', ...
          'Location', 'SouthEast')
   pause;
end

%  plot the two weight functions for the second derivative

subplot(1,1,1)
lhdl = plot(fdatime, b2vecX, 'b-', fdatime, b2vecY, 'g-', ...
            [0, 2300], [b2meanX, b2meanX], 'b--', ...
            [0, 2300], [b2meanY, b2meanY], 'g--');
set(lhdl, 'LineWidth', 2);
hold on
plotrange = [0,6e-3];
for k=1:length(cycle)
    lhdl = plot([cycle(k),cycle(k)],plotrange,'r:');
    set(lhdl, 'LineWidth', 1')
end
for j=1:length(featureindex)
    k = featureindex(j);
    lhdl = plot([cycle(k),cycle(k)],plotrange,'m-');
    set(lhdl, 'LineWidth', 2')
end
hold off
axis([0, 2300, 0, 6e-3])

%  Conclusions:

%  For each coordinate, a single differential equation was
%  estimated to describe all 20 coordinate records 
%  simulataneously.
%  The differential equation solutions are readable.
%  Moreover, they accommodate a good deal of the curve-to-curve
%  variation in the scripts.  
%  The cusp first in the "a" is not handled as well as the other cusps.
%  This may be due to poor registration at this point, or to the fact
%  that this cusp has a lot of variation from one replication to another.

% Combined analysis; time-varying second order equation

bfd     = fd(zeros(nwbasis,1), wbasis);
bfdPar  = fdPar(bfd, 1, 0);
bwtcell = cell(2,2,2);
bwtcell(:) = {bfdPar};

yfdcell = cell(2,1);
yfdcell{1} = fdafd(:,1);
yfdcell{2} = fdafd(:,2);

%  carry out principal differential analysis

[bestwtcell, aestwtcell, resfdcell] = ...
    pda_fd(yfdcell, bwtcell, [], [], 2,501);



