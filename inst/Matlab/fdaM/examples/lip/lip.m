% Need 'matlab\fdaM' in the path, 
% e.g., by installing the 'fda' package for R and running 'fdaMatlabPath' 
% or by creating a copy where the following will work:  
% addpath ('c:\Program Files\matlab\fdaM')
% addpath ('c:\Program Files\matlab\fdaM\examples\lip')

%  Last modified 2008.06.15;  previously modified 24 July 2006

%  -----------------------------------------------------------------------
%                       Lip Movement Data
%  -----------------------------------------------------------------------

%%
% 0.  input the data  

fid = fopen('lip.dat','rt');
lipmat = reshape(fscanf(fid,'%f'), [51, 20]);
nobs = size(lipmat,2);  %  number of replications

liptime  = (0:0.007:0.35)';  %  sampling points for each curve
%%
% 1.  Create an 'fd' object 'lipfd'
%
%  ----------  set up the b-spline basis object  ------------
%       use order 6 splines so we can look at acceleration

lipbasis = create_bspline_basis([0,0.35], 31, 6);

%  -----  create the fd object  ---------

%  set up the functional parameter object

Lfdobj = int2Lfd(4);  %  penalize fourth derivative
lambda = 1e-12;       %  use light smoothing
lipfdPar = fdPar(lipbasis, Lfdobj, lambda);

%  carry out the smoothing
 
lipfd = smooth_basis(liptime, lipmat, lipfdPar);

%  add names to dimensions

lipfd_fdnames{1} = 'time (seconds)';
lipfd_fdnames{2} = 'Replications';
lipfd_fdnames{3} = 'mm';
lipfd = putnames(lipfd, lipfd_fdnames);

%  ------------  plot data and fit  ----------------

plotfit_fd(lipmat, liptime, lipfd)

%  ------------- plot residuals  ------------

plotfit_fd(lipmat, liptime, lipfd, [], [], 1, 1)

%  ---------  summarize and plot the functions and their accelerations  -----

subplot(2,1,1)
plot(lipfd);
title('Lip position')

subplot(2,1,2)
plot(lipfd, 1, 1, 2);
ylabel('mm/t/t')
title('Lip acceleration')
%%
% 2.  register the data 

% using the minimum and the elbow as landmarks  ------
%               manually identify these points in each curve

%  there are two landmarks, in addition to the beginning and end

nmarks   = 2;

% if the landmarks have already been identified and saved as below,
%  skip over the following to the "load landmarks" statement

%  set up the matrix of acceleration values for identifying landmarks

D2lipmat = eval_fd(liptime, lipfd, int2Lfd(2));

%  plot each acceleration curve, and click on the two
%    maxima, near t = .4, the other near t = .75

%  locate landmarks manually by clicking the minimum and elbow
%  for each curve

lipmarks = zeros(nobs,nmarks);
index    = zeros(nmarks,1);
subplot(1,1,1)
for i = 1:nobs
  plot(liptime, D2lipmat(:,i), 'o', [0,1], [0,0], ':')
  title(['Curve ',num2str(i)])
  for j = 1:nmarks
    [x y] = ginput(1);  %  input two clicks on points here
    index(j) = round(x*51);
  end
  lipmarks(i,:) = liptime(index)';
end

% save lipmarks lipmarks  % save lip marks in case you need them for future work

%  load the landmarks if they have already been identified

load lipmarks

lipmarks = lipmarks*0.35;
lipmeanmarks = mean(lipmarks);

%  -----------   register the curves  ---------------------

%  In the first registration, the monotone smoothing option 
%  is used to compute the warping function.  This is slower,
%  but gives better results because the warping functions
%  are smoother.  See  below for a smooth not using the 
%  monotone smoothing option.

%  create a basis object for the warping function
%  it has order 4 (piecewise cubic) and two interior knots
%  positioned at the mean landmark values since
%  NBASIS = NORDER + # interior knots

nbasis = 6;
norder = 4;
breaks = [0,lipmeanmarks,0.35];
warpbasis = create_bspline_basis([0,0.35], nbasis, norder, breaks);
WfdPar    = fdPar(warpbasis,int2Lfd(2),0);

%  plot the basis

subplot(1,1,1)
plot(warpbasis)  %  plot of B-spline basis functions
line([lipmeanmarks(1),lipmeanmarks(1)],[0,1])  % first knot
line([lipmeanmarks(2),lipmeanmarks(2)],[0,1])  % second knot

%  call landmark registration function to set up struct LMRKSTR

[lipregfd, lipwarpfd, Wfd] = ...
     landmarkreg(lipfd, lipmarks, lipmeanmarks, WfdPar,0,0);

%  plot Wfd, the functions defining the warping functions

plot(lipwarpfd)

%  plot unregistered and registered curves

subplot(1,2,1)
plot(lipfd)
title('Unregistered')

subplot(1,2,2)
plot(lipregfd)
title('Registered')

%  plot unregistered and registered accelerations

subplot(1,2,1)
plot(lipfd,int2Lfd(2),1,1)
title('Unregistered')

subplot(1,2,2)
plot(lipregfd,int2Lfd(2),1,1)
title('Registered')

%  plot warping functions

subplot(1,1,1)
plot(lipwarpfd)
title('Warping Functions')

%  plot deformation functions: warp(t) - t

lipwarpmat = eval_fd(lipwarpfd,liptime);
lipdefmat  = lipwarpmat - liptime*ones(1,nobs);
plot(liptime, lipdefmat, '-', [0,1], [0,0], ':')
title('Deformation Functions')

%  In this second registration, the relation between curve
%  landmarks and marker landmarks is simply interpolated.  
%  This is faster.  

%  Create a basis object for the warping function.
%  it has order 4 (piecewise linear) and no interior knots,
%  and is therefore simply a cubic polynomial.
%  NBASIS = NORDER + # interior knots
%  We are lucky to get away with this because no curves
%  result in a failure of monotonicity. 

nbasis = 4;
norder = 4;
breaks = [0,lipmeanmarks,0.35];
warpbasis = create_bspline_basis([0,0.35], nbasis, norder, breaks);
WfdPar = fdPar(warpbasis,int2Lfd(0),1e-1);

%  plot the basis

subplot(1,1,1)
plot(warpbasis)  %  plot of B-spline basis functions
line([lipmeanmarks(1),lipmeanmarks(1)],[0,0.35])  % first knot
line([lipmeanmarks(2),lipmeanmarks(2)],[0,0.35])  % second knot

%  call landmark registration function to set up struct LMRKSTR

[lipregfd, lipwarpfd, Wfd] = ...
     landmarkreg(lipfd, lipmarks, lipmeanmarks, WfdPar, 0, 0);

%  plot unregistered and registered curves

subplot(1,2,1)
plot(lipfd)
title('Unregistered')

subplot(1,2,2)
plot(lipregfd)
title('Registered')

%  plot unregistered and registered accelerations

subplot(1,2,1)
plot(lipfd,int2Lfd(2),1,1)
title('Unregistered')

subplot(1,2,2)
plot(lipregfd,int2Lfd(2),1,1)
title('Registered')

%  plot warping functions

subplot(1,1,1)
plot(lipwarpfd)
title('Warping Functions')

%  plot deformation functions: warp(t) - t

lipwarpmat = eval_fd(lipwarpfd,liptime);
lipdefmat  = lipwarpmat - liptime*ones(1,nobs);
plot(liptime, lipdefmat, '-', [0,1], [0,0], ':')
title('Deformation Functions')
%%
% 3.  principal component analysis --------------------------

nharm = 4;
lippcastr = pca_fd(lipfd, nharm);

%  plot unrotated harmonics

subplot(1,1,1)
plot_pca(lippcastr)

%  rotate harmonics

lippcastr = varmx_pca(lippcastr);

%  plot rotated harmonics

plot_pca(lippcastr)
  
%  plot log eigenvalues

lipeigvals = lippcastr.values;
plot(1:19,log10(lipeigvals(1:19)),'-o')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')
%%
% 4.  principal differential analysis  -------------------

%  compute the weight functions by the basis expansion method

difeorder  = 2;  %  order of equation

awtcell = {};
ufdcell = {};

nwbasis = 21;
wbasis = create_bspline_basis([0,0.35], nwbasis);
wcoef0 = zeros(nwbasis,1);
betafd = fd(wcoef0, wbasis); 
betafdPar = fdPar(betafd);
bwtcell{1,1} = betafdPar;
bwtcell{1,2} = betafdPar;

xfdcell{1} = lipfd;

%  carry out principal differential analysis

[bfdcell, resfdcell, afdcell] = ...
    pda_fd(xfdcell, bwtcell, awtcell, ufdcell, difeorder);

%  plot the weight functions

for j=1:2
    subplot(2,1,j)
    plot(getfd(bfdcell{1,j}));
    ylabel(['Weight function ',num2str(j-1)]);
end

% overlay plot

pda_overlay(bfdcell)

%  set up a linear differential operator 

wcoef  = [getcoef(getfd(bfdcell{1})), getcoef(getfd(bfdcell{2}))];
wfd    = fd(wcoef,wbasis);
bwtcell = fd2cell(wfd);
lipLfd = Lfd(difeorder, bwtcell);

%  compute forcing functions

force = eval_fd(liptime, lipfd, lipLfd);
%  plot the forcing functions for each curve
subplot(1,1,1)
plot(liptime, force);
% axis([0,0.35,-1e3,1e3]);

%  plot the mean forcing function along with second deriv.

forcemean = mean(force,2);
D2mean = eval_fd(mean(lipregfd),liptime,int2Lfd(2));

subplot(1,1,1)
plot(liptime, forcemean, '-', liptime, D2mean, ':')
axis([0,1,-700,700])

%  solve equation

ystart = eye(2);
[tp1, yp1] = ode45(@derivs, liptime, ystart(:,1), [], wfd);
[tp2, yp2] = ode45(@derivs, liptime, ystart(:,2), [], wfd);

%  plot the two solutions

umat = [yp1(:,1),yp2(:,1)];
subplot(2,1,1)
plot(liptime, umat, liptime, zeros(51,1), ':'), title('Function');
Dumat = [yp1(:,2),yp2(:,2)];
subplot(2,1,2)
plot(liptime, Dumat, liptime, zeros(51,1), ':'), title('Derivative');

%  plot fit to each curve ...  hit any key after each plot

index = 1:nobs;
lipmat   = eval_fd(liptime, lipfd);
D2lipmat = eval_fd(liptime, lipfd, 2);
forcemat = eval_fd(liptime, lipfd, lipLfd);
for i = index
   subplot(2,1,1)
   %  solid line is forcing function, dashed 2nd deriv.
   plot(liptime, forcemat(:,i),    '-',  ...
        liptime, D2lipmat(:,i), '--', ...
        [0,0.35], [0,0],   ':')
   axis([0,0.35,-10000,10000])
   title(['Record ',num2str(i),' Forcing Fn.'])
   xhat = umat * (umat\lipmat(:,i));
   subplot(2,1,2)
   %  solid is fit, dashed is actual
   plot(liptime, xhat, '-', liptime, lipmat(:,i), '--')
   title('Function')
   pause;
end

