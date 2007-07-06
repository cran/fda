%  -----------------------------------------------------------------------
%                Pinch force data
%  -----------------------------------------------------------------------

%  Last modified 26 July 2006

addpath ('c:\Program Files\matlab\fdaM')
addpath ('c:\Program Files\matlab\fdaM\examples\pinch')

%  ------------------  input the data  --------------------

fid = fopen('pinch.dat','rt');
pinchvec = fscanf(fid,'%f');
pinchmat = reshape(pinchvec, [20,151])';

pinchtime  = linspace(0,150,151)'./600;

pinchbasis = create_bspline_basis([0,0.25], 153, 4);
pinchfdPar = fdPar(pinchbasis, 2, 1e-6);

%  -----------  create fd object (no smoothing)  --------------------

pinchfd = smooth_basis(pinchtime, pinchmat, pinchfdPar);
pinchfd_fdnames{1} = 'Seconds';
pinchfd_fdnames{2} = 'Replications';
pinchfd_fdnames{3} = 'Force (N)';
pinchfd = putnames(pinchfd, pinchfd_fdnames);

%  plot all curves

subplot(1,1,1)
plot(pinchfd)
title('Pinch Force Curves')

%  plot each curve along with the data

plotfit_fd(pinchmat, pinchtime, pinchfd)

%  plot the residuals, with cases sorted by size of mean squared residuals

casenames = [];
varnames  = [];
residual  = 1;
sortwrd   = 1;

plotfit_fd(pinchmat, pinchtime, pinchfd, casenames, varnames, ...
           residual, sortwrd)

%  ---------------------  do a PCA (light smoothing)  --------------

lambda = 1e-4;
pinchfdPar = fdPar(pinchbasis, 2, lambda);

nharm  = 3;
pinchpcastr = pca_fd(pinchfd, nharm, pinchfdPar);

plot_pca(pinchpcastr)

pincheigvals = pinchpcastr.values;
x = ones(17,2);
x(:,2) = reshape((3:19),[17,1]);
y = log10(pincheigvals(3:19));
c = x\y;
subplot(1,1,1)
plot(1:19,log10(pincheigvals(1:19)),'-o', ...
     1:19, c(1)+ c(2).*(1:19), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

