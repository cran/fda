% Add paths to functional data functions and to the data

addpath ('c:\Program Files\matlab\fdaM')
addpath ('c:\Program Files\matlab\fdaM\examples\melanoma')

%  Last modified  26 July 2006

%  -----------------------------------------------------------------------
%                      Melanoma Incidence data
%  -----------------------------------------------------------------------

%  input the data  

tempmat = load('melanoma.dat');

year  = tempmat(:,2);
mela  = tempmat(:,3);
nyear = length(year);
rng   = [min(year),max(year)];

xmat1 = [ones(37,1),year];
coef1 = xmat1\mela;
melahat1 = xmat1*coef1;

plot(year, mela, 'o', year, melahat1, '--')

xmat2 = [ones(37,1),year,sin(0.65*year)];
coef2 = xmat2\mela;
melahat2 = xmat2*coef2;

plot(year, mela, 'ko', year, melahat1, 'k--', year, melahat2, 'k-')
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Cases of Melanoma per 100,000')

% print -dps2 'c:/MyFiles/P651/melanoma.ps'

lnmela = log(mela);

xmat1 = [ones(37,1),year];
coef1 = xmat1\lnmela;
lnmelahat1 = xmat1*coef1;
sse1 = sum((lnmela-lnmelahat1).^2);

plot(year, lnmela, 'o', year, lnmelahat1, '--')

xmat2 = [ones(37,1),year,sin(0.65*year)];
coef2 = xmat2\lnmela;
lnmelahat2 = xmat2*coef2;
sse2 = sum((lnmela-lnmelahat2).^2);

plot(year, lnmela, 'ko', year, lnmelahat1, 'k--', year, lnmelahat2, 'k-')
xlabel('\fontsize{16} Year')
ylabel('\fontsize{16} Log Cases of Melanoma per 100,000')

% print -dps2 'c:/MyFiles/P651/melanoma.ps'

%  -----------------------------------------------------------------------
%             smooth data using B-splines
%  -----------------------------------------------------------------------

%  set up the basis with a knot at every year

knots    = year';
nbasis   = nyear + 4;
norder   = 6;
basisobj = create_bspline_basis(rng, nbasis, norder, knots);

%  smooth the data by penalizing the second derivative

Lfdobj = 2;
lambda = 1;
melafdPar = fdPar(basisobj, Lfdobj, lambda);

melafd = smooth_basis(year, mela, melafdPar);

%  plot the data and the smooth

subplot(1,1,1)
plotfit_fd(mela, year, melafd)

%  plot the residuals

plotfit_fd(mela, year, melafd, [], [], 1, 1)

%  set up operator to remove sinusoid plus trend

omega  = 0.65;
Lbasis = create_constant_basis(rng);
Lcoef  = [0, 0, omega^2, 0];
wfd          = fd(Lcoef, Lbasis);
wfdcell      = fd2cell(wfd);     % convert the FD object to a cell object
melaLfd      = Lfd(4, wfdcell);  %  define the operator object
lambda       = 1e-2;
melafdPar    = fdPar(basisobj, melaLfd, lambda);

%  smooth the data

melafd = smooth_basis(year, mela, melafdPar);

%  plot the results

plotfit_fd(mela, year, melafd)






