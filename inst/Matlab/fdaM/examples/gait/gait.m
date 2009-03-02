addpath ('c:\Program Files\matlab\fdaM')

%  Last modified 26 July 2006

%  -----------------------------------------------------------------------
%                        Gait data
%  -----------------------------------------------------------------------

addpath ('c:\Program Files\matlab\fdaM\examples\gait')

%  -------------  input the data for the two measures  ---------------

fid  = fopen('hip.dat','rt');
hip  = reshape(fscanf(fid,'%f'), [20,39]);
fid  = fopen('knee.dat','rt');
knee = reshape(fscanf(fid,'%f'), [20,39]);

gaittime = linspace(0.025, 0.975, 20)';

gaitarray = zeros(20, 39, 2);
gaitarray(:,:,1) = hip;
gaitarray(:,:,2) = knee;

%  ---------------  set up the fourier basis  ------------------------

gaitbasis = create_fourier_basis([0,1], 21);

%  -----------  create the fd object (no smoothing)  -----------------

gaitfd = data2fd(gaitarray, gaittime,  gaitbasis);
gaitfd_fdnames{1} = 'Normalized time';
gaitfd_fdnames{2} = 'Boys';
gaitfd_fdnames{3} = 'Angle (deg.)';
gaitfd = putnames(gaitfd, gaitfd_fdnames);

% -----------  set up the harmonic acceleration operator  ----------

Lbasisobj    = create_constant_basis([0,1]);
Lcoef        = [0, (2*pi)^2, 0];
wfd          = fd(Lcoef, Lbasisobj);
wfdcell      = fd2cell(wfd);     % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

%   smooth the data a bit penalizing harmonic acceleration

%  set up functional parameter object to define amount of smoothingt

lambda   = 1e-9;
fdParobj = fdPar(gaitbasis, harmaccelLfd, lambda);

%  smooth the functional data object

gaitfd = smooth_fd(gaitfd, fdParobj);

%  plot the curves in two panels

subplot(2,1,1)
plot(gaitfd(:,1))
xlabel('')
title('\fontsize{12} Knee Angle')
subplot(2,1,2)
plot(gaitfd(:,2))
title('\fontsize{12} Hip Angle')
xlabel('')

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/gaitangles.ps'

%  --------  plot curves and their first derivatives  ----------------

%  plot each pair of curves interactively

subplot(1,1,1)
casenames = [];
varnames  = ['Knee angle'; 'Hip angle ';];
plotfit_fd(gaitarray, gaittime, gaitfd, casenames, varnames)

%  plot the residuals, sorting cases by residual sum of squares

residual = 1;
sortwrd  = 1;
plotfit_fd(gaitarray, gaittime, gaitfd, casenames, varnames, ...
           residual, sortwrd)

%  plot first derivative of all curves

plot(gaitfd, 1)

%  -----  plot curves as cycles  --------

subplot(1,1,1);
cycleplot(gaitfd, 0);

%  -----  plot the mean functions and their first two derivatives

gaitmeanfd = mean(gaitfd);

plot(gaitmeanfd, 0)

plot(gaitmeanfd, 1)

plot(gaitmeanfd, 2)

%  plot of gait cycle for FDA lecture

gaitvec     = squeeze(eval_fd(gaitfd(1,:), gaittime));
gaitmeanvec = squeeze(eval_fd(gaitmeanfd,  gaittime));
gaitlet = ['A', 'B', 'C', 'D', 'E'];
gaitind = [1,4,7,12,16];

subplot(1,1,1)
plot(gaitvec(:,1), gaitvec(:,2), '.-', ...
     gaitmeanvec(:,1), gaitmeanvec(:,2), '.--')
xlabel('\fontsize{12} Knee Angle')
xlabel('\fontsize{12} Hip Angle')
axis([0,50,0,80])
hold on
for i=1:5
    text(gaitvec(gaitind(i),1),     gaitvec(gaitind(i),2),     gaitlet(i))
    text(gaitmeanvec(gaitind(i),1), gaitmeanvec(gaitind(i),2), gaitlet(i))
end
hold off

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/gaitloop.ps'

% ---------------  do a PCA of gait data  -------------------------------

%  do the PCA with varimax rotation

nharm   = 4;
lambda  = 1e-9;
gaitfdPar = fdPar(gaitfd, harmaccelLfd, lambda);
gaitpcastr = pca_fd(gaitfd, nharm, gaitfdPar);
gaitpcastr = varmx_pca(gaitpcastr);

%  plot harmonics

subplot(1,1,1)
plot_pca(gaitpcastr);

%  plot eigenvalues

gaiteigvals = gaitpca.values;
x = ones(16,2);
x(:,2) = reshape((5:20),[16,1]);
y = log10(gaiteigvals(5:20));
c = x\y;
subplot(1,1,1)
plot(1:20,log10(gaiteigvals(1:20)),'-o', ...
     1:20, c(1)+ c(2).*(1:20), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

gaitharmfd  = gaitpca.harmfd;
gaitharmmat = eval_fd(gaittime, gaitharmfd);
gaitvarprop = gaitpca.varprop;
gaitmeanvec = squeeze(eval_fd(gaittime, gaitmeanfd));

con = 5.*ones(1,4);
for j=1:4
    subplot(2,2,j)
    yplus = gaitmeanvec + con(j).*squeeze(gaitharmmat(:,j,:));
    plot(gaitmeanvec(:,1),gaitmeanvec(:,2),'g.')
    hold on
    for i=1:20
        plot([gaitmeanvec(i,1),yplus(i,1)],...
             [gaitmeanvec(i,2),yplus(i,2)],'b-')
    end
    hold off
    xlabel('Hip Angle')
    ylabel('Knee Angle')
    title(['PC ',num2str(j),' (',num2str(round(gaitvarprop(j)*1000)/10),'%)'])
    axis([-20,60,0,80])
end

print -dpsc2 'c:/Myfiles/talks/fdacourse/figs/gaitpca.ps'

%  ------  do a canonical correlation analysis of knee-hip curves  ------

%  first penalize the second derivative to get the results in the book

ncan    = 3;
lambda  = 7e-4;
gaitfdPar  = fdPar(gaitfd, 2, lambda);

gaitccastr = cca_fd(gaitfd(:,1), gaitfd(:,2), ncan, gaitfdPar, gaitfdPar);

subplot(2,1,1)
plot(gaitccastr.wtfd1)
xlabel('')
title('\fontsize{12} Hip canonical weight functions')
subplot(2,1,2)
plot(gaitccastr.wtfd2)
title('\fontsize{12} Knee canonical weight functions')

gaitccastr.corrs(1:ncan)

%  now penalize the harmonic acceleration

lambda  = 1e-6;
gaitfdPar  = fdPar(gaitfd, harmaccelLfd, lambda);

gaitccastr = cca_fd(gaitfd(:,1), gaitfd(:,2), ncan, gaitfdPar, gaitfdPar);

subplot(2,1,1)
plot(gaitccastr.wtfd1)
xlabel('')
title('\fontsize{12} Hip canonical weight functions')
subplot(2,1,2)
plot(gaitccastr.wtfd2)
title('\fontsize{12} Knee canonical weight functions')

gaitccastr.corrs(1:ncan)

%  ----------  compute the variance and covariance functions  -------

gaitvarbifd = var(gaitfd);

gaitvararray = eval_bifd(gaitvarbifd, gaittime, gaittime);

subplot(2,3,1)
contour(gaitvararray(:,:,1,1))
title('Knee - Knee')

subplot(2,3,2)
contour(gaitvararray(:,:,1,2))
title('Knee - Hip')

subplot(2,3,3)
contour(gaitvararray(:,:,1,3))
title('Hip - Hip')

subplot(2,3,4)
surf(gaitvararray(:,:,1,1))
title('Knee - Knee')

subplot(2,3,5)
surf(gaitvararray(:,:,1,2))
title('Knee - Hip')

subplot(2,3,6)
surf(gaitvararray(:,:,1,3))
title('Hip - Hip')

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/gaitcorr.ps'

%  ----  register the first derivative of the gait data  

%  set up basis for warping function

nbasis = 7;
wbasis = create_fourier_basis([0,1],nbasis);
Lfdobj = 3;
lambda = 1e-3;
WfdPar = fdPar(wbasis, Lfdobj, lambda);

index  = 1:39;  %  curves to be registered

%  set up target for registration of first derivatives

Dgaitfd = deriv(gaitfd(index,:),1);
y0fd    = mean(Dgaitfd);
yfd     = Dgaitfd(index);
xfine   = linspace(0,1,101)';
ofine   = ones(101,1);
y0vec   = squeeze(eval_fd(y0fd, xfine));
yvec    = eval_fd(yfd, xfine);

%  set parameters for registerfd

periodic = 1;  %  data are periodic
crit     = 2;
conv     = 1e-3;
iterlim  = 50;

%  carry out the registration

[regfd, Wfd, shift] = registerfd(y0fd, yfd, WfdPar, periodic, ...
                                 crit, conv, iterlim);

%  compute registered function and warping function values

yregmat = eval_fd(xfine, regfd);
warpmat = monfn(xfine, Wfd);
warpmat = ofine*shift' + warpmat./(ofine*warpmat(101,:));

%  plot the registered gait functions

plot(regfd)

%  plot the knee and hip derivatives for each case:
%  blue  solid:   unregistered 
%  green dashed:  mean unregistered 
%  red   solid:   registered 

for i = index
   subplot(1,2,1)
   plot(xfine, yvec(:,i,1), '-', xfine, y0vec(:,1), '--', xfine, yregmat(:,i,1), '-');
   axis('square')
   title(['Knee derivative ',num2str(index(i))])
   subplot(1,2,2)
   plot(xfine, yvec(:,i,2), '-', ...
        xfine, y0vec(:,2), '--', ...
        xfine, yregmat(:,i,2), '-');
   axis('square')
   title(['Hip derivative ',num2str(index(i))])
   pause
end

%  Plot the warping functions and display shifts
%  Note case 4, for which knee remains way out of
%    of phase with target, and case 31 with a large shift.

subplot(1,1,1)
for i = index
   plot(xfine, warpmat(:,i),   '-',  ...
        xfine, xfine+shift(i), '--', ...
        [0,1], [0,1],          ':')
   axis([-.05,1.05,-.05,1.05])
   axis('square')
   title(['Case ',num2str(index(i)),' shift = ',num2str(shift(i))])
   pause
end

