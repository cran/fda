addpath ('c:\matlab\fdaM.dir')
addpath ('c:\matlab\fdaM.dir\examples\growth')

%  Last modified 11 January 2001

%  -----------------------------------------------------------------------
%                    Berkeley Growth Data
%  -----------------------------------------------------------------------

%  ------------------------  input the data  -----------------------

% male data
fid = fopen('hgtm.dat','rt');
hgtmmat = reshape(fscanf(fid,'%f'),[31,39]);
% female data
fid = fopen('hgtf.dat','rt');
hgtfmat = reshape(fscanf(fid,'%f'),[31,54]);

age = [ 1:0.25:2, 3:8, 8.5:0.5:18 ]';
nage = length(age);
ncasem  = size(hgtmmat,2);
ncasef  = size(hgtfmat,2);

%  --------------  Smooth the data nonmonotonically  --------------
%  This smooth uses the usual smoothing methods to smooth the data,
%  but is not guaranteed to produce a monotone fit.  This may not
%  matter much for the estimate of the height function, but it can
%  have much more serious consequences for the velocity and
%  accelerations.  See the monotone smoothing method below for a
%  better solution, but one with a much heavier calculation overhead.

%  -----------  create fd object   --------------------
%  We use b-spline basis functions of order 6
%  Knots are positioned at the ages of observation.
%  A B-spline basis with knots at age values and order 6 is used

%  -----------  Create fd objects   ----------------------------
%  A B-spline basis with knots at age values and order 6 is used

knots  = age;
norder = 6;
nbasis = length(knots) + norder - 2;
hgtbasis = create_bspline_basis([1,18], nbasis, norder, knots);

hgtmfd = data2fd(hgtmmat, age, hgtbasis, {'Age', 'Boys',  'cm'});

hgtffd = data2fd(hgtfmat, age, hgtbasis, {'Age', 'Girls', 'cm'});

%  --- Smooth these objects, penalizing the 4th derivative  --
%  This gives a smoother estimate of the acceleration functions

Lfd    = 4;
lambda = 1e-2;

hgtmfd = smooth(hgtmfd, lambda, Lfd);
hgtffd = smooth(hgtffd, lambda, Lfd);

%  plot data and smooth, residuals, velocity, and acceleration

agefine = linspace(1,18,101)';

%  Males:

hgtmfit = eval_fd(age,     hgtmfd);
hgtmhat = eval_fd(agefine, hgtmfd);
velmhat = eval_fd(agefine, hgtmfd, 1);
accmhat = eval_fd(agefine, hgtmfd, 2);

for i = 1:ncasem
  subplot(2,2,1)
  plot(age, hgtmmat(:,i), 'o', agefine, hgtmhat(:,i), '-')
  axis([1, 18, 60, 200])
  xlabel('Years')
  title(['Height for male ',num2str(i)])
  subplot(2,2,2)
  resi = hgtmmat(:,i) - hgtmfit(:,i);
  ind  = find(resi >= -.7 & resi <= .7);
  plot(age(ind), resi(ind), 'o-', [1,18], [0,0], '--')
  axis([0,18,-.7,.7]);
  xlabel('Years')
  title('Residuals')
  subplot(2,2,3)
  ind = find(velmhat(:,i) >= 0 & velmhat(:,i) <= 20);
  plot(agefine(ind), velmhat(ind,i), '-', [1,18], [0,0], '--')
  axis([1,18,0,20])
  xlabel('Years')
  title('Velocity')
  subplot(2,2,4)
  ind = find(accmhat(:,i) >= -6 & accmhat(:,i) <= 6);
  plot(agefine(ind), accmhat(ind,i), '-', [1,18], [0,0], '--')
  axis([1,18,-6,6])
  xlabel('Years')
  title('Acceleration')
  pause
end

% Females:

hgtffit = eval_fd(age,     hgtffd);
hgtfhat = eval_fd(agefine, hgtffd);
velfhat = eval_fd(agefine, hgtffd, 1);
accfhat = eval_fd(agefine, hgtffd, 2);

for i = 1:ncasem
  subplot(2,2,1)
  plot(age, hgtfmat(:,i), 'o', agefine, hgtfhat(:,i), '-')
  axis([1, 18, 60, 200])
  xlabel('Years')
  title(['Height for female ',num2str(i)])
  subplot(2,2,2)
  resi = hgtfmat(:,i) - hgtffit(:,i);
  ind  = find(resi >= -.7 & resi <= .7);
  plot(age(ind), resi(ind), 'o-', [1,18], [0,0], '--')
  axis([0,18,-.7,.7]);
  xlabel('Years')
  title('Residuals')
  subplot(2,2,3)
  ind = find(velfhat(:,i) >= 0 & velfhat(:,i) <= 20);
  plot(agefine(ind), velfhat(ind,i), '-', [1,18], [0,0], '--')
  axis([1,18,0,20])
  xlabel('Years')
  title('Velocity')
  subplot(2,2,4)
  ind = find(accfhat(:,i) >= -6 & accfhat(:,i) <= 6);
  plot(agefine(ind), accfhat(ind,i), '-', [1,18], [0,0], '--')
  axis([1,18,-6,6])
  xlabel('Years')
  title('Acceleration')
  pause
end 

%  -------  Compute monotone smooths of the data  -----------

%  These analyses use a function called smooth_monotone 
%  that fits the data with a function of the form
%                   f(x) = b_0 + b_1 D^{-1} exp W(x)
%     where  W  is a function defined over the same range as X,
%                 W + ln b_1 = log Df and w = D W = D^2f/Df.
%  The constant term b_0 in turn can be a linear combinations of covariates:
%                         b_0 = zmat * c.
%  The fitting criterion is penalized mean squared error:
%    PENSSE(lambda) = \sum [y_i - f(x_i)]^2 +
%                     \lambda * \int [L W(x)]^2 dx
%  where L is a linear differential operator defined in argument LFD.
%  The function W(x) is expanded by the basis in functional data object
%  These are best quality fits that I and my colleagues, notably
%  R. D. Bock, have been able to achieve to date.

%  ------  First set up a basis for monotone smooth   --------
%  We use b-spline basis functions of order 6
%  Knots are positioned at the ages of observation.

norder = 6;
nbasis = nage + norder - 2;
wbasis = create_bspline_basis([1,18], nbasis, norder, age');

%  starting values for coefficient

cvec0 = zeros(nbasis,1);  
Wfd0 = fd(cvec0, wbasis);

zmat  = ones(nage,1);
wgt   = zmat;

Lfd    = 3;           %  penalize curvature of acceleration
lambda = 10.^(-0.5);  %  smoothing parameter

% -----------------  Male data  --------------------

cvecm = zeros(nbasis, ncasem);
betam = zeros(2, ncasem);
RMSEm = zeros(1, ncasem);

for icase=1:ncasem
   hgt    = hgtmmat(:,icase);
   [Wfd, beta, Fstr, iternum, iterhist] = ...
       smooth_monotone(age, hgt, wgt, Wfd0, zmat, Lfd, lambda);
   cvecm(1:nbasis,icase) = getcoef(Wfd);
   betam(:, icase) = beta;
   hgthat = beta(1) + beta(2).*monfn(age, Wfd);
   RMSE   = sqrt(mean((hgt - hgthat).^2.*wgt)/mean(wgt));
   RMSEm(icase) = RMSE;
   fprintf('%5.f %5.f %12.5f %10.4f\n', [icase, iternum, Fstr.f, RMSE])
end

% -----------------  Female data  --------------------

cvecf = zeros(nbasis, ncasef);
betaf = zeros(2, ncasef);
RMSEf = zeros(1, ncasef);

for icase=1:ncasef
   hgt    = hgtfmat(:,icase);
   [Wfd, beta, Fstr, iternum, iterhist] = ...
       smooth_monotone(age, hgt, wgt, Wfd0, zmat, Lfd, lambda);
   cvecf(1:nbasis,icase) = getcoef(Wfd);
   betaf(:, icase) = beta;
   hgthat = beta(1) + beta(2).*monfn(age, Wfd);
   RMSE   = sqrt(mean((hgt - hgthat).^2.*wgt)/mean(wgt));
   RMSEf(icase) = RMSE;
   fprintf('%5.f %5.f %12.5f %10.4f\n', [icase, iternum, Fstr.f, RMSE])
end

%  -----------------------------------------------------------
%                 Registration of Velocity Curves            
%  -----------------------------------------------------------

%  set up warping function basis

nbasisw = 5;
norder  = 4;
basisw  = create_bspline_basis([1,18], nbasisw, norder);
coef0   = zeros(nbasisw,1);
Wfd0    = fd(coef0, basisw);

%  parameters for registerfd

Lfd      = 2;
lambda   = 1;
iterlim  = 10;
dbglev   = 1;
conv     = 1e-2;
crit     = 2;
periodic = 0;

%  -----------------------------------------------------------
%              register the velocity curves for the girls
%    The results for the monotone smooth above are used here  
%  -----------------------------------------------------------

Wfdf    = fd(cvecf,wbasis);
xfine   = linspace(1,18,101)';
hgtmat  = monfn(xfine, Wfdf);
hgtffd  = data2fd(hgtmat, xfine, wbasis);
Dhgtffd = deriv(hgtffd);
y0fd    = deriv(mean(hgtffd),1);
y0vec   = eval_fd(y0fd, xfine);

index = 1:ncasef;

%  with plotting

for icase = index
   disp(['Case ',num2str(icase)])
   yfd = Dhgtffd(icase);
   subplot(2,2,1)
   plot(xfine, eval(yfd,xfine), '-', xfine, y0vec, '--');
   axis('square')
   title(['Unreg. Female ',num2str(icase)])
   regstr = registerfd(y0fd, yfd, Wfd0, ...
                 Lfd, lambda, periodic, iterlim, dbglev, conv, crit);
   yregfd = regstr.regfd;
   subplot(2,2,2)
   plot(xfine, eval(yregfd,xfine), '-', xfine, y0vec, '--');
   axis('square')
   title(['Reg. Female ',num2str(icase)])
   Wfd = regstr.Wfd;
   subplot(2,2,3)
   warpvec = eval_mon(xfine, Wfd);
   warpvec = 1 + 17.*warpvec./warpvec(101);
   plot(xfine, warpvec, '-', xfine, xfine, '--');
   axis('square')
   title('Warping Fn.')
   subplot(2,2,4)
   plot(xfine, warpvec- xfine, '-', [1,18], [0,0], '--')
   axis('square')
   title('Deformation Fn.')
   pause
end

%  without plotting

cvecfreg = cvecf;
for icase = index
   disp(['Case ',num2str(icase)])
   yfd = Dhgtffd(icase);
   regstr = registerfd(y0fd, yfd, Wfd0, ...
                 Lfd, lambda, periodic, iterlim, dbglev, conv, crit);
   cvecfreg(:,icase) = getcoef(regstr.regfd);
end
              
%  -----------------------------------------------------------
%  register the velocity curves for the boys
%  -----------------------------------------------------------

Wfdm     = fd(cvecm,wbasis);
xfine    = linspace(1,18,101)';
hgtmat   = monfn(xfine, Wfdm);
hgtmfd   = data2fd(hgtmat, xfine, wbasis);
Dhgtmfd  = deriv(hgtmfd);
y0fd     = deriv(mean(hgtmfd),1);
y0vec    = eval(y0fd, xfine);

%  interactively, with plotting

index = 1:ncasem;

for icase = index
   disp(['Case ',num2str(icase)])
   yfd = Dhgtmfd(icase);
   subplot(2,2,1)
   plot(xfine, eval(yfd,xfine), '-', xfine, y0vec, '--');
   axis('square')
   title(['Unreg. Male ',num2str(icase)])
   regstr = registerfd(y0fd, yfd, Wfd0, ...
                 Lfd, lambda, periodic, iterlim, dbglev, conv, crit);
   yregfd = regstr.regfd;
   subplot(2,2,2)
   plot(xfine, eval(yregfd,xfine), '-', xfine, y0vec, '--');
   axis('square')
   title(['Reg. Male ',num2str(icase)])
   Wfd = regstr.Wfd;
   warpvec = eval_mon(xfine, Wfd);
   warpvec = 1 + 17.*warpvec./warpvec(101);
   subplot(2,2,3)
   plot(xfine, warpvec, '-', xfine, xfine, '--');
   axis('square')
   title('Warping Fn.')
   subplot(2,2,4)
   plot(xfine, warpvec- xfine, '-', [1,18], [0,0], '--')
   axis('square')
   title('Deformation Fn.')
   pause
end

%  without plotting

cvecmreg = cvecm;
for icase = index
   disp(['Case ',num2str(icase)])
   yfd = Dhgtmfd(icase);
   regstr = registerfd(y0fd, yfd, Wfd0, ...
                 Lfd, lambda, periodic, iterlim, dbglev, conv, crit);
   cvecmreg(:,icase) = getcoef(regstr.regfd);
end

%  ---------------------------------------------------------------------
%        Monotone smooth of short term height measurements
%  ---------------------------------------------------------------------

%  ---------------- input the data  ----------------------------------

clear;
fid = fopen('onechild.dat','rt');
temp = fscanf(fid,'%f');
n = 83;
data = reshape(temp, [n, 2]);
day  = data(:,1);
hgt  = data(:,2);
wgt  = ones(n,1);
zmat = wgt;


%  set up the basis

nbasis   = 43;
norder   = 4;
hgtbasis = create_bspline_basis([day(1), day(n)], nbasis, norder);

%  set up the functional data object for W = log Dh

cvec0 = zeros(nbasis,1);
Wfd   = fd(cvec0, hgtbasis);

%  set parameters for the monotone smooth

Lfd     = 2;
lambda  = 1e-1;

%  carry out the monotone smooth
 
[Wfd, beta, Fstr, iternum, iterhist] = ...
           smooth_monotone(day, hgt, wgt, Wfd, zmat, Lfd, lambda);

%  plot the function W = log Dh

subplot(1,1,1)
plot(Wfd);

%  plot the data plus smooth

dayfine  = linspace(day(1),day(n),151)';
yhat     = beta(1) + beta(2).*eval_mon(day, Wfd);
yhatfine = beta(1) + beta(2).*eval_mon(dayfine, Wfd);
plot(day, hgt, 'o', dayfine, yhatfine, '-')
xlabel('\fontsize{12} Day')
ylabel('\fontsize{12} Centimeters')

%  plot growth velocity

Dhgt = beta(2).*eval_mon(dayfine, Wfd, 1);
plot(dayfine, Dhgt)
xlabel('\fontsize{12} Days')
ylabel('\fontsize{12} Centimeters/day')

