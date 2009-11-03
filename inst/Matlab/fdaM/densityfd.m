function [Wfdobj, C, Fstr, iter, iterhist] = ...
    densityfd(x, WfdPar, conv, iterlim, active, dbglev)
%DENSITYFD estimates the density of a sample of scalar observations.

%  Arguments are:
%  X       ... data value array
%  WFDPAR  ... functional parameter object specifying the initial log
%              density, the linear differential operator used to smooth
%              smooth it, and the smoothing parameter.
%  CONV    ... convergence criterion
%  ITERLIM ... iteration limit for scoring iterations
%  ACTIVE  ... indices among 1:NBASIS of parameters to optimize
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFDOBJ   ... functional data basis object defining final density
%  C        ... normalizing constant for density p = exp(Wfdobj)/C
%  FSTR     ... Struct object containing
%               FSTR.F    ... final log likelihood
%               FSTR.NORM ... final norm of gradient
%  ITERNUM  ... Number of iterations
%  ITERHIST ... History of iterations

%  last modified 12 January 2004

if nargin < 2
    error('WFDPAR is not supplied.');
end

%  check WfdPar

if ~isa_fdPar(WfdPar) 
    if isa_fd(WfdPar) | isa_basis(WfdPar)
        WfdPar = fdPar(WfdPar);
    else
        error(['WFDPAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
    end
end

%  set up WFDOBJ

Wfdobj = getfd(WfdPar);

%  set up LFDOBJ

Lfdobj = getLfd(WfdPar);
Lfdobj = int2Lfd(Lfdobj);
nderiv = getnderiv(Lfdobj);

%  set up LAMBDA

lambda = getlambda(WfdPar);

%  set up BASIS

basis  = getbasis(Wfdobj);
nbasis = getnbasis(basis);
rangex = getbasisrange(basis);

[N,m] = size(x);

if m > 2 & N > 2
    error('Argument X must have either one or two columns.');
end

if (N == 1 | N == 2) & m > 1
    x = x';
    n = N;
    N = m;
    m = n;
end

if m == 1
    f = ones(N,1);
else
    f    = x(:,2);
    fsum = sum(f);
    f    = f./fsum;
    x    = x(:,1);
end

inrng = find(x >= rangex(1) & x <= rangex(2));
if (length(inrng) ~= N)
    warning('Some values in X out of range and not used.');
end

x     = x(inrng);
f     = f(inrng);
nobs  = length(x);

%  set some default arguments and constants

if nargin < 6, dbglev  = 1;        end
if nargin < 5, active  = 2:nbasis; end
if nargin < 4, iterlim = 20;       end
if nargin < 3, conv    = 1e-2;     end

%  set up some arrays

climit    = [-50,0;0,400]*ones(2,nbasis);
cvec0     = getcoef(Wfdobj);
hmat      = zeros(nbasis,nbasis);
inactive  = ones(1,nbasis);
inactive(active) = 0;
inactive  = find(inactive);
ninactive = length(inactive);
dbgwrd    = dbglev > 1;

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda.*eval_penalty(basis, Lfdobj);
end

%  evaluate log likelihood
%    and its derivatives with respect to these coefficients

[logl, Dlogl] = loglfnden(x, f, basis, cvec0);

%  compute initial badness of fit measures

Foldstr.f  = -logl;
gvec       = -Dlogl;
if lambda > 0
    gvec      = gvec      + 2.*(Kmat * cvec0);
    Foldstr.f = Foldstr.f + cvec0' * Kmat * cvec0;
end
if ninactive > 0, gvec(inactive) = 0; end
Foldstr.norm = sqrt(mean(gvec.^2));

%  compute the initial expected Hessian

hmat = Varfnden(x, basis, cvec0);
if lambda > 0
    hmat = hmat + 2.*Kmat;
end
if ninactive > 0
    hmat(inactive,:) = 0;
    hmat(:,inactive) = 0;
    hmat(inactive,inactive) = eye(ninactive);
end

%  evaluate the initial update vector for correcting the initial bmat

deltac   = -hmat\gvec;
cosangle = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Foldstr.f, -logl, Foldstr.norm];
if dbglev > 0
    fprintf('\nIteration  Criterion  Neg. Log L  Grad. Norm\n')
    fprintf('\n%5.f     %10.4f %10.4f %10.4f\n', status);
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:)  = status;

%  quit if ITERLIM == 0

if iterlim == 0
    Fstr = Foldstr;
    iterhist = iterhist(1,:);
    C = normalize_phi(basis, cvec0);
    return;
end

%  -------  Begin iterations  -----------

STEPMAX = 5;
MAXSTEP = 400;
trial   = 1;
cvec    = cvec0;
linemat = zeros(3,5);

for iter = 1:iterlim
    iternum = iternum + 1;
    %  take optimal stepsize
    dblwrd = [0,0]; limwrd = [0,0]; stpwrd = 0; ind = 0;
    %  compute slope
    Fstr = Foldstr;
    linemat(2,1) = sum(deltac.*gvec);
    %  normalize search direction vector
    sdg     = sqrt(sum(deltac.^2));
    deltac  = deltac./sdg;
    dgsum   = sum(deltac);
    linemat(2,1) = linemat(2,1)/sdg;
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        disp('Initial slope nonnegative.')
        ind = 3;
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -1e-5;
        if dbglev>1, disp('Initial slope too small'); end
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  load up initial search matrix 
    linemat(1,1:4) = 0;
    linemat(2,1:4) = linemat(2,1);
    linemat(3,1:4) = Foldstr.f;
    %  output initial results for stepsize 0
    stepiter  = 0;
    if dbglev>1
        fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
            [stepiter, linemat(:,1)']);
    end
    ips = 0;
    %  first step set to trial
    linemat(1,5)  = trial;
    %  Main iteration loop for linesrch
    for stepiter = 1:STEPMAX
        %  ensure that step does not go beyond limits on parameters
        limflg  = 0;
        %  check the step size
        [linemat(1,5),ind,limwrd] = ...
            stepchk(linemat(1,5), cvec, deltac, limwrd, ind, ...
            climit, active, dbgwrd);
        if linemat(1,5) <= 1e-9
            %  Current step size too small ... terminate
            Fstr    = Foldstr;
            cvecnew = cvec;
            gvecnew = gvec;
            if dbglev > 1
                fprintf('Stepsize too small:  %10.4f\n', linemat(1,5));
            end
            if limflg
                ind = 1;
            else
                ind = 4;
            end
            break;
        end
        cvecnew = cvec + linemat(1,5).*deltac;
        %  compute new function value and gradient
        [logl, Dlogl] = loglfnden(x, f, basis, cvecnew);
        Fstr.f  = -logl;
        gvecnew = -Dlogl;
        if lambda > 0
            gvecnew = gvecnew + 2.*Kmat * cvecnew;
            Fstr.f = Fstr.f + cvecnew' * Kmat * cvecnew;
        end
        if ninactive > 0, gvecnew(inactive) = 0;  end
        Fstr.norm = sqrt(mean(gvecnew.^2));
        linemat(3,5) = Fstr.f;
        %  compute new directional derivative
        linemat(2,5) = sum(deltac.*gvecnew);
        linemat(3,5) = Fstr.f;
        %  output current results
        if dbglev > 1
            fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
                [stepiter, linemat(:,5)']);
        end
        %  compute next step
        [linemat,ips,ind,dblwrd] = ...
            stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbgwrd);
        trial  = linemat(1,5);
        %  ind == 0 implies convergence
        if ind == 0 | ind == 5, break; end
        %  end iteration loop
    end
    
    %  update current parameter vectors
    
    cvec   = cvecnew;
    gvec   = gvecnew;
    Wfdobj = putcoef(Wfdobj, cvec);
    status = [iternum, Fstr.f, -logl, Fstr.norm];
    iterhist(iter+1,:) = status;
    if dbglev > 0
        fprintf('%5.f     %10.4f %10.4f %10.4f\n', status);
    end
    %  test for convergence
    if abs(Fstr.f-Foldstr.f) < conv
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  exit loop if convergence
    if Fstr.f >= Foldstr.f, break; end
    %  compute the new Hessian
    hmat = Varfnden(x, basis, cvec);
    if lambda > 0
        hmat = hmat + 2.*Kmat;
    end
    if ninactive > 0
        hmat(inactive,:) = 0;
        hmat(:,inactive) = 0;
        hmat(inactive,inactive) = eye(ninactive);
    end
    %  evaluate the update vector
    deltac    = -hmat\gvec;
    cosangle  = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));
    if cosangle < 0
        if dbglev > 1, disp('cos(angle) negative'); end
        deltac = -gvec;
    end
    Foldstr = Fstr;
end
C = normalize_phi(basis, cvec);

%  ---------------------------------------------------------------

function [logl, Dlogl] = loglfnden(x, f, basis, cvec)
N = length(x);
nbasis  = getnbasis(basis);
fmat    = f*ones(1,nbasis);
fsum    = sum(f);
nobs    = length(x);
oneobs  = ones(nobs,1);
phimat  = getbasismatrix(x, basis);
cval    = normalize_phi(basis, cvec);
logl    = sum((phimat * cvec) .* f - fsum*log(cval)/N);
EDW     = oneobs*expect_phi(basis, cvec, cval);
Dlogl   = sum((phimat - EDW).*fmat)';

%  ---------------------------------------------------------------

function  Varphi = Varfnden(x, basis, cvec)
nbasis  = getnbasis(basis);
nobs    = length(x);
oneobs  = ones(nobs,1);
cval    = normalize_phi(basis, cvec);
EDw     = oneobs*expect_phi(basis, cvec, cval);
EDwDwt  = nobs.*squeeze(expect_phiphit(basis, cvec, cval));
Varphi  = (EDwDwt - EDw'*EDw);

%  ---------------------------------------------------------------

function ss = normalize_phi(basisfd, cvec, JMAX, EPS)

%  Computes integrals of
%      p(x) = exp phi'(x) * cvec
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISFD  ...  Basis function object with basis functions phi.
%  CVEC ... coefficient vector defining density, of length NBASIS
%  MU   ... mean values to be subtracted from variates
%  SIGMA .. standard deviation to define u = (x - mu)/sigma
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place.  Multiply a standard interval
%           like (-5,5) by sigma to make it scale free
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  The integral of the function.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisfd), 'basis')
    error('First argument must be a basis function object.');
end

nbasis = getnbasis(basisfd);
oneb   = ones(1,nbasis);
rng    = getbasisrange(basisfd);

%  set default arguments

if nargin < 4, EPS = 1e-7;     end
if nargin < 3, JMAX = 15;      end

%  set up first iteration

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h = ones(JMAXP,1);
h(2) = 0.25;
%  matrix SMAT contains the history of discrete approximations to the integral
smat = zeros(JMAXP,1);
%  the first iteration uses just the endpoints
x  = rng;
nx = length(x);
ox = ones(nx,1);
fx = getbasismatrix(x, basisfd);
wx = fx * cvec;
wx(wx < -50) = -50;
px = exp(wx);
smat(1)  = width.*sum(px)./2;
tnm = 0.5;
j   = 1;

%  now iterate to convergence
for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if j == 2
        x = (rng(1) + rng(2))/2;
    else
        x = rng(1)+del/2 : del : rng(2);
    end
    fx = getbasismatrix(x, basisfd);
    wx = fx * cvec;
    wx(wx < -50) = -50;
    px = exp(wx);
    smat(j) = (smat(j-1) + width.*sum(px)./tnm)./2;
    if j >= 5
        ind = (j-4):j;
        [ss, dss] = polintarray(h(ind),smat(ind),0);
        if ~any(abs(dss) >= EPS.*max(abs(ss)))
            %  successful convergence
            return;
        end
    end
    smat(j+1) = smat(j);
    h(j+1)   = 0.25*h(j);
end
warning(['No convergence after ',num2str(JMAX),' steps in NORMALIZE.PHI'])

%  ---------------------------------------------------------------

function ss = expect_phi(basisfd, cvec, Cval, nderiv, ...
    rng, JMAX, EPS)
%  Computes expectations of basis functions with respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISFD  ...  A basis function
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density, of length NBASIS
%  CVAL ... normalizing constant defining density
%  MU   ... mean value to be subtracted from variates
%  SIGMA .. standard deviation to define u = (x - mu)/sigma
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place
%  NDERIV . order of derivative required for basis function expectation
%  UWRD ... If T, expectation is of (D PHI)*U
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A vector SS of length NBASIS of integrals of functions.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisfd),'basis')
    error('First argument must be a basis function object.');
end

nbasis = getnbasis(basisfd);
oneb   = ones(1,nbasis);
rangeval = getbasisrange(basisfd);

%  set default arguments

if nargin < 7, EPS = 1e-7;     end
if nargin < 6, JMAX = 15;      end
if nargin < 5, rng = rangeval; end
if nargin < 4, nderiv = 0;     end
if nargin < 3, Cval = 1;       end

%  set up first iteration

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h = ones(JMAXP,1);
h(2) = 0.25;
%  matrix SMAT contains the history of discrete approximations to the integral
smat = zeros(JMAXP,nbasis);
sumj = zeros(1,nbasis);
%  the first iteration uses just the endpoints
x  = rng';
nx = length(x);
ox = ones(nx);
fx = getbasismatrix(x, basisfd);
wx = fx * cvec;
wx(wx < -50) = -50;
px = exp(wx)./Cval;
if nderiv == 0
    Dfx = fx;
else
    Dfx = getbasismatrix(x, basisfd, 1);
end
sumj = Dfx' * px;
smat(1,:)  = width.*sumj'./2;
tnm = 0.5;
j   = 1;

%  now iterate to convergence

for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if j == 2
        x = (rng(1) + rng(2))/2;
    else
        x = (rng(1)+del/2 : del : rng(2))';
    end
    nx = length(x);
    fx = getbasismatrix(x, basisfd);
    wx = fx * cvec;
    wx(wx < -50) = -50;
    px = exp(wx)./Cval;
    if nderiv == 0
        Dfx = fx;
    else
        Dfx = getbasismatrix(x, basisfd, 1);
    end
    sumj = Dfx' * px;
    smat(j,:) = (smat(j-1,:) + width.*sumj'./tnm)./2;
    if j >= 5
        ind = (j-4):j;
        temp = squeeze(smat(ind,:));
        [ss, dss] = polintarray(h(ind),temp,0);
        if ~any(abs(dss) > EPS*max(abs(ss)))
            %  successful convergence
            return;
        end
    end
    smat(j+1,:) = smat(j,:);
    h(j+1) = 0.25*h(j);
end
warning(['No convergence after ',num2str(JMAX),' steps in EXPECT.PHI'])

%  ---------------------------------------------------------------

function ss = expect_phiphit(basisfd, cvec, Cval, nderiv1, nderiv2, ...
    rng, JMAX, EPS)

%  Computes expectations of cross product of basis functions with
%  respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISFD  ...  A basis function
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density
%  CVAL ... normalizing constant defining density
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A matrix of order NBASIS of integrals of functions.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisfd),'basis')
    error('First argument must be a basis function object.');
end

nbasis = getnbasis(basisfd);
oneb   = ones(1,nbasis);
rangeval = getbasisrange(basisfd);

%  set default arguments

if nargin < 8, EPS = 1e-7;       end
if nargin < 7, JMAX = 9;         end
if nargin < 6, rng = rangeval;   end
if nargin < 5, nderiv2 = 0;      end
if nargin < 4, nderiv1 = 0;      end
if nargin < 3, Cval = ones(n,1); end

%  set up first iteration

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h = ones(JMAXP,1);
h(2) = 0.25;
%  matrix SMAT contains the history of discrete approximations to the integral
smat = zeros([JMAXP,nbasis,nbasis]);
%  the first iteration uses just the endpoints
x  = rng';
nx = length(x);
fx = getbasismatrix(x, basisfd);
wx = fx * cvec;
wx(wx < -50) = -50;
px = exp(wx)./Cval;
if nderiv1 == 0
    Dfx1 = fx;
else
    Dfx1 = getbasismatrix(x, basisfd, 1);
end
if nderiv2 == 0
    Dfx2 = fx;
else
    Dfx2 = getbasismatrix(x, basisfd, 1);
end
sumj = Dfx1' * ((px * oneb) .* Dfx2);
smat(1,:,:)  = width.*sumj./2;
tnm = 0.5;
j   = 1;

%  now iterate to convergence
for j = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if j == 2
        x = (rng(1) + rng(2))/2;
    else
        x = (rng(1)+del/2 : del : rng(2))';
    end
    nx = length(x);
    fx = getbasismatrix(x, basisfd);
    wx = fx * cvec;
    wx(wx < -50) = -50;
    px = exp(wx)./Cval;
    if nderiv1 == 0
        Dfx1 = fx;
    else
        Dfx1 = getbasismatrix(x, basisfd, 1);
    end
    if nderiv2 == 0
        Dfx2 = fx;
    else
        Dfx2 = getbasismatrix(x, basisfd, 1);
    end
    sumj = Dfx1' * ((px * oneb) .* Dfx2);
    smat(j,:,:) = (squeeze(smat(j-1,:,:)) + width.*sumj./tnm)./2;
    if j >= 5
        ind = (j-4):j;
        temp = squeeze(smat(ind,:,:));
        [ss, dss] = polintarray(h(ind),temp,0);
        if ~any(abs(dss) > EPS.*max(max(abs(ss))))
            %  successful convergence
            return;
        end
    end
    smat(j+1,:,:) = smat(j,:,:);
    h(j+1) = 0.25*h(j);
end

%  ---------------------------------------------------------------

function [y,dy] = polintarray(xa, ya, x)
%  YA is an array with up to 4 dimensions
%     with 1st dim the same length same as the vector XA
n     = length(xa);
yadim = size(ya);
nydim = length(yadim);
if yadim(2) == 1, nydim = 1; end
if yadim(1) ~= n, error('First dimension of YA must match XA'); end
difx = xa - x;
absxmxa = abs(difx);
tmp = 1:n;
ns = min(tmp(absxmxa == min(absxmxa)));
cs = ya;
ds = ya;
if nydim == 1, y = ya(ns);  end
if nydim == 2, y = ya(ns,:);  end
if nydim == 3, y = ya(ns,:,:);  end
if nydim == 4, y = ya(ns,:,:,:);  end
ns = ns - 1;
for m = 1:(n-1)
    if nydim == 1
        for i = 1:(n-m)
            ho      = difx(i);
            hp      = difx(i+m);
            w       = (cs(i+1) - ds(i))./(ho - hp);
            ds(i) = hp.*w;
            cs(i) = ho.*w;
        end
        if 2*ns < n-m
            dy = cs(ns+1);
        else
            dy = ds(ns);
            ns = ns - 1;
        end
    end
    if nydim == 2
        for i = 1:(n-m)
            ho      = difx(i);
            hp      = difx(i+m);
            w       = (cs(i+1,:) - ds(i,:))./(ho - hp);
            ds(i,:) = hp.*w;
            cs(i,:) = ho.*w;
        end
        if 2*ns < n-m
            dy = cs(ns+1,:);
        else
            dy = ds(ns,:);
            ns = ns - 1;
        end
    end
    if nydim == 3
        for i = 1:(n-m)
            ho      = difx(i);
            hp      = difx(i+m);
            w       = (cs(i+1,:,:) - ds(i,:,:))./(ho - hp);
            ds(i,:,:) = hp.*w;
            cs(i,:,:) = ho.*w;
        end
        if 2*ns < n-m
            dy = cs(ns+1,:,:);
        else
            dy = ds(ns,:,:);
            ns = ns - 1;
        end
    end
    if nydim == 4
        for i = 1:(n-m)
            ho      = difx(i);
            hp      = difx(i+m);
            w       = (cs(i+1,:,:,:) - ds(i,:,:,:))./(ho - hp);
            ds(i,:,:,:) = hp.*w;
            cs(i,:,:,:) = ho.*w;
        end
        if 2*ns < n-m
            dy = cs(ns+1,:,:,:);
        else
            dy = ds(ns,:,:,:);
            ns = ns - 1;
        end
    end
    y = y + dy;
end
