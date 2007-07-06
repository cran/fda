function [Wfdobj, C, hmat, Fstr, iter, iterhist] = ...
    density_fd(x, Wfd0Par, conv, iterlim, dbglev)
%DENSITY_FD estimates the density p(x) of a sample of scalar observations.
%  These observations may be one of two forms:
%   1.  a vector of observatons x_i 
%   2.  a two-column matrix, with the observations x_i in the 
%       first column, and frequencies f_i in the second.  
%   Option 1. corresponds to all f_i = 1.

%  Arguments are:
%  X       ... data value array, either a vector or a two-column
%              matrix.
%  WFDPAR  ... functional parameter object specifying the initial log
%              density, the linear differential operator used to smooth
%              smooth it, and the smoothing parameter.
%  CONV    ... convergence criterion
%  ITERLIM ... iteration limit for scoring iterations
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFDOBJ   ... functional data basis object defining final density
%  C        ... normalizing constant for density p = exp(Wfdobj)/C
%  FSTR     ... Struct object containing
%               FSTR.F    ... final log likelihood
%               FSTR.NORM ... final norm of gradient
%  ITERNUM  ... Number of iterations
%  ITERHIST ... History of iterations

%  To plot the density function or to evaluate it, evaluate WFDOBJ,
%  exponentiate the resulting vector, and then divide by the normalizing
%  constant C.

%  last modified 16 June 2007

if nargin < 2
    error('WFDPAR is not supplied.');
end

%  check Wfd0Par

if ~isa_fdPar(Wfd0Par) 
    if isa_fd(Wfd0Par) || isa_basis(Wfd0Par)
        Wfd0Par = fdPar(Wfd0Par);
    else
        error(['WFDPAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
    end
end

%  set up WFDOBJ

Wfdobj = getfd(Wfd0Par);

%  set up LFDOBJ

Lfdobj = getLfd(Wfd0Par);
Lfdobj = int2Lfd(Lfdobj);

%  set up LAMBDA

lambda = getlambda(Wfd0Par);

%  set up BASIS

basisobj = getbasis(Wfdobj);
nbasis   = getnbasis(basisobj);
rangex   = getbasisrange(basisobj);
active   = 1:nbasis;

[N,m] = size(x);

if m > 2 && N > 2
    error('Argument X must have either one or two columns.');
end

if (N == 1 || N == 2) && m > 1
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
    warning('Wid1:range', ...
            'Some values in X out of range and not used.');
end

x     = x(inrng);
f     = f(inrng);

%  set some default arguments and constants

if nargin < 5, dbglev  = 1;        end
if nargin < 4, iterlim = 50;       end
if nargin < 3, conv    = 1e-4;     end

%  set up some arrays

climit = [-50,0;0,400]*ones(2,nbasis);
cvec0  = getcoef(Wfdobj);
dbgwrd = dbglev > 1;

zeromat = zerobasis(nbasis);

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda.*eval_penalty(basisobj, Lfdobj);
end

%  evaluate log likelihood
%    and its derivatives with respect to these coefficients

[logl, Dlogl, Cval0, EDw0] = loglfnden(x, f, basisobj, cvec0);

%  compute initial badness of fit measures

Foldstr.f  = -logl;
gvec       = -Dlogl;
if lambda > 0
    gvec      = gvec      + 2.*(Kmat * cvec0);
    Foldstr.f = Foldstr.f + cvec0' * Kmat * cvec0;
end
Foldstr.norm = sqrt(mean(gvec.^2));
gvec0 = zeromat'*gvec;

%  compute the initial expected Hessian

hmat = Varfnden(x, basisobj, cvec0, Cval0, EDw0);
if lambda > 0
    hmat = hmat + 2.*Kmat;
end

hmat0 = zeromat'*hmat*zeromat;

%  evaluate the initial update vector for correcting the initial bmat

deltac0  = -hmat0\gvec0;
deltac   = zeromat*deltac0;

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
    C = normalize_phi(basisobj, cvec0);
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
    dblwrd = [0,0]; limwrd = [0,0]; ind = 0;
    %  compute slope
    Fstr = Foldstr;
    linemat(2,1) = sum(deltac.*gvec);
    %  normalize search direction vector
    sdg     = sqrt(sum(deltac.^2));
    deltac  = deltac./sdg;
    linemat(2,1) = linemat(2,1)/sdg;
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        disp('Initial slope nonnegative.')
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
            break;
        end
        cvecnew = cvec + linemat(1,5).*deltac;
        %  compute new function value and gradient
        [logl, Dlogl, Cval, EDw] = ...
            loglfnden(x, f, basisobj, cvecnew);
        Fstr.f  = -logl;
        gvecnew = -Dlogl;
        if lambda > 0
            gvecnew = gvecnew + 2.*Kmat * cvecnew;
            Fstr.f = Fstr.f + cvecnew' * Kmat * cvecnew;
        end
        gvecnew0  = zeromat'*gvecnew;
        Fstr.norm = sqrt(mean(gvecnew0.^2));
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
                             stepit(linemat, ips, dblwrd, MAXSTEP);
        trial  = linemat(1,5);
        %  ind == 0 implies convergence
        if ind == 0 || ind == 5, break; end
        %  end iteration loop
    end
    
    %  update current parameter vectors
    
    cvec   = cvecnew;
    gvec   = gvecnew;
    gvec0  = zeromat'*gvec;
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
    hmat = Varfnden(x, basisobj, cvec, Cval, EDw);
    if lambda > 0
        hmat = hmat + 2.*Kmat;
    end
    hmat0 = zeromat'*hmat*zeromat;
    %  evaluate the update vector
    deltac0   = -hmat0\gvec0;
    cosangle  = -gvec0'*deltac0/ ...
                     sqrt(sum(gvec0.^2)*sum(deltac0.^2));
    if cosangle < 0
        if dbglev > 1, disp('cos(angle) negative');  end
        deltac0 = -gvec0;
    end
    deltac = zeromat*deltac0;
    Foldstr = Fstr;
end
C = normalize_phi(basisobj, cvec);

%  ---------------------------------------------------------------

function [logl, Dlogl, Cval, EDw] = loglfnden(x, f, basisobj, cvec)
nbasis  = getnbasis(basisobj);
fmat    = f*ones(1,nbasis);
fsum    = sum(f);
nobs    = length(x);
oneobs  = ones(nobs,1);
phimat  = getbasismatrix(x, basisobj);
Cval    = normalize_phi(basisobj, cvec);
logl    = sum((phimat * cvec) .* f - fsum*log(Cval)/nobs);
EDw     = expect_phi(basisobj, cvec, Cval);
Dlogl   = sum((phimat - oneobs*EDw').*fmat)';

%  ---------------------------------------------------------------

function  Varphi = Varfnden(x, basisobj, cvec, Cval, EDw)
nobs    = length(x);
EDwDwt  = expect_phiphit(basisobj, cvec, Cval);
Varphi  = nobs.*(EDwDwt - EDw*EDw');

%  ---------------------------------------------------------------

function Cval = normalize_phi(basisobj, cvec)

%  Computes integrals of
%      p(x) = exp phi'(x) * cvec
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISOBJ ... Basis function object with basis functions phi.
%  CVEC ... coefficient vector defining density, of length NBASIS

%  Return:
%  The integral of the function.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj), 'basis')
    error('First argument must be a basis function object.');
end

Cval = funcint(@sumpxfun, cvec, basisobj);

%  ---------------------------------------------------------------

function Ephi = expect_phi(basisobj, cvec, Cval)
%  Computes expectations of basis functions with respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISOBJ ...  A basis function object.  
%  CVEC ... coefficient vector defining density, of length NBASIS
%  CVAL ... normalizing constant defining density

%  Return:
%  A vector SS of length NBASIS of integrals of functions.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj),'basis')
    error('First argument must be a basis function object.');
end

Ephi = funcint(@fxpxfun, cvec, basisobj)./Cval;

%  ---------------------------------------------------------------

function Ephiphi = expect_phiphit(basisobj, cvec, Cval)

%  Computes expectations of cross product of basis functions with
%  respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISOBJ ... A basis function object.  
%  CVEC ... coefficient vector defining density
%  CVAL ... normalizing constant defining density

%  Return:
%  A matrix of order NBASIS of integrals of functions.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj),'basis')
    error('First argument must be a basis function object.');
end

Ephiphi = funcint(@fxpxfxfun, cvec, basisobj)./Cval;

%  ---------------------------------------------------------------

function sumpx = sumpxfun(x, cvec, basisobj, nderiv1, nderiv2)
if nargin < 4, nderiv1 = 0; end
if nargin < 5, nderiv2 = 0; end
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px = exp(wx);
sumpx = sum(px);

%  ---------------------------------------------------------------

function fxpx = fxpxfun(x, cvec, basisobj, nderiv1, nderiv2)
if nargin < 4, nderiv1 = 0; end
if nargin < 5, nderiv2 = 0; end
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px   = exp(wx);
fxpx = fx' * px;

%  ---------------------------------------------------------------

function fxpxfx = fxpxfxfun(x, cvec, basisobj, nderiv1, nderiv2)
if nargin < 4, nderiv1 = 0; end
if nargin < 5, nderiv2 = 0; end
nbasis = getnbasis(basisobj);
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px   = exp(wx);
fxpxfx = fx'*((px*ones(1,nbasis)).*fx);

