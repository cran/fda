function [Wfdobj, beta, Cval, res, hmat, Fstr, iternum, iterhist] = ...
    LMdens_fd(y, WfdPar, zmat, beta0, sigma0, ...
              conv, iterlim, dbglev)
%  LMDENS_FD estimates a regression and the density of the residuals.
%  If betaO is [], or empty, then only the density is estimated. 

%  Arguments are:
%  Y       ... vector of dependent variable values
%  WFDPAR  ... functional parameter object specifying the initial log
%              density, the linear differential operator used to smooth
%              smooth it, and the smoothing parameter.
%  ZMAT    ... matrix of covariates
%  BETA0   ... initial model   coefficients
%  SIGMA0  ... initial standard error
%  CONV    ... convergence criterion
%  ITERLIM ... iteration limit for scoring iterations
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFD      ... functional data basis object defining final density
%  BETA     ... final model coefficients
%  FSTR     ... Struct object containing
%               FSTR.F    ... final log likelihood
%               FSTR.NORM ... final norm of gradient
%  ITERNUM  ... Number of iterations
%  ITERHIST ... History of iterations

%  last modified 16 June 2007

if nargin < 2
    error('WFDPAR is not supplied.');
end

%  check WfdPar

if ~isa_fdPar(WfdPar) 
    if isa_fd(WfdPar) || isa_basis(WfdPar)
        WfdPar = fdPar(WfdPar);
    else
        error(['WFDPAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
    end
end

%  set up WFDOBJ

Wfdobj = getfd(WfdPar);
cvec0  = getcoef(Wfdobj);
cvec0  = cvec0 - mean(cvec0);

%  set up LFDOBJ

Lfdobj = getLfd(WfdPar);
Lfdobj = int2Lfd(Lfdobj);

%  set up BASIS

basisobj = getbasis(Wfdobj);
nbasis   = getnbasis(basisobj);
rangey   = getbasisrange(basisobj);
active   = 1:nbasis;

nobs   = length(y);

%  set some default arguments

if nargin < 8, dbglev  = 1;     end
if nargin < 7, iterlim = 50;    end
if nargin < 6, conv    = 1e-4;  end

%  initialize some arrays

climit    = [-50,0;0,50]*ones(2,nbasis);

zeromat = zerobasis(nbasis);

ind1 = 1:nbasis;

%  Set up linear model

covwrd = ~isempty(beta0);
if covwrd
    %  Independent variables are present ... estimate linear model
    if size(zmat,1) ~= nobs 
        error('ZMAT must have as many rows as length(X)');
    end
    ncov = size(zmat,2);
    res0 = (y - zmat * beta0);
    ind2 = (nbasis+1):(nbasis+ncov);
    zeromat = [zeromat,              zeros(nbasis, ncov);
               zeros(ncov,nbasis-1), eye(ncov)          ];
else
    %  No independent variables are present
    zmat = [];
    res0 = y;
end

%  bring residuals out of range to range and set 
%    U0 = res./sigma0;

[res0, U0, indlo] = reschk(res0, rangey, sigma0);

%  initialize matrix Kmat defining penalty term

lambda = getlambda(WfdPar);
if lambda > 0 
    Kmat = lambda.*eval_penalty(basisobj, Lfdobj);
end

%  evaluate log likelihood
%    and its derivatives with respect to these coefficients

[logl, Dlogl, Cval, Ephi] = ...
               loglfnLM(basisobj, cvec0, U0, zmat, sigma0);

%  compute initial badness of fit measures

Foldstr.f  =  -logl;
gvec       = -Dlogl;
if lambda > 0 
    gvec(ind1) = gvec(ind1) + 2.*(Kmat * cvec0);
    Foldstr.f  = Foldstr.f  + cvec0' * Kmat * cvec0;
end
gvec0 = zeromat'*gvec;
Foldstr.norm = sqrt(mean(gvec0.^2));

%  compute the initial expected Hessian

hmat = -EHessfnLM(basisobj, cvec0, zmat, Cval, Ephi, sigma0);

% disp(hmat)

if lambda > 0 
    hmat(ind1,ind1) = hmat(ind1,ind1) + 2.*Kmat;
end

hmat0 = zeromat'*hmat*zeromat;

%  evaluate the initial update vector for correcting the initial bmat

deltac0 = -hmat0\gvec0;
deltac  = zeromat*deltac0;

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Foldstr.f, -logl, Foldstr.norm];
if dbglev > 0
    fprintf('\nIteration  Criterion  Neg. Log L  Grad. Norm\n')  
    fprintf('\n%5.f     %10.4f %10.4f %10.4f\n', status);
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:) = status;

%  quit if ITERLIM == 0

if iterlim == 0
    Fstr = Foldstr;
    iterhist = iterhist(1,:);
    beta = beta0;
    res  = res0;
    if length(indlo) > 0
        warning('Wid1:Ltrim', ...
            [num2str(length(indlo)),' lower residuals trimmed']);
    end
    indhi = find(res > rangey(2)*sigma0);
    if length(indhi) > 0
        warning('Wid2:Utrim', ...
            [num2str(length(indhi)),' upper residuals trimmed']);
    end
    return;
end

%  -------  Begin iterations  -----------

STEPMAX = 5;
MAXSTEP = 100;
trial   = 1;
cvec    = cvec0;
beta    = beta0;
res     = res0;
linemat = zeros(3,5);

for iter = 1:iterlim
    iternum = iternum + 1;
    Fstr    = Foldstr;
    %  set initial switches
    dblwrd = [0,0]; limwrd = [0,0]; ind = 0; ips = 0;
    %  normalize search direction vector
    sdg     = sqrt(sum(deltac.^2));
    deltac  = deltac./sdg;
    %  compute initial slope
    linemat(2,1) = sum(deltac.*gvec);
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        fprintf('Initial slope nonnegative.\n');
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -1e-5;
        if dbglev > 1, fprintf('Initial slope too small\n'); end
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  load up initial search matrix 
    linemat(1,1:4) = 0;
    linemat(2,1:4) = linemat(2,1);
    linemat(3,1:4) = Foldstr.f;
    %  output initial results for stepsize 0
    stepiter  = 0;
    if dbglev > 1
        fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
            [stepiter, linemat(:,1)']); 
    end
    %  first step set to trial
    linemat(1,5)  = trial;
    %  Main iteration loop for linesrch
    for stepiter = 1:STEPMAX
        %  ensure that step does not go beyond limits on parameters
        %  check the step size
        [linemat(1,5),ind,limwrd] = ...
            stepchk(linemat(1,5), cvec, deltac, limwrd, ind, ...
            climit, active, dbglev);
        if linemat(1,5) <= 1e-9 
            %  Current step size too small ... terminate
            Fstr    = Foldstr;
            cvecnew = cvec;
            betanew = beta;
            gvecnew = gvec;
            if dbglev > 1
                fprintf('Stepsize too small: %10.4f', linemat(1,5));
            end
            break;
        end
        cvecnew = cvec + linemat(1,5).*deltac(ind1);
        if covwrd
            betanew = beta + linemat(1,5).*deltac(ind2);
            resnew  = y - zmat * betanew;
        else
            resnew  = y;
        end
        %  compute new function value and gradient
        [resnew, Unew] = reschk(resnew, rangey, sigma0);
        [logl, Dlogl, Cval, Ephi]  = ...
              loglfnLM(basisobj, cvecnew, Unew, zmat, sigma0);
        Fstr.f  =  -logl;
        gvecnew = -Dlogl;
        if lambda > 0 
            gvecnew(ind1) = gvecnew(ind1) + 2.*Kmat * cvecnew;
            Fstr.f = Fstr.f + cvecnew' * Kmat * cvecnew;
        end
        gvecnew0  = zeromat'*gvecnew;
        Fstr.norm = sqrt(mean(gvecnew0.^2));
        %  update search matrix
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
    beta   = betanew;
    Wfdobj = putcoef(Wfdobj, cvec);
    if covwrd
        beta = betanew;
        res  = y - zmat * beta;
    else
        res = y;
    end
    %  check residuals and truncate if needed
    [res, U, indlo, indhi] = reschk(res, rangey, sigma0);
    %  update and output iteration status
    status = [iternum, Fstr.f, -logl, Fstr.norm];
    iterhist(iter+1,:) = status;
    fprintf('%5.f     %10.4f %10.4f %10.4f\n', status);
    %  test for convergence
    if abs(Fstr.f-Foldstr.f) < conv
        iterhist = iterhist(1:(iternum+1),:);
        if length(indlo) > 0
            warning('Wid1:Ltrim', ...
                [num2str(length(indlo)),' lower residuals trimmed']);
        end
        if length(indhi) > 0
            warning('Wid2:Utrim', ...
                [num2str(length(indhi)),' upper residuals trimmed']);
        end
        break;
    end
    %  exit loop if convergence
    if Fstr.f >= Foldstr.f,  
        break;  
    end
    %  compute the new Hessian
    hmat = EHessfnLM(basisobj, cvec, zmat, Cval, Ephi, sigma0);
    if lambda > 0
        hmat(ind1,ind1) = hmat(ind1,ind1) + 2.*Kmat;
    end
    hmat0 = zeromat'*hmat*zeromat;
    %  evaluate the update vector
    deltac0    = -hmat0\gvec0;
    cosangle  = -gvec0'*deltac0/ ...
                     sqrt(sum(gvec0.^2)*sum(deltac0.^2));
    if cosangle < 0
        if dbglev > 1, disp('cos(angle) negative');  end
        deltac0 = -gvec0;
    end
    deltac = zeromat*deltac0;
    Foldstr = Fstr;    
end

%  ------------------------------------------------------

function [logl, Dlogl, Cval, Ephi] = ...
               loglfnLM(basisobj, cvec, U, zmat, sigma)
%  U is a vector containing standardized residuals.
%  The first NBASIS elements in the gradient are the 
%  derivatives with respect to CVEC defining the density.
%  The next NCOV elements, if required, are the derivatives
%  with respect to BETA defining the linear model
nobs   = length(U);
%  CVEC derivatives
phimat = getbasismatrix(U, basisobj);
Cval   = normalize_phi(basisobj, cvec);
logl   = sum(phimat*cvec  - log(Cval));
Ephi   = expect_phi(basisobj, cvec, Cval);
Dlogl  = sum(phimat)' - nobs.*Ephi;
%  BETA derivatives
if ~isempty(zmat)
    Dphimat = getbasismatrix(U, basisobj, 1);
    DWvec   = Dphimat*cvec;
    Dlogl   = [Dlogl; -zmat'*DWvec./sigma];
else
    Dlogl = [];
end

%  ------------------------------------------------------

function  EHess = EHessfnLM(basisobj, cvec, zmat,  Cval, Ephi, sigma) 
%  The upper left order NBASIS sub-matrix contains
%  the Hessian with respect to CVEC defining the density.
%  The lower right order NCOV sub-matrix, if required, 
%  is the expected Hessian with respect to BETA defining 
%  the linear model.  The off-diagonal submatrices 
%  contain the cross derivatives
nobs    = size(zmat,1);
Ephiphi = expect_phiphit(basisobj, cvec, Cval);
EHess   = -nobs*(Ephiphi - Ephi*Ephi');
if ~isempty(zmat)
    onesn = ones(nobs,1);
%  This code for exact hessian
%     D1phimat = getbasismatrix(U, basisobj, 1);
%     D2phimat = getbasismatrix(U, basisobj, 2);
%     D2Wvec   = D2phimat * cvec;
%     DcdlnL   = -zmat'*D1phimat./sigma;
%     DddlnL   = (zmat.*(D2Wvec*ones(1,ncov)))'*zmat./sigma^2;
%     EHess    = [EHess, DcdlnL'; DcdlnL, DddlnL];
%  This code is for expected hessian
    EDphi   = expect_phi(basisobj, cvec, Cval, 1);
    EDcdlnL = -zmat'*onesn*EDphi'./sigma;
    ED2phi  = expect_phi(basisobj, cvec, Cval, 2);
    EDddlnL = zmat'*zmat.*(ED2phi'*cvec)./sigma^2;
    EHess   = [EHess, EDcdlnL'; EDcdlnL, EDddlnL];
end

%  ------------------------------------------------------

function  [res, U, indlo, indhi] = reschk(res, rangey, sigma0) 
% RESCHK brings residuals outside of limits to limits
%  First center the residuals.  See discussion of zmat
%  about the conditions that zmat must satisfy to make this
%  operation legitimate.
res = res - mean(res);
%  Look for residuals below lower boundary
indlo = find(res < rangey(1)*sigma0);
if length(indlo) > 0
    res(indlo) = rangey(1)*sigma0;
end
%  Look for residuals above upper boundary
indhi = find(res > rangey(2)*sigma0);
if length(indhi) > 0
    res(indhi) = rangey(2)*sigma0;
end
%  Compute standardized residuals
U   = res./sigma0;

%  ---------------------------------------------------------------

function Cval = normalize_phi(basisobj, cvec)

%  Computes integrals of
%      p(x) = exp phi'(x) * cvec
%  by numerical integration using Romberg integration

%  Arguments:
%  basisobj  ...  Basis function object with basis functions phi.
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

if ~strcmp(class(basisobj), 'basis')
    error('First argument must be a basis function object.');
end

Cval = funcint(@sumpxfun, cvec, basisobj, 0, 0);

%  ---------------------------------------------------------------

function Ephi = expect_phi(basisobj, cvec, Cval, nderiv)
%  Computes expectations of basis functions with respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  basisobj  ...  A basis function
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
%  UWRD ... If T, expectation is of (D PHI)*U
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A vector SS of length NBASIS of integrals of functions.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj),'basis')
    error('First argument must be a basis function object.');
end

if nargin < 4, nderiv = 0;  end

Ephi = funcint(@fxpxfun, cvec, basisobj, nderiv, 0)./Cval;

%  ---------------------------------------------------------------

function Ephiphi = expect_phiphit(basisobj, cvec, Cval)

%  Computes expectations of cross product of basis functions with
%  respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  basisobj  ...  A basis function
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

if ~strcmp(class(basisobj),'basis')
    error('First argument must be a basis function object.');
end

Ephiphi = funcint(@fxpxfxfun, cvec, basisobj, 0, 0)./Cval;

%  ---------------------------------------------------------------

function sumpx = sumpxfun(x, cvec, basisobj, nderiv1, nderiv2)
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px = exp(wx);
sumpx = sum(px);

%  ---------------------------------------------------------------

function fxpx = fxpxfun(x, cvec, basisobj, nderiv1, nderiv2)
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px   = exp(wx);
if nderiv1 == 0
    Dfx = fx;
else
    Dfx = getbasismatrix(x, basisobj, 1);
end
fxpx = Dfx' * px;

%  ---------------------------------------------------------------

function fxpxfx = fxpxfxfun(x, cvec, basisobj, nderiv1, nderiv2)
nbasis = getnbasis(basisobj);
phimat = full(eval_basis(x, basisobj));
wx = phimat * cvec;
wx(wx < -50) = -50;
px = exp(wx);
fxpxfx = (phimat.*(px*ones(1,nbasis)))' * phimat;



