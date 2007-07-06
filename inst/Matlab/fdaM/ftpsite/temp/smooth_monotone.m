function [Wfd, beta, Fstr, iternum, iterhist, y2cMap] = ...
          smooth_monotone(argvals, y, fdParobj, zmat, wt, ...
                          conv, iterlim, active, dbglev)
%SMOOTH_MONOTONE smooths the relationship of Y to X by fitting 
%     a monotone function of the form
%                   f(x) = b_0 + b_1 D^{-1} exp W(x)
%     where  W  is a function defined over the same range as X,
%  W + ln b_1 = log Df and w = D W = D^2f/Df.
%  The constant term b_0 in turn can be a linear combinations of covariates:
%     b_0 = zmat * c
%  The fitting criterion is penalized mean squared error:
%    PENSSE(lambda) = \sum w_i[y_i - f(x_i)]^2 +
%                     \lambda * \int [L W(x)]^2 dx
%  The function W(x) is expanded by the basis in functional data object
%    Wfd.   The coefficients of this expansion are called "coefficients"
%    in the comments, while the b's are called "regression coefficients"
%
%  Arguments are ...
%  ARGVALS ...  argument value array
%  Y       ...  function value array (the values to be fit)
%  FDPAROBJ ... A functional parameter or fdPar object.  This object 
%               contains the specifications for the functional data
%               object to be estimated by smoothing the data.  See
%               comment lines in function fdPar for details.
%               The functional data object WFD in FDPAROBJ is used
%               to initialize the optimization process.
%               It's coefficient array has a single column, and these 
%               are the starting values for the iterative minimization 
%               of mean squared error.
%               This argument may also be either a FD object, or a 
%               BASIS object.  In this case, the smoothing parameter 
%               LAMBDA is set to 0.
%  ZMAT    ...  a matrix of covariate values for the constant term.
%               It defaults to a column of one's;
%  WT      ...  a vector of weights
%  CONV    ...  convergence criterion, 0.0001 by default
%  ITERLIM ...  maximum number of iterations, 20 by default
%  ACTIVE  ... indices among 1:NBASIS of parameters to optimize 
%  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
%               no output, if 1, output at each iteration, if higher, output
%               at each line search iteration. 1 by default.
%
%  Returns are:
%  WFD     ...  Functional data object for W.  It's coefficient vector
%               contains the optimized coefficients.
%  BETA    ...  The regression coefficients b_0 and b_1.
%  FNEW    ...  The minimized penalized mean squared error.
%  MSG     ...  The mean of the squares of the components of the gradient
%  ITERNUM ...  The total number of iterations taken.
%  ITERHIST...  An array with number of rows equal to number of iterations,
%               and columns corresponding to it. number, function value,
%               mean squared gradient, and beta values.

%  Last modified 14 August 2006

if nargin < 3
    error('The first three arguments are not specified.');
end

%  check ARGVALS

if ~strcmp(class(argvals), 'double')
    error('ARGVALS is not of class double.');
end

if size(argvals,1) == 1
    argvals = argvals';
end

[nobs, ncl] = size(argvals);  %  number of observations

if ncl > 1
    error('ARGVALS is not a vector.')
end
if nobs < 3
    error('ARGVALS does not contain at least three values.');
end

%  check Y

if ~strcmp(class(y), 'double')
    error('Y is not of class double.');
end

if size(y,1) == 1
    y = y';
end

[nrw, ncl] = size(y);  %  number of observations

if ncl > 1
    error('Y is not a vector.')
end
if nrw ~= nobs
    error('Y is not the same length as ARGVALS.');
end

%  check fdParobj

if ~isa_fdPar(fdParobj) 
    if isa_fd(fdParobj) || isa_basis(fdParobj)
        fdParobj = fdPar(fdParobj);
    else
        error(['FDPAROBJ is not a functional parameter object, ', ...
               'not a functional data object, and ', ...
               'not a basis object.']);
    end
end

%  set up LFDOBJ

Lfdobj = getLfd(fdParobj);
Lfdobj = int2Lfd(Lfdobj);

%  set up LAMBDA

lambda = getlambda(fdParobj);

%  set up BASIS

Wfd      = getfd(fdParobj);
basisobj = getbasis(Wfd);

%  set up some intial arrays

oneobs = ones(nobs,1);           %  vector of one's
cvec   = getcoef(Wfd);           %  initial coefficients
nbasis = getnbasis(basisobj);    %  no. basis functions

%  initialize arguments that are not included

if nargin < 9, dbglev = 1;      end
if nargin < 8, active=1:nbasis; end
if nargin < 7, iterlim = 50;    end
if nargin < 6, conv = 0.0001;   end
if nargin < 5, wt = oneobs;     end
if nargin < 4, zmat  = oneobs;  end

%  set default for ZMAT

if isempty(zmat)
    zmat = oneobs;
end

%  check ZMAT

zdim = size(zmat);
if zdim(1) ~= nobs
    error('First dimension of ZMAT not correct.')
end

%  check WT

if any(wt) < 0
    error('One or more weights are negative.');
end

%  set up some variables

ncov    = zdim(2);                %  number of covariates
ncovp1  = ncov + 1;
wtroot  = sqrt(wt);
climit  = 100.0.*([-1; 1]*ones(1,nbasis));

%  set an indicator array for the fixed coefficients

inactive  = ones(1,nbasis);
inactive(active) = 0;
inactive  = find(inactive);

%  set up cell for storing basis function values

basiscell = cell(1,15);

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda*eval_penalty(basisobj, Lfdobj);
else
    Kmat  = zeros(nbasis,nbasis);
end

%  Compute initial function and gradient values

[Fstr, beta, Dyhat, basiscell] = ...
    fngrad(y, argvals, zmat, wt, Wfd, lambda, Kmat, basiscell, inactive);

%  compute the initial expected Hessian

hessmat = hesscal(beta, Dyhat, wtroot, lambda, Kmat, inactive);

%  evaluate the initial line search direction vector

deltac = linesearch(Fstr, hessmat, dbglev);

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Fstr.f, Fstr.norm, beta'];
if dbglev >= 1
    fprintf('\nIter.   PENSSE   Grad Length Intercept   Slope\n')
    fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
        [status(1:4),beta(ncovp1)]);
end
if dbglev > 2
    for i = 1:nbasis, fprintf('%10.4f%', cvec(i)); end
    fprintf('\n');
    for i = 1:3,      fprintf('%10.4f%', beta(i)); end
    fprintf('\n');
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:)  = status;
if iterlim == 0 
    y2cMap = [];
    return;  
end

%  -------  Begin main iterations  -----------

MAXSTEPITER = 10;
MAXSTEP = 100;
trial   = 1;
reset   = 0;
linemat = zeros(3,5);
betaold = beta;
cvecold = cvec;
Foldstr = Fstr;
dbgwrd  = dbglev >= 2;
%  ---------------  beginning of optimization loop  -----------
for iter = 1:iterlim
    iternum = iternum + 1;
    %  initialize logical variables controlling line search
    dblwrd = [0,0];  limwrd = [0,0];  ind = 0; ips = 0;
    %  compute slope at 0 for line search
    linemat(2,1) = sum(deltac.*Fstr.grad);
    %  normalize search direction vector
    sdg          = sqrt(sum(deltac.^2));
    deltac       = deltac./sdg;
    linemat(2,1) = linemat(2,1)/sdg;
    % initialize line search vectors
    linemat(:,1:4) = [0; linemat(2,1); Fstr.f]*ones(1,4);
    stepiter  = 0;
    if dbglev >= 2
        fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
            [stepiter, linemat(:,1)']);
    end
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        if dbglev >= 2, disp('Initial slope nonnegative.'); end
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -1e-7;
        if dbglev >= 2, disp('Initial slope too small'); end
        break;
    end
    %  first step set to trial
    linemat(1,5)  = trial;
    %  -------  Begin line search iterations  -----------
    cvecnew = cvec;
    Wfdnew  = Wfd;
    for stepiter = 1:MAXSTEPITER
        %  check that step size does not go beyond limits on parameters
        [linemat(1,5), ind, limwrd] = ...
            stepchk(linemat(1,5), cvec, deltac, limwrd, ind, ...
            climit, dbgwrd);
        if ind == 1, break; end  % break of limit hit twice in a row
        if linemat(1,5) <= 1e-7
            %  Current step size too small ... terminate
            if dbglev >= 2
                fprintf('Stepsize too small: %15.7f\n', linemat(1,5));
            end
        end
        %  compute new function value and gradient
        cvecnew = cvec + linemat(1,5).*deltac;  %  update coefficients
        Wfdnew  = putcoef(Wfd, cvecnew);   %  update function W
        [Fstr, beta, Dyhat, basiscell] = ...
            fngrad(y, argvals, zmat, wt, Wfdnew, lambda, ...
                   Kmat, basiscell, inactive);
        linemat(3,5) = Fstr.f;
        %  compute new directional derivative
        linemat(2,5) = sum(deltac.*Fstr.grad);
        if dbglev >= 2
            fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
                [stepiter, linemat(:,5)']);
        end
        %  compute next line search step, also testing for convergence
        [linemat, ips, ind, dblwrd] = ...
            stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbglev);
        trial  = linemat(1,5);
        if trial == MAXSTEP, break; end
        %  ind == 0 means convergence
        if ind == 0 || ind == 5, break; end
    end
    %  -------  End of line search iterations  -----------
    cvec = cvecnew;
    Wfd  = Wfdnew;
    %  check that function value has not increased
    if Fstr.f > Foldstr.f
        %  if it has, terminate iterations with a warning
        if dbglev >= 2
            fprintf('Criterion increased:');
            fprintf('%10.4f %10.4f\n',[Foldstr.f, Fstr.f]);
        end
        %  reset parameters and fit
        beta   = betaold;
        cvec   = cvecold;
        Wfd    = putcoef(Wfd, cvecold);
        Fstr   = Foldstr;
        deltac = -Fstr.grad;
        if dbglev > 2
            for i = 1:nbasis, fprintf('%10.4f%', cvec(i)); end
            fprintf('\n');
            for i = 1:3,      fprintf('%10.4f%', beta(i)); end
            fprintf('\n');
        end
        if reset == 1
            %  This is the second time in a row that this
            %     has happened ...  quit
            if dbglev >= 2
                fprintf('Reset twice, terminating.\n');
            end
            y2cMap = [];
            return;
        else
            reset = 1;
        end
    else
        %  function value has not increased,  check for convergence
        if abs(Foldstr.f-Fstr.f) < conv
            status = [iternum, Fstr.f, Fstr.norm, beta'];
            iterhist(iter+1,:) = status;
            if dbglev >= 1
                fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
                    [status(1:4),beta(ncovp1)]);
            end
            break;
        end
        %  update old parameter vectors and fit structure
        cvecold = cvec;
        betaold = beta;
        Foldstr = Fstr;
        %  update the expected Hessian
        hessmat = hesscal(beta, Dyhat, wtroot, lambda, Kmat, inactive);
        %  update the line search direction vector
        deltac = linesearch(Fstr, hessmat, dbglev);
        reset = 0;
    end
    %  store iteration status
    status = [iternum, Fstr.f, Fstr.norm, beta'];
    iterhist(iter+1,:) = status;
    if dbglev >= 1
        fprintf('%3.f %10.4f %10.4f %10.4f %10.4f\n', ...
            [status(1:4),beta(ncovp1)]);
    end
end
y2cMap = (inv(Dyhat' * Dyhat + lambda*Kmat)*Dyhat')./sqrt(nobs);

%  ----------------------------------------------------------------

function [deltac, cosangle] = linesearch(Fstr, hessmat, dbglev)
deltac = -hessmat\Fstr.grad;
%deltac = -symsolve(hessmat,Fstr.grad);
cosangle  = -sum(Fstr.grad.*deltac)./sqrt(sum(Fstr.grad.^2)*sum(deltac.^2));
if dbglev >= 2
    fprintf('Cos(angle) = %8.4f\n', cosangle);
end
if cosangle < 1e-7
    if dbglev >=2, fprintf('\nCosine of angle too small\n'); end
    deltac = -Fstr.grad;
end

%  ----------------------------------------------------------------

function [Fstr, beta, Dyhat, basiscell] = ...
    fngrad(y, argvals, zmat, wt, Wfd, lambda, Kmat, basiscell, inactive)
ncov   = size(zmat,2) + 1;
nobs   = length(argvals);
cvec   = getcoef(Wfd);
nbasis = size(cvec,1);
[f,     basiscell]  =   monfn(argvals, Wfd, basiscell);
[Dyhat, basiscell]  = mongrad(argvals, Wfd, basiscell);
xmat   = [zmat,f];
Dxmat  = zeros([nobs,ncov,nbasis]);
Dxmat(:,ncov,:) = Dyhat;
wtroot = sqrt(wt);
wtrtmt = wtroot*ones(1,ncov);
yroot  = y.*wtroot;
xroot  = xmat.*wtrtmt;
ninactive = length(inactive);
%  compute regression coefs.
beta = LSfit(yroot, zmat, f, wtrtmt);
%  update fitted values
yhat   = xmat*beta;
%  update residuals and function values
res    = y - yhat;
Fstr.f = mean(res.^2.*wt);
grad   = zeros(nbasis,1);
for j=1:nbasis
    Dxroot = squeeze(Dxmat(:,:,j)).*wtrtmt;
    yDx = yroot'*Dxroot*beta;
    xDx = xroot'*Dxroot;
    grad(j)    = beta'*(xDx+xDx')*beta - 2*yDx;
end
Fstr.grad = grad/nobs;
if lambda > 0
    Fstr.grad = Fstr.grad +    2 .* Kmat * cvec;
    Fstr.f    = Fstr.f    + cvec' * Kmat * cvec;
end
if ninactive > 0, Fstr.grad(inactive) = 0; end
Fstr.norm = sqrt(sum(Fstr.grad.^2));   %  gradient norm

%  ----------------------------------------------------------------

function hessmat = hesscal(beta, Dyhat, wtroot, lambda, Kmat, inactive)
nbet = length(beta);
[nobs, nbasis] = size(Dyhat);
temp = beta(nbet).*Dyhat;
temp = temp.*(wtroot*ones(1,nbasis));
hessmat = 2.*temp'*temp./nobs;
%  adjust for penalty
if lambda > 0, hessmat = hessmat + 2.*Kmat; end
%  adjust for inactive coefficients
ninactive = length(inactive);
if ninactive > 0
    hessmat(inactive,:    ) = 0;
    hessmat(:    ,inactive) = 0;
    hessmat(inactive,inactive) = eye(ninactive);
end

%  ----------------------------------------------------------------

function beta = LSfit(yroot, zmat, f, wtrtmt)
xmat  = [zmat,f];
ncol  = size(xmat,2);
xroot = xmat.*wtrtmt;
[Q, R, E] = qr(xroot);
tol       = size(xmat,1)*3e-16*abs(R(1,1));
y         = Q'*yroot;
for j=1:ncol
    if abs(R(j,j)) < tol
        if R(j,j) < 0, R(j,j) = -tol;  else R(j,j) = tol; end
    end
end
beta      = R\y;
beta      = E*beta;

