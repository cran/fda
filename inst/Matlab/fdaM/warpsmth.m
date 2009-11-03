function [Wfd, Fstr, iternum, iterhist] = ...
          warpsmth(x, y, wt, Wfd, Lfdobj, lambda, ...
                   conv, iterlim, active, dbglev)
%WARPSMTH smooths the relationship of Y to X 
%  by fitting a monotone fn.  f(x) = b_0 + b_1 D^{-1} exp W(x)
%     where  W  is a function defined over the same range as X,
%  W + ln b_1 = log Df and w = D W = D^2f/Df.
%  b_0 and b_1 are chosen so that f(x_1) = y_1 and f(x_n) = y_n.
%  The fitting criterion is penalized mean squared error:
%    PENSSE(lambda) = \sum [y_i - f(x_i)]^2 +
%                     \lambda * \int [L W(x)]^2 dx
%  W(x) is expanded by the basis in functional data object Wfd.
%
%  Arguments are ...
%  X       ...  A vector of argument values.
%  Y       ...  A vector of function value.
%  WT      ...  a vector of weights.
%  WFD     ...  functional data object for W.  It's coefficient array
%               has a single column, and these are the starting values
%               for the iterative minimization of mean squared error.
%  LFDOBJ  ...  linear differential opr defining roughness penalty to
%               be applied to W.  This may be either a functional data
%               object defining a linear differential operator, or a 
%               nonnegative integer.  If the latter, it specifies the
%               order of derivative to be penalized.
%               Lfdobj = 1 by default, corresponding to L = D.
%  LAMBDA  ...  smoothing parameter determining the amount of penalty, 
%               0 by default.
%  CONV    ...  convergence criterion, 0.0001 by default
%  ITERLIM ...  maximum number of iterations, 20 by default
%  ACTIVE  ...  vector of 1's and 0's indicating which coefficients
%               are to be optimized (1) or remain fixed (0).  All values
%               are 1 by default, except that if a B-spline basis is used,
%               the first value is set to 0.
%  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
%               no output, if 1, output at each iteration, if higher, output
%               at each line search iteration. 1 by default.
%
%  Returns are:
%  WFD     ...  Functional data object for W.  It's coefficient vector
%               contains the optimized coefficients.
%  FSTR    ...  Structure containing final criterion value, gradient, and
%               gradient norm.
%  ITERNUM ...  The total number of iterations taken.
%  ITERHIST...  An array with number of rows equal to number of iterations,
%               and columns corresponding to it. number, function value,
%               and mean squared gradient.

%  Last modified 16 June 2007

%  ensure that X and Y are column vectors

if size(x,1) > 1 && size(x,2) > 1
    error('Argument X is not a vector.');
end
if size(y,1) > 1 && size(y,2) > 1
    error('Argument Y is not a vector.');
end
x = x(:);
y = y(:);

nobs     = length(x);              %  number of observations
cvec     = getcoef(Wfd);           %  initial coefficients
basis    = getbasis(Wfd);          %  basis for Wfd
nbasis   = getnbasis(basis);    %  no. basis functions
type     = getbasistype(basis); %  type of basis
rangeval = getbasisrange(basis);

%  initialize arguments that are not included

if nargin < 10
    dbglev = 0;
end
if nargin < 9
    if strcmp(type, 'bspline') || strcmp(wtype, 'fourier')
        active = 2:nbasis;
    else
        active = 1:nbasis;
    end
end
if nargin < 8, iterlim = 20;         end
if nargin < 7, conv    = 0.0001;     end
if nargin < 6, lambda  = 0;          end
if nargin < 5, Lfdobj  = int2Lfd(1); end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  check some arguments

if any(diff(x) <= 0)
    error('Values in X not strictly increasing.');
end
if x(1) < rangeval(1) || x(nobs) > rangeval(2)
    error('Values in X are out of bounds.');
end
if any(wt) < 0
    error('One or more weights are negative.');
end

%  set up some variables

wtroot  = sqrt(wt);
climit  = 100.0.*([-1; 1]*ones(1,nbasis));     
inactive  = ones(1,nbasis);
inactive(active) = 0;
inactive  = find(inactive);

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda*eval_penalty(basis, Lfdobj);
else
    Kmat  = zeros(nbasis,nbasis);
end

%  Compute initial function and gradient values

[Fstr, Dyhat] = fngrad(y, x, wt, Wfd, lambda, Kmat, inactive);

%  compute the initial expected Hessian

hessmat = hesscal(Dyhat, wtroot, lambda, Kmat, inactive);  
if dbglev > 2
    %  output eigenvalues of Hessian
    eigval = sort(eig(hessmat));
    fprintf('%10.4e ', eigval);
    fprintf('\n');
end

%  evaluate the initial line search direction vector
deltac = linesearch(Fstr, hessmat, dbglev);

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Fstr.f, Fstr.norm'];
if dbglev >= 1
    fprintf('\nIter.   PENSSE   Grad Length\n')
    fprintf('%3.f %10.4f %10.4f\n', status);
end
if dbglev > 2
    for i = 1:nbasis, fprintf('%10.4f%', cvec(i)); end
    fprintf('\n');
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:)  = status;
if iterlim == 0, return;  end

%  -------  Begin main iterations  -----------

MAXSTEPITER = 10;
MAXSTEP = 100;
trial   = 1;
reset   = 0;
linemat = zeros(3,5);
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
            stepchk(linemat(1,5), cvec, deltac, limwrd, ind, climit, dbgwrd);
        if ind == 1, break; end  % break of limit hit twice in a row
        if linemat(1,5) <= 1e-7 
            %  Current step size too small ... terminate
            if dbglev >= 2
                fprintf('Stepsize too small: %15.7f\n', linemat(1,5));
            end
            break;
        end
        %  compute new function value and gradient
        cvecnew = cvec + linemat(1,5).*deltac;  %  update coefficients
        Wfdnew  = putcoef(Wfd, cvecnew);   %  update function W
        [Fstr, Dyhat] = fngrad(y, x, wt, Wfdnew, lambda, Kmat, inactive);
        linemat(3,5) = Fstr.f;
        %  compute directional derivative
        linemat(2,5) = sum(deltac.*Fstr.grad);
        if dbglev >= 2
            fprintf('                 %3.f %10.4f %12.6f %12.6f\n', ...
                [stepiter, linemat(:,5)']);
        end
        %  compute next line search step, also testing for convergence
        [linemat, ips, ind, dblwrd] = ...
                              stepit(linemat, ips, dblwrd, MAXSTEP);
        trial  = linemat(1,5);
        if trial == MAXSTEP, break; end
        %  ind == 0 means convergence
        if ind == 0 || ind == 5, break; end
    end
    %  -------  End of line search iterations  -----------
    %  complete fit to data by LS regression of y on f
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
        cvec   = cvecold;
        Wfd    = putcoef(Wfd, cvecold);
        Fstr   = Foldstr;
        deltac = -Fstr.grad;
        if dbglev > 2
            for i = 1:nbasis, fprintf('%10.4f%', cvec(i)); end
            fprintf('\n');
        end
        if reset == 1
            %  This is the second time in a row that this
            %     has happened ...  quit
            if dbglev >= 2
                fprintf('Reset twice, terminating.\n');
            end
            return;
        else
            reset = 1;
        end
    else
        %  function value has not increased,  check for convergence
        if abs(Foldstr.f-Fstr.f) < conv
            break;
        end
        %  update old parameter vectors and fit structure
        cvecold = cvec;
        Foldstr = Fstr;
        %  update the expected Hessian
        hessmat = hesscal(Dyhat, wtroot, lambda, Kmat, inactive);
        if dbglev > 2
            eigval = sort(eig(hessmat));
            fprintf('%10.4e ', eigval);
            fprintf('\n');
        end
        %  update the line search direction vector
        deltac = linesearch(Fstr, hessmat, dbglev);
        reset = 0;
    end
    %  store iteration status
    status = [iternum, Fstr.f, Fstr.norm'];
    iterhist(iter+1,:) = status;
    if dbglev >= 1
        fprintf('%3.f %10.4f %10.4f\n', status);
    end
end

%  ----------------------------------------------------------------

function [deltac, cosangle] = linesearch(Fstr, hessmat, dbglev)
deltac = -hessmat\Fstr.grad;
cosangle  = -sum(Fstr.grad.*deltac)./sqrt(sum(Fstr.grad.^2)*sum(deltac.^2));
if dbglev >= 2
    fprintf('Cos(angle) = %8.4f\n', cosangle);
end
if sum(cosangle) < 1e-7 
    if dbglev >=2, fprintf('\nCosine of angle too small\n'); end
    deltac = -Fstr.grad;
end

%  ----------------------------------------------------------------

function [Fstr, Dyhat] = fngrad(y, x, wt, Wfd, lambda, Kmat, inactive)
nobs   = length(x);
basis  = getbasis(Wfd);
rangex = getbasisrange(basis);
width  = rangex(2) - rangex(1);
cvec   = getcoef(Wfd);
%  get unnormalized function and gradient values
f      = monfn(x, Wfd);
Dyhat  = mongrad(x, Wfd);
%  adjust functions and derivatives for normalization
fmax   = f(nobs);
Dymax  = Dyhat(nobs,:);
Dyhat  = width.*(fmax.*Dyhat - f*Dymax)./fmax^2;
f      = rangex(1) + width*f/fmax;
%  update residuals and function values
res    = y - f;
Fstr.f = mean(res.^2.*wt);
nbasis = size(Dyhat,2);
temp   = Dyhat.*(wt*ones(1,nbasis));
grad   = -2.*temp'*res./nobs;
Fstr.grad = grad;
if lambda > 0
    Fstr.grad = Fstr.grad + 2 .* Kmat * cvec;
    Fstr.f    = Fstr.f + cvec' * Kmat * cvec;
end
if length(inactive) > 0, Fstr.grad(inactive) = 0; end
Fstr.norm = sqrt(sum(Fstr.grad.^2));   %  gradient norm

%  ----------------------------------------------------------------

function hessmat = hesscal(Dyhat, wtroot, lambda, Kmat, inactive)
[nobs, nbasis] = size(Dyhat);
temp = Dyhat.*(wtroot*ones(1,nbasis));
hessmat = 2.*temp'*temp./nobs; 
%  adjust for penalty
if lambda > 0, hessmat = hessmat + 2.*Kmat; end
%  adjust for inactive coefficients
if length(inactive) > 0
    hessmat(inactive,:    ) = 0;
    hessmat(:    ,inactive) = 0;
    hessmat(inactive,inactive) = eye(sum(inactive));
end

