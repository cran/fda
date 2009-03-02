function [gval, basiscell, tval] = ...
                      mongrad(x, Wfd, basiscell, EPS, JMIN, JMAX)
% MONGRAD evaluates the gradient of a single monotone function of the form
%             h(x) = (D^{-1} exp Wfd)(x)
%  where  D^{-1} means taking the indefinite integral.
%  The interval over which the integration takes places is defined in
%  the basis object in Wfd.
%  A major overhead in this function and MONFN is the evaluation
%  of the basis functions at the dyadic argument sequence required for
%  the trapezoidal rule integration.  This is reduced by storing these
%  values in the cells in cell array BASISCELL, each cell corresponding
%  to an iteration of the trapezoidal rule, with a maximum of 15.  

%  Arguments:
%  X         ...  Vector of argument values at which gradient is evaluated
%  WFD       ...  Functional data object defining monotone function
%  BASISCELL ...  A cell array object of length 15
%                 containing basis function values.
%  EPS       ...  Relative error needed for convergence
%  JMIN      ...  Minimum number of step halving steps
%  JMAX      ...  Maximum number of step halving steps
%
%  Returns:
%  GVAL  ... values of derivatives in NBASIS cols
%  TVAL  ... Arguments used for trapezoidal approximation to integral
%  FVAL  ... Values of exp Wfd corresponding to TVAL

%  Last modified 14 August 2006

if nargin < 2,  error('There are less than two arguments');  end

%  set some constants

if nargin < 6,  EPS  = 1e-5;  end
if nargin < 5,  JMIN = 11;    end
if nargin < 4,  JMAX = 15;    end
if JMIN   < 5,  JMIN =  5;    end

%  get coefficient matrix and check it

coef  = getcoef(Wfd);
coefd = size(coef);
ndim  = length(coefd);
if ndim > 1 && coefd(2) ~= 1 
    error('WFD is not a single function');
end

basis  = getbasis(Wfd);
rng    = getbasisrange(basis);
nbasis = getnbasis(basis);
onebas = ones(1,nbasis);
width  = rng(2) - rng(1);

%  set up first iteration

JMAXP = JMAX + 1;
h     = ones(JMAXP,1);
h(2)  = 0.25;
%  matrix SMAT contains the history of discrete approximations to the
%    integral
smat = zeros(JMAXP,nbasis);
%  array TVAL contains the argument values used = the approximation
%  array FVAL contains the integral values at these argument values,
%     rows corresponding to argument values
%  the first iteration uses just the endpoints
iter  = 1;
xiter = rng';
tval  = xiter;
if nargin == 3
    if isempty(basiscell{iter})
        bmat = eval_basis(basis, xiter);
        basiscell{iter} = bmat;
    else
        bmat = basiscell{iter};
    end
else
    bmat = eval_basis(basis, xiter);
end
fx   = exp(bmat*coef);
grad = (fx * onebas).*bmat;
fval = grad;
smat(1,:)  = width*sum(grad)./2;
tnm = 0.5;
%  now iterate to convergence
nx = 1;
for iter = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if iter == 2
        xiter = (rng(1)+rng(2))/2;
    else
        xiter = linspace(rng(1)+del/2, rng(2)-del/2, nx)';
    end
    tval = [tval; xiter];
    if nargin == 3
        if isempty(basiscell{iter})
            bmat = eval_basis(basis, xiter);
            basiscell{iter} = bmat;
        else
            bmat = basiscell{iter};
        end
    else
        bmat = eval_basis(basis, xiter);
    end
    fx   = exp(bmat*coef);
    grad = (fx*onebas).*bmat;
    fval = [fval; grad];
    if iter == 2
        smat(iter,:) = (smat(iter-1,:) + grad)./2;
    else
        smat(iter,:) = (smat(iter-1,:) + del.*sum(grad))./2;
    end
    if iter >= max([JMIN,5])
        ind = (iter-4):iter;
        [ss,dss] = polintmat(h(ind),smat(ind,:),0);
        if all(abs(dss) < EPS*max(abs(ss))) || iter == JMAX
            % successful convergence
            % sort argument values and corresponding function values
            [tval,ordind] = sort(tval);
            fval  = fval(ordind,:);
            % set up partial integral values
            integfval = (tval(2) - tval(1)).*cumtrapz(fval);
            gval = interp1(tval, integfval, x, 'cubic');
            return;
        end
    end
    h(iter+1) = 0.25*h(iter);
    nx = nx*2;
end
warning(['No convergence after ',num2str(JMAX),' steps in MONGRAD']);
