function [hval, basiscell, ifval] = monfn(x, Wfd, basiscell, EPS, JMIN, JMAX)
%MONFN evaluates a monotone function  h(x) = (D^{-1} exp Wfd)(x)
%  where  D^{-1} means taking the indefinite integral.
%  The interval over which the integration takes places is defined in
%       the basis object in WFD.
%  A major overhead in this function and MONGRAD is the evaluation
%  of the basis functions at the dyadic argument sequence required for
%  the trapezoidal rule integration.  This is reduced by storing these
%  values in the cells in cell array BASISCELL, each cell corresponding
%  to an iteration of the trapezoidal rule, with a maximum of 15.  

%  Arguments:
%  X         ... argument values at which function and derivatives are evaluated
%  WFD       ... a functional data object
%  BASISCELL ...  A cell array object of length 15
%                 containing basis function values.
%  EPS       ...  Relative error needed for convergence. Default 1e-5
%  JMIN      ...  Minimum number of step halving steps.  Default 11
%  JMAX      ...  Maximum number of step halving steps.  Default 15
%
%  Returns:
%  HVAL      ... value of h at input argument array X in first column.
%  BASISCELL ... An updated cell object for BASISCELL.
%  TVAL      ... Arguments used for trapezoidal approximation to integral

%  Last modified 31 October 2007

if nargin < 2,  error('There are less than two arguments');  end

%  set some constants

if nargin < 6,  EPS  = 1e-4;  end
if nargin < 5,  JMIN = 11;    end
if nargin < 4,  JMAX = 15;    end
if JMIN   < 5,  JMIN =  5;    end

%  get coefficient matrix and check it

coef  = getcoef(Wfd);
coefd = size(coef);
ndim  = length(coefd);
if ndim > 2  
    error('WFD is not a univariate function');
end
ncurve = coefd(2);

%  get the basis

basis  = getbasis(Wfd);
rng    = getbasisrange(basis);
width  = rng(2) - rng(1);

%  return linear values if all coefficients 0

if all(coef == 0)
    tval = rng;
    hval = interp1(tval,[0,width],x);
    return;
end

%  set up first iteration

JMAXP = JMAX + 1;
h     = ones(JMAXP,1);
h(2)  = 0.25;
%  matrix SMAT contains the history of discrete approximations to the
%    integral
smath = zeros(JMAXP,ncurve);
%  array TVAL contains the argument values used = the approximation
%  array FVAL contains the integral values at these argument values,
%     rows corresponding to argument values
%  the first iteration uses just the endpoints
iter    = 1;
xiter   = rng';
tval = xiter;
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
fval = fx;
smath(1,:) = (width/2).*sum(fx);
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
    fval = [fval; fx];
    smath(iter,:) = ( smath(iter-1,:) + del.*sum(fx) )./2;
    if iter >= max([JMIN,5])
        ind = (iter-4):iter;
        [ss,dss] = polintmat(h(ind),smath(ind,:),0);
        if all(abs(dss) < EPS*max(abs(ss)))
            % successful convergence
            % sort argument values and corresponding function values
            [tval,ordind] = sort(tval);
            % set up partial integral values
            ifval = (tval(2) - tval(1)).*cumtrapz(fval(ordind,:));
            hval  = interp1(tval, ifval, x, 'cubic');
            return;
        end
    end
    h(iter+1)       = 0.25*h(iter);
    nx = nx*2;
end
warning('Wid:convergence',...
    ['No convergence after ',num2str(JMAX),' steps in MONFN']);
[tval,ordind] = sort(tval);
% set up partial integral values
ifval = (tval(2) - tval(1)).*cumtrapz(fval(ordind,:));
hval  = interp1(tval, ifval, x, 'cubic');
