function [hessmat, tval] = monhess(x, Wfd, basiscell)
%MONHESS evaluates the second derivative of monotone fn. wrt coefficients
%  The function is of the form h(x) = (D^{-1} exp Wfd)(x)
%  where  D^{-1} means taking the indefinite integral.
%  The interval over which the integration takes places is defined in
%       the basis object = WFD.
%  The derivatives with respect to the coefficients in WFD up to order
%       NDERIV are also computed, max(NDERIV) = 2.
%  Arguments:
%  X      ... argument values at which function and derivatives are evaluated
%             x(1) must be at lower limit, and x(n) at upper limit.
%  WFD    ... a functional data object
%  Returns:
%  D2H  ... values of D2 h wrt c
%  TVAL ... Arguments used for trapezoidal approximation to integral

%  Last modified 26 February 2012

%  set some constants

EPS  = 1e-4;
JMIN =  7;
JMAX = 11;

%  get coefficient matrix and check it

coef  = getcoef(Wfd);
coefd = size(squeeze(coef));
ndim  = length(coefd);
if ndim > 1 && coefd(2) ~= 1
    error('WFD is not a single function');
end

%  get the basis

basisobj = getbasis(Wfd);
rng      = getbasisrange(basisobj);
nbasis   = getnbasis(basisobj);
params   = getbasispar(basisobj);
norder   = nbasis - length(params);

%  Compute the number of active pairs.
%  Note that the first basis is not active, but we still
%    need its space in the array.
nbaspr = nbasis*norder - norder*(norder-1)/2;

%  set up first iteration

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h     = ones(JMAXP,1);
h(2)  = 0.25;
%  matrix SMAT contains the history of discrete approximations to the
%    integral
smat = zeros(JMAXP,nbaspr);
%  array TVAL contains the argument values used = the approximation
%  array FVAL contains the integral values at these argument values,
%     rows corresponding to argument values
%  the first iteration uses just the endpoints
iter    = 1;
xiter   = rng';
tval    = xiter;
if nargin == 3
    if isempty(basiscell{iter})
        bmat = getbasismatrix(xiter, basisobj);
        basiscell{iter} = bmat;
    else
        bmat = basiscell{iter};
    end
else
    bmat = getbasismatrix(xiter, basisobj);
end
fx   = exp(bmat*coef);
D2fx = zeros(2,nbaspr);
m = 0;
for ib=1:nbasis
    for jb=ib:min(nbasis,ib+norder-1)
        m=m+1;
        D2fx(:,m) = fx.*bmat(:,ib).*bmat(:,jb);
    end
end
D2fval = D2fx;
smat(1,:) = width*sum(D2fx)./2;
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
            bmat = getbasismatrix(xiter, basisobj);
            basiscell{iter} = bmat;
        else
            bmat = basiscell{iter};
        end
    else
        bmat = getbasismatrix(xiter, basisobj);
    end
    fx   = exp(bmat*coef);
    D2fx = zeros(length(xiter),nbaspr);
    m = 0;
    for ib=1:nbasis
        for jb=ib:min(nbasis,ib+norder-1)
            m = m + 1;
            D2fx(:,m) = fx.*bmat(:,ib).*bmat(:,jb);
        end
    end
    D2fval = [D2fval; D2fx];
    if iter == 2
        smat(iter,:) = (smat(iter-1,:) + D2fx)./2;
    else
        smat(iter,:) = (smat(iter-1,:) + del.*sum(D2fx))./2;
    end
    if iter >= max([JMIN,5])
        ind = (iter-4):iter;
        [D2ss, D2dss] = polintmat(h(ind),smat(ind,:),0);
        if all(abs(D2dss) < EPS*max(abs(D2ss))) || iter == JMAX
            % successful convergence
            % sort argument values and corresponding function values
            [tval,ordind] = sort(tval);
            % set up partial integral values
            del     = tval(2) - tval(1);
            D2ifval = del.*cumtrapz(D2fval(ordind,:));
            hessmat = safeinterp(tval, D2ifval, x);
            return;
        end
    end
    h(iter+1) = 0.25*h(iter);
    nx = nx*2;
end
%err
