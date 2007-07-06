%  ---------------------------------------------------------------

function ss = funcint(fhdl, cvec, basisobj, nderiv1, nderiv2, ...
                     JMAX, EPS)
%  integrates a function of event times and basis values
%  Arguments:
%  FHDL     ... a function handle for the function
%  CVEC     ... a coefficient vector
%  BASISOBJ ... a functional data object for log intensity function w(t)
%  NDERIV1  ... a derivative order 
%  NDERIV2  ... another derivative order

%  Last modified 20 July 2006

if ~strcmp(class(basisobj),'basis')
    error('Fourth argument must be a basis function object.');
end

%  set default arguments

if nargin < 7, EPS  = 1e-7; end
if nargin < 6, JMAX = 9;    end
if nargin < 5, nderiv2 = 0; end
if nargin < 4, nderiv1 = 0; end

rng = getbasisrange(basisobj);

%  set up first iteration

width = rng(2) - rng(1);
JMAXP = JMAX + 1;
h     = ones(JMAXP,1);
h(2) = 0.25;
%  matrix SMAT contains the history of discrete approximations to the integral
%  the first iteration uses just the endpoints
x  = rng';
if nargin == 3
    fval = feval(fhdl, x, cvec, basisobj);
elseif nargin == 4
    fval = feval(fhdl, x, cvec, basisobj, nderiv1);
else
    fval = feval(fhdl, x, cvec, basisobj, nderiv1, nderiv2);
end
scell{1} = width.*fval./2;
tnm = 0.5;

%  now iterate to convergence

for iter = 2:JMAX
    tnm  = tnm*2;
    del  = width/tnm;
    if iter == 2
        x = (rng(1) + rng(2))/2;
    else
        x = (rng(1)+del/2 : del : rng(2))';
    end
    if nargin == 3
        fval = feval(fhdl, x, cvec, basisobj);
    elseif nargin == 4
        fval = feval(fhdl, x, cvec, basisobj, nderiv1);
    else
        fval = feval(fhdl, x, cvec, basisobj, nderiv1, nderiv2);
    end
    scell{iter} = (scell{iter-1} + width.*fval./tnm)./2;
    if iter >= 5
        ind = (iter-4):iter;
        temp = { scell{iter-4}, scell{iter-3}, scell{iter-2}, ...
                 scell{iter-1}, scell{iter} }; 
        [ss, dss] = polintarray(h(ind),temp,0);
        if ~any(abs(dss) > EPS.*max(max(abs(ss))))
            %  successful convergence
            ss = squeeze(ss);
            return;
        end
    end
    scell{iter+1} = scell{iter};
    h(iter+1) = 0.25*h(iter);
end

%  ---------------------------------------------------------------

function [y,dy] = polintarray(xa, ya, x)
%  YA is an array with up to 4 dimensions
%     with 1st dim the same length same as the vector XA
n       = length(xa);
difx    = xa - x;
absxmxa = abs(difx);
tmp = 1:n;
ns = min(tmp(absxmxa == min(absxmxa)));
cs = ya;
ds = ya;
y  = ya{ns};
ns = ns - 1;
for m = 1:(n-1)
    for i = 1:(n-m)
        ho      = difx(i);
        hp      = difx(i+m);
        w       = (cs{i+1} - ds{i})./(ho - hp);
        ds{i} = hp.*w;
        cs{i} = ho.*w;
    end
    if 2*ns < n-m
        dy = cs{ns+1};
    else
        dy = ds{ns};
        ns = ns - 1;
    end
    y = y + dy;
end

