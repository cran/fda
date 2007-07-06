function df = lambda2df(argvals, basisobj, wtvec, Lfdobj, lambda)
%  LAMBDA2DF computes degrees of freedom for a regularized basis smooth
%    by computing the trace of the hat matrix.

%  Arguments for this function:
%
%  ARGVALS  ... A set of argument values, set by default to equally spaced on
%               the unit interval (0,1).
%  BASISOBJ  ... A basis.fd object created by function create_basis.fd.
%  WTVEC    ... A vector of N weights, set to one by default, that can
%               be used to differentially weight observations = the
%               smoothing phase
%  LFDOBJ   ... The order of derivative or a linear differential
%               operator to be penalized.
%               By default Lfdobj is set to 2.
%  LAMBDA   ... The smoothing parameter determining the weight to be
%               placed on the size of the derivative = smoothing.  This
%               is 0 by default.
%  Returns:
%  DF    ...  a degrees of freedom measure

%  last modified 2 December 2006

%  set default values

n = length(argvals);
if nargin < 5
    lambda = 0;
end

if nargin < 4
    Lfdobj = int2Lfd(2);
end;

if nargin < 3
    wtvec = ones(n,1);
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

nbasis = getnbasis(basisobj);
if lambda == 0
    df = nbasis;
    return;
end

sizew = size(wtvec);
if (length(sizew) > 1 && sizew(1) > 1 && sizew(2) > 1) || ...
        length(sizew) > 2
    error ('WTVEC must be a vector.');
end
if length(sizew) == 2 && sizew(1) == 1
    wtvec = wtvec';
end
if length(wtvec) ~= n
    error('WTVEC of wrong length');
end
if min(wtvec) <= 0
    error('All values of WTVEC must be positive.');
end
if lambda < 0
    warning ('Wid:negative', ...
             'Value of LAMBDA was negative, and 0 used instead.');
    lambda = 0;
end
basismat  = getbasismatrix(argvals, basisobj);
basisw = basismat .* (wtvec * ones(1,nbasis));
Bmat    = basisw' * basismat;
penmat  = eval_penalty(basisobj, Lfdobj);
Bnorm   = sqrt(sum(sum(Bmat.^2)));
pennorm = sqrt(sum(sum(penmat.^2)));
condno  = pennorm/Bnorm;
if lambda*condno > 1e12
    lambda = 1e12/condno;
%     disp(['lambda reduced to ',num2str(lambda),' to prevent overflow']);
end
Cmat   = Bmat + lambda .* penmat;
if is_diag(Cmat)
    Cmatinv = diag(1/diag(Cmat));
else
    Cmatinv = inv(Cmat);
end
df = sum(diag(Cmatinv * Bmat));
