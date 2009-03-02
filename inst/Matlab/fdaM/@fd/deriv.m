function derivfd = deriv(fdobj, Lfdobj)
%  DERIV  Applies linear differential operator object LFDOBJ 
%  to functional data object FDOBJ.
%  LFDOBJ is either a positive integer or a
%    a linear differential operator.

%  last modified 3 March 2009

%  check the linear differential operator object LFDOBJ

if nargin < 2
    Lfdobj = int2Lfd(1);
else
    Lfdobj = int2Lfd(Lfdobj);
end

%  get basis information

basisobj = getbasis(fdobj);
nbasis   = getnbasis(basisobj);
rangeval = getbasisrange(basisobj);

%  evaluate FDOBJ for a fine mesh of argument values

nfine    = max([201, 10*nbasis+1]);
evalarg  = linspace(rangeval(1), rangeval(2), nfine)';
Lfdmat   = eval_fd(evalarg, fdobj, Lfdobj);

%  coefficient matrix for derivative functional data object

Lfdcoef  = project_basis(Lfdmat, evalarg, basisobj);

%  set up the derivative object

Dfdnames = getnames(fdobj);
%  Name and labels for variables
if iscell(Dfdnames{3})
    Dfdnames{3}{1} = ['L-',Dfdnames{3}{1}];
else
    if ischar(Dfdnames{3}) && size(Dfdnames{3},1) == 1
        Dfdnames{3} = ['L-',Dfdnames{3}];
    else
        Dfdnames{3} = 'L-function';
    end
end

derivfd = fd(Lfdcoef, basisobj, Dfdnames);


