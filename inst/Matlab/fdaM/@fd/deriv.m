function derivfd = deriv(fdobj, Lfdobj)
%  DERIV  Applies linear differential operator object LFDOBJ 
%  to functional data object FDOBJ.
%  LFDOBJ is either a positive integer or a
%    a linear differential operator.

%  last modified 28 January 2003

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

evalarg  = linspace(rangeval(1), rangeval(2), 10*nbasis+1)';
Lfdmat   = eval_fd(evalarg, fdobj, Lfdobj);

%  coefficient matrix for derivative functional data object

Lfdcoef  = project_basis(Lfdmat, evalarg, basisobj);

%  set up the derivative object

Dfdnames    = getnames(fdobj);
Dfdnames{3} = ['D',Dfdnames{3}];

derivfd = fd(Lfdcoef, basisobj, Dfdnames);


