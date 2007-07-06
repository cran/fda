function  smthfd = smooth_fd(fdobj, fdParobj, rebase)
%SMOOTH_FD smooths a functional data object.
%
%  Arguments for this function:
%
%  FDOBJ    ... A functional data object.
%  FDPAROBJ ... A functional parameter object.
%
%  If rebase=1 and the basis type is 'polyg' then the basis
%    is changed to a cubic bspline  basis and before smoothing
%
%  Returns a functional data object containing a smoothed version
%    of the input functional data object
%
%  Last modified:  2 December 2006

%
% Rebase to default B spline basis if rebase is T and basistype is
%    polygonal.  Then test to see if any smoothing is actually required.
%

if nargin < 2
    error('There are less than two arguments.');
end

%  set default values

if nargin < 3, rebase = 1;  end

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

%  check LFD

Lfdobj = getLfd(fdParobj);
Lfdobj = int2Lfd(Lfdobj);

%  check FDOBJ

if ~isa_fd(fdobj)
    error('FDOBJ is not a functional data object.');
end

%  set up basis

basisobj = getbasis(fdobj);
type  = getbasistype(basisobj);
if rebase == 1 && strcmp(type,'polyg')
    fdnames = getnames(fdobj);
    params = getbasispar(basisobj);
    fdobj = data2fd(getcoef(fdobj), params, fdnames);
    basisobj = getbasis(fdobj);
end

%  check lambda

lambda = getlambda(fdParobj);
if lambda <= 0
    warning('Wid:positive', ...
        'LAMBDA was not positive. No smoothing carried out.');
    smthfd = fdobj;
    return;
end
%
%  Main smoothing step
%
coef  = getcoef(fdobj);
coefd = size(coef);
ndim  = length(coefd);
if ndim == 3
     nvar = coefd(3);
else
     nvar = 1;
end
Bmat  = inprod_basis(basisobj);
%
%  set up coefficient matrix for normal equations
%
penmat = eval_penalty(basisobj, Lfdobj);
Cmat   = Bmat + lambda .* penmat;
%
%  solve normal equations for each observation
%
if(ndim < 3)
     Dmat = inprod(basisobj, fdobj);
     coef = symsolve(Cmat, Dmat);
else
    for ivar = 1:nvar
        Dmat = inprod(basisobj, fdobj(:,ivar));
        coef(:,:,ivar) = symsolve(Cmat, Dmat);
    end
end
%
%  replace coefficient matrix = fdobj, leaving other properties alone
%
fdnames = getnames(fdobj);
smthfd = fd(coef, basisobj, fdnames);


