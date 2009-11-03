function [evalarray, basisobj] = eval_basis(evalarg, basisobj, Lfdobj)
%  EVAL_BASIS evaluates a basis at argument values EVALARG.
%
%  LFDOBJ is a functional data object defining the order m 
%  HOMOGENEOUS linear differential operator of the form
%  Lx(t) = w_0(t) x(t) + ... + w_{m-1}(t) D^{m-1}x(t) + D^m x(t) + 
%
%  Arguments:
%  EVALARG ... A vector of values at which all functions are to 
%              evaluated.
%  BASISOBJ ... A basis object
%  LFDOBJ   ... A linear differential operator object
%              applied to the functions that are evaluated.
%
%  Note that the first two arguments may be interchanged.
%
%  Returns:  An array of function values corresponding to the evaluation
%              arguments in EVALARG

%  Last modified 9 June 2010

if nargin < 2
    error('Number of arguments is less than 2.');
end

%  set default LFDOBJ to 0

if nargin < 3     
    Lfdobj = int2Lfd(0); 
end

%  check LFDOBJ

Lfdobj = int2Lfd(Lfdobj);

%  Exchange the first two arguments if the first is a BASIS
%    and the second numeric

if isnumeric(basisobj) && isa_basis(evalarg)
    temp     = basisobj;
    basisobj = evalarg;
    evalarg  = temp;
end

%  check EVALARG

sizeevalarg = size(evalarg);
type = getbasistype(basisobj);
if ~strcmp(type, 'FEM')
    if sizeevalarg(1) > 1 && sizeevalarg(2) > 1
        error('Argument EVALARG is not a vector.');
    end
    evalarg = evalarg(:);
else
    if size(evalarg,2) ~= 2
        error('Argument EVALARG is not a matrix with two columns.');
    end
end

%  check BASISOBJ

if ~isa_basis(basisobj)
    error ('Argument BASISOBJ is not a basis object.');
end

%  determine the highest order of derivative NDERIV required

nderiv = getnderiv(Lfdobj);

%  get weight coefficient functions

wfdcell = getwfdcell(Lfdobj);
    
%  get highest order of basis matrix

if nargout == 1
    evalarray = getbasismatrix(evalarg, basisobj, nderiv);
else
    [evalarray, basisobj] = getbasismatrix(evalarg, basisobj, nderiv);
end
nbasis    = size(evalarray,2);
onerow    = ones(1,nbasis);

%  Compute the weighted combination of derivatives is 
%  evaluated here if the operator is not defined by an 
%  integer and the order of derivative is positive.
%  Only the homogeneous part of the operator, defined by
%  cell object WFDCELL is used.

if nderiv > 0 && ~isinteger(Lfdobj)
    %  In this version, only a scalar operator is allowed.
    if size(wfdcell,1) ~= 1 && size(wfdcell,2) ~= 1
        error('WFDCELL has more than one row.');
    end
    for j = 1:nderiv
        wfdParj = wfdcell{j};
        wfdj    = getfd(wfdParj);
        wcoef   = getcoef(wfdj);
        if ~all(all(wcoef == 0.0))
            wjarray   = eval_fd(evalarg, wfdj);
            evalarray = evalarray + (wjarray*onerow).* ...
                getbasismatrix(evalarg, basisobj, j-1);
        end
    end
end


