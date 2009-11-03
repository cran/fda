function [Sigma, basisobj] = eval_fdVar_fd(evalarg, fdobj, Lfdobj)
%  EVAL_FD evaluates a functional data observation at argument
%  values EVALARG for a fdVariance type basis object
%
%  LFDOBJ is a functional data object defining the order m
%  HOMOGENEOUS linear differential operator of the form
%  Lx(t) = w_0(t) x(t) + ... + w_{m-1}(t) D^{m-1}x(t) +
%          \exp[w_m(t)] D^m x(t) + ...
%
%  Arguments:
%  EVALARG ... A vector of values at which all functions are to
%              evaluated.
%  FDOBJ   ... Functional data object
%  LFDOBJ  ... A linear differential operator object
%              applied to the functions that are evaluated.
%
%  Note that the first two arguments may be interchanged.
%
%  Returns:  An cell array of function values corresponding
%              to the evaluation arguments in EVALARG

%  Last modified 27 September 2011

%  Check arguments

if nargin < 2
    error('Number of arguments is less than 2.');
end

%  Set default arguments

if nargin == 3
    if Lfdobj ~= 0 || ~isempty(Lfdobj)
        error('Lfdobj is not either 0 or empty for fdVariance object.');
    end
end

%  Exchange the first two arguments if the first is an FD object
%    and the second numeric

if isnumeric(fdobj) && isa_fd(evalarg)
    temp    = fdobj;
    fdobj   = evalarg;
    evalarg = temp;
end

%  check EVALARG

if isnumeric(evalarg)
sizeevalarg = size(evalarg);
if sizeevalarg(1) > 1 && sizeevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
evalarg = evalarg(:);
elseif isstruct(evalarg)
    if ~isfield(evalarg, 'pts')
        error('Argument EVALARG does not contain a field pts.');
    end
    if ~isfield(evalarg, 'z')
        error('Argument EVALARG does not contain a field z.');
    end
else
    error('Argument evalarg is neither numeric nor a struct object.');
end

%  check FDOBJ

if ~isa_fd(fdobj)
    error('Argument FD is not a functional data object.');
end

%  Extract information about the basis

basisobj = getbasis(fdobj);
type     = getbasistype(basisobj);

if ~strcmp(type, 'fdVariance')
    error('BASISOBJ is not of type fdVariance');
else
    cvec    = getcoef(fdobj);
    if nargout == 1
        RstCell = eval_basis(evalarg, basisobj);
    else
        [RstCell, basisobj] = eval_basis(evalarg, basisobj);
    end
    Sigma   = fdVar_Sigma(cvec, RstCell);      
end
