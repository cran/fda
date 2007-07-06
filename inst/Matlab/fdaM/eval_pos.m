function evalarray = eval_pos(evalarg, Wfdobj, Lfdobj)
%  Evaluates a value or a derivative of a positive functional  
%  data object. 
%  A positive functional data object h  is = the form
%           h(x) = (exp Wfdobj)(x)
%  Note that the first two arguments may be interchanged.
%
%  Arguments:
%  EVALARG ... A vector of values at which all functions are to 
%              evaluated.
%  WFDOBJ  ... Functional data object.  It must define a single
%              functional data observation.
%  LFDOBJ  ... A linear differential operator object
%              applied to the functions that are evaluated.
%
%  Returns:  An array of function values corresponding to the 
%              argument values in EVALARG

%  This function is identical to EVAL_POS.

%  Last modified 20 July 2006

if nargin < 2
    error('Number of arguments is less than 2.');
end

%  check LFDOBJ and convert an integer to Lfd if needed.

if nargin < 3 
    %  set default LFDOBJ to 0
    Lfdobj = int2Lfd(0); 
end

Lfdobj = int2Lfd(Lfdobj);

%  Exchange the first two arguments if the first is an FD object
%    and the second numeric

if isnumeric(Wfdobj) && isa_fd(evalarg)
    temp    = Wfdobj;
    Wfdobj   = evalarg;
    evalarg = temp;
end

%  Check the arguments

if ~(isnumeric(evalarg))
    error('Argument EVALARG is not numeric.');
end

%  transpose EVALARG if necessary to make it a column vector

evaldim = size(evalarg);
if evaldim(1) == 1 && evaldim(2) > 1  
    evalarg = evalarg';  
end

%  check EVALARG

sizeevalarg = size(evalarg);
if sizeevalarg(1) > 1 && sizeevalarg(2) > 1
    error('Argument EVALARG is not a vector.');
end
evalarg = evalarg(:);

%  check FDOBJ

if ~isa_fd(Wfdobj)
    error('Argument FD is not a functional data object.');
end

%  Extract information about the basis

basisfd  = getbasis(Wfdobj);
rangeval = getbasisrange(basisfd);

%  determine the highest order of derivative NDERIV required

nderiv = getnderiv(Lfdobj);

%  Set up coefficient array for FD

coef  = getcoef(Wfdobj);

%  Case where EVALARG is a vector of values to be used for all curves

evalarg(evalarg < rangeval(1)-1e-10) = NaN;
evalarg(evalarg > rangeval(2)+1e-10) = NaN;
basismat = getbasismatrix(evalarg, basisfd);
fdmat    = exp(basismat*coef);

%  If a differential operator has been defined in LFDOBJ, compute
%  the weighted combination of derivatives

if nderiv > 0
    basismat  = eval_basis(evalarg, basisfd, Lfdobj);
    evalarray = fdmat.*(basismat*coef);
else
    evalarray = fdmat;
end

