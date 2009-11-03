function outvec = eval_mon(evalarg, Wfdobj, Lfdobj)
%  Evaluates a single monotone functional object, 
%    or one of its derivatives. 
%  The monotone functional data object h  is = the form
%           h(x) = (D^{-1} exp Wfdobj)(x)
%  where  D^{-1} means taking the indefinite integral.
%  Note that the linear differential operator object LFDOBJ
%  MUST be an integer in the range 0 to 3.
%  Note that the first two arguments may be interchanged.
%
%  Arguments:
%  EVALARG ... A vector of values at which all functions are to 
%              evaluated.
%  WFDOBJ   ... Functional data object.  It must define a single
%              functional data observation. 
%  LFDOBJ  ... A linear differential operator object
%              applied to the functions that are evaluated.
%
%  Returns:  An array of function values corresponding to the evaluation
%              arguments in EVALARG

%  This function is identical to EVAL_MONFD.

%  Last Modified 14 August 2006

if nargin < 2
    error('Number of arguments is less than 2.');
end

%  check LFDOBJ and convert an integer to Lfd if needed.

if nargin < 3 
    %  set default LFDOBJ to 0
    Lfdobj = int2Lfd(0); 
else
    %  check LFDOBJ
    if isnumeric(Lfdobj)
        % if integer, check and convert to Lfd
        nderiv = Lfdobj;
        if nderiv ~= round(Lfdobj)
            error('LFDOBJ numeric but not an integer.');
        end
        if nderiv < 0
            error('LFDOBJ an integer but negative.');
        end
        Lfdobj = int2Lfd(nderiv);
    else
        if ~isa_Lfd(Lfdobj)
            error (['Argument LFDOBJ is neither a functional data object', ...
                    ' nor an integer.']);
        end
        if ~isinteger(Lfdobj)
            error('LFD is not a D^m operator.');
        end        
    end
end

%  Exchange the first two arguments if the first is an FD object
%    and the second numeric

if isnumeric(Wfdobj) && isa_fd(evalarg)
    temp    = Wfdobj;
    Wfdobj   = evalarg;
    evalarg = temp;
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

%  check LFDOBJ

if ~(isa_Lfd(Lfdobj))
    error('LFD is not linear differential operator object.');
end

%  Check FDOBJ

coef  = getcoef(Wfdobj);
coefd = size(coef);
ndim  = length(coefd);
if ndim > 1 && coefd(2) ~= 1  
    error('FDOBJ is not a single function');
end

nderiv = getnderiv(Lfdobj);

if nderiv == 0
    hval = monfn(evalarg, Wfdobj);
    outvec = hval;
    return;
end

if nderiv == 1
    Dhval = exp(eval_fd(Wfdobj, evalarg));
    outvec = Dhval;
    return;
end

if nderiv == 2
    basisobj = getbasis(Wfdobj);
    Dwmat    = getbasismatrix(evalarg, basisobj, 1);
    D2hval   = (Dwmat * coef) .* exp(eval_fd(Wfdobj, evalarg));
    outvec   = D2hval;
    return;
end

if nderiv == 3
    basisobj = getbasis(Wfdobj);
    Dwmat    = getbasismatrix(evalarg, basisobj, 1);
    D2wmat   = getbasismatrix(evalarg, basisobj, 2);
    D3hval   = ((D2wmat .* coef) + (Dwmat * coef).^2) * ...
        exp(eval_fd(Wfdobj, evalarg));
    outvec = D3hval;
    return;
end

if nderiv > 3 
    error ('Derivatives higher than 3 not implemented.');
end


    
