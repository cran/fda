function fdParobj = fdPar(fdobj, Lfdobj, lambda, estimate, penmat)
% Sets up a functional parameter object
%  Arguments:
%  FDOBJ    ... A functional data object.  
%               The basis for this object is used to define 
%               the functional parameter, or functional 
%               parameters of FDOBJ has replications.  
%               When an initial value is required for iterative 
%               estimation of a functional parameter, the coefficients
%               will give the initial values for the iteration.
%  LFDOBJ   ... A linear differential operator value or a derivative
%               value for penalizing the roughness of the object.
%               By default, this is 0.
%  LAMBDA   ... The penalty parameter controlling the smoothness of
%               the estimated parameter.  By default this is 0.
%  ESTIMATE ... If nonzero, the parameter is estimated; if zero, the
%               parameter is held fixed at this value.
%               By default, this is 1.
%  PENMAT   ... The penalty matrix.  
%               In repeated calls to SMOOTH_BASIS, if this is
%               saved, then the penalty does not need evaluating
%               repeatedly.  Don't use, though, if LFDOBJ or LAMBDA
%               are changed in the calculation.
%
%  An alternative argument list:
%  The first argument can also be a basis object.  In this case, an
%  FD object is set up with an empty coefficient matrix.  
%  For many purposes, the coefficient array is either not needed, or
%  supplied later.
%
%  Return:
%  FDPAROBJ ... A functional parameter object

%  last modified 1 November 2007

superiorto('double', 'struct', 'cell', 'char', ...
           'inline', 'basis');

if nargin == 0
    %  case of no argument
    fdobj    = fd;
    Lfdobj   = int2Lfd(0);
    lambda   = 0;
    estimate = 1;
    penmat   = 0;
    
else
    
    if isa_basis(fdobj)
        %  if the first argument is a basis object, convert it to
        %  a default FD object with an empty coefficient matrix.
        nbasis  = getnbasis(fdobj);
        coef    = zeros(nbasis,1);
        fdobj   = fd(coef, fdobj);
    end
    
    if isa_fd(fdobj)
        basisobj = getbasis(fdobj);
        nbasis   = getnbasis(getbasis(fdobj));
        if nargin < 5;  penmat   = [];             end
        if nargin < 4;  estimate = 1;              end
        if nargin < 3;  lambda   = 0;              end
        if nargin < 2;  Lfdobj   = int2Lfd(0);     end
        
    else
        error(['First argument is neither a functional data object nor ', ...
                'a basis object.']);
    end
    
    %  check Lfdobj
    
    Lfdobj = int2Lfd(Lfdobj);
    if ~isa_Lfd(Lfdobj)
        error('LFDOBJ is not a linear differential operator object.');
    end
    
    %  check lambda
    
    if ~isnumeric(lambda)
        error('LAMBDA is not numeric.');
    end
    if lambda < 0
        error('LAMBDA is negative.');
    end
    
    %  check estimate
    
    if ~isnumeric(estimate)
        error('ESTIMATE is not numeric.');
    end
    
    %  check penmat
    
    if ~isempty(penmat)
        if ~isnumeric(penmat)
            error('PENMAT is not numeric.');
        end
        penmatsize = size(penmat);
        if any(penmatsize ~= nbasis)
            error('Dimensions of PENMAT are not correct.');
        end
    end
    
end

%  set up the fdPar object

fdParobj.fd       = fdobj;
fdParobj.Lfd      = Lfdobj;
fdParobj.lambda   = lambda;
fdParobj.estimate = estimate;
fdParobj.penmat   = penmat;

fdParobj = class(fdParobj, 'fdPar');
