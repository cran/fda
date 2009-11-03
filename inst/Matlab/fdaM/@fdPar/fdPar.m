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
%  PENMAT   ... The penalty matrix.  Ordinarily, smooth_basis will
%               compute this matrix by a call to function eval_penalty.
%               For problems of small to medium size, for which the
%               number of sampling points is roughly 500 or less, the 
%               time required for this computation is of no great 
%               consequence.  However, for large problems where 
%               repeated smoothing if required, often in determining
%               the desired level of smoothing defined by the value of
%               lambda, one will want to avoid this overhead.  If PENMAT 
%               is supplied along with the other aspects of a functional 
%               parameter object, then function smooth_basis and the many 
%               other functions having fdPar objects as arguments can
%               avoid this overhead.  Do NOT do this, however, if you
%               intend to change LFDOBJ in any way; each such change
%               requires a new penalty matrix computation.
%
%  An alternative argument list:
%  The first argument can also be a basis object.  In this case, an
%  FD object is set up with an empty coefficient matrix.  
%  For many purposes, the coefficient array is either not needed, or
%  supplied later.
%
%  Return:
%  FDPAROBJ ... A functional parameter object

%  last modified 21 July 2011

superiorto('double', 'struct', 'cell', 'char', 'inline', 'basis');

if nargin == 0
    %  case of no argument
    fdobj    = fd;
    Lfdobj   = int2Lfd(0);
    lambda   = 0;
    estimate = true;
    penmat   = [];
    
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
        if nargin < 4;  estimate = true;           end
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
    
    if ~islogical(estimate) 
        if estimate ~= 0 && estimate ~= 1
            error('ESTIMATE is not logical.');
        end
    end
    
    %  check penmat
    
    if ~isempty(penmat)
        if ~isnumeric(penmat)
            error('PENMAT is not numeric.');
        end
        penmatsize = size(penmat);
        dropind = getdropind(basisobj);
        if any(penmatsize ~= nbasis - length(dropind))
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
