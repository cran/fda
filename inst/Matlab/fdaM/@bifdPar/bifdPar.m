function bifdParobj = bifdPar(bifdobj, estimate, lambdas, lambdat, Lfdobjs, Lfdobjt)
% Sets up a functional parameter object
%  Arguments:
%  BIFDOBJ  ... A bivariate functional data object.  The basis for this object 
%               is used to define the bivariate functional parameter.
%               When an initial value is required for iterative 
%               estimation of a bivariate functional parameter, the coefficients
%               will give the initial values for the iteration.
%  ESTIMATE ... If nonzero, the parameter is estimated; if zero, the
%                parameter is held fixed at this value.
%                By default, this is 1.
%  LAMBDAS  ... The penalty parameter controlling the smoothness of
%               the estimated parameter with respect to the first argument s.  
%               By default this is 0.
%  LAMBDAT  ... The penalty parameter controlling the smoothness of
%               the estimated parameter with respect to the second argument t.  
%               By default this is 0.
%  LFDOBJS  ... A linear differential operator value or a derivative
%               value for penalizing the roughness of the object
%               with respect to the first argument s.
%               By default, this is 2.
%  LFDOBJT  ... A linear differential operator value or a derivative
%               value for penalizing the roughness of the object
%               with respect to the second argument t.
%               By default, this is 2.

%  last modified 20 July 2006

superiorto('double', 'sparse', 'struct', 'cell', 'char', ...
    'inline', 'basis');

%  check BIFDOBJ

if ~isa_bifd(bifdobj)
    error('BIFDOBJ is not a bivariate functional data object.');
end

%  set some default argument values

if nargin < 6;  Lfdobjt   = int2Lfd(2);  end
if nargin < 5;  Lfdobjs   = int2Lfd(2);  end
if nargin < 4;  lambdat   = 0;           end
if nargin < 3;  lambdas   = 0;           end
if nargin < 2;  estimate  = 1;           end

%  check the linear differential operators

Lfdobjs = int2Lfd(Lfdobjs);
Lfdobjt = int2Lfd(Lfdobjt);

if ~isa_Lfd(Lfdobjs)
    error('LFDOBJS is not a linear differential operator object.');
end
if ~isa_Lfd(Lfdobjt)
    error('LFDOBJT is not a linear differential operator object.');
end

%  check the roughness penalty parameters

if ~isnumeric(lambdas)
    error('LAMBDAS is not numeric.');
end
if lambdas < 0
    warning('LAMBDAS is negative, and is set to zero.');
    lambdas = 0;
end
if ~isnumeric(lambdat)
    error('LAMBDAT is not numeric.');
end
if lambdat < 0
    warning('LAMBDAT is negative, and is set to zero.');
    lambdat = 0;
end

if ~isnumeric(estimate)
    error('ESTIMATE is not numeric.');
end

%  set up the bifdPar object

bifdParobj.bifd     = bifdobj;
bifdParobj.estimate = estimate;
bifdParobj.lambdas  = lambdas;
bifdParobj.Lfds     = Lfdobjs;
bifdParobj.lambdat  = lambdat;
bifdParobj.Lfdt     = Lfdobjt;

bifdParobj = class(bifdParobj, 'bifdPar');
