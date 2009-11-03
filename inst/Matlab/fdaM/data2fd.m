function [fdobj, df, gcv, coef, SSE] = ...
                 data2fd(argvals, y, basisobj, nderiv, lambda, fdnames)
%  DATA2FD Converts an array Y of function values plus an array
%    ARGVALS of argument values into a functional data object.
%
%  A functional data object is a sample of one or more functions, called 
%    functional data observations.
%  A functional data observation consists of one or more
%    curves, each curve corresponding to a variable.
%  For example, a functional data object can be a sample of size 35 of 
%    temperature functions, one for each of 35 Canadian weather stations.
%    In this case, each observations consists of a single temperature 
%    function.
%  Or, for example, a functional data object can be a sample of size 35
%    of temperature and precipitation functional observations.  In this case
%    each observations consists of two curves, one for the temperature 
%    and one for the precipitation variable.
%  All functional objects have a one-dimensional argument.  In the above
%    examples, this argument is time measured in months or days.
%    
%  It is now possible to call data2fd with an argument sequence 
%  that permits a penalization of the size of a derivative that
%  can be specified.  That is, this gives data2fd some of the 
%  capability of smooth.basis, except for the possibility of
%  linear differential operators other than D^m.  
%
%  Arguments for this function are as follows.  The first three are necessary
%    and the fourth is optional.
%
%  ARGVALS  ... (necessary)
%  A set of argument values.  These are common to all
%    observations, and ARGVALS will be a one-dimensional vector, with one
%    element per observation.  These values need not be increasing.
%    In the weather station example for monthly data, ARGVALS is 
%    a vector of length 12, with values 0.5, 1.5,..., 11.5.
%    Argument values falling outside of the range specified in the
%    BASIS and their corresponding values in Y will not be used, 
%    but if this happens, a warning message is displayed.
%  Argument ARGVALS is necessary, and there is no default value.
%  In the original release of data2fd, it was possible to input 
%  arguments as a matrix, permitting different argument values and
%  different numbers of arguments for curves.  This option has been
%  discontinued.
%    
%  Y ... (necessary)
%  The array Y stores curve values used to create functional data object FDOBJ.
%  Y can have one, two, or three dimensions according to whether whether
%    the sample size, the number of variables in each observation.  Its 
%    dimensions are:
%     1.  argument values  ...  size = no. argument values in ARGVAL
%     2.  replications     ...  size = sample size
%     3.  variables        ...  size = no. variables per observation
%  If Y is a one-way array, either as a vector or a matrix with one column,
%     it's single non-trivial dimension = no. argument values.  If Y
%     is two-dimensional, each observation is assumed to have one variable.
%     If Y is three-dimensional, each observation is assumed to have
%     multiple variables.  Note:  a single multivariate observation must
%     be an array Y with three dimensions, the middle of which is of length 1.
%  Example:  For monthly temperature data for 35 weather stations, 
%     Y will be 12 by 35.  For both temperature and precipitation observations, 
%     Y will be 12 by 35 by 2.  For temperature/precipitation data at Montreal 
%     only, Y will be 12 by 1 by 2.
%  This argument is necessary, and there is no default value.
%
%  BASISOBJ  ...  (necessary)
%    A functional data basis object created by function CREATE_BASIS_FD 
%    or one of its specialized version, such as CREATE_BSPLINE_BASIS or
%    CREATE_FOURIER_BASIS.  The functional data basis object specifies
%    a basis type (eg. 'fourier' or 'bspline'), a range or argument values,
%    the number of basis functions, and fixed parameters determining these
%    basis functions (eg. period for 'fourier' bases or knots for 'bspline'
%    bases.  
%    In most applications, BASIS will be supplied.  If BASIS is supplied, 
%    the next three arguments are ignored.  
%    If BASIS is an essential argument, and there no default value.  But
%    see function MAKE_BASIS for a simplified technique for defining this
%    basis.  For example, function call 
%         MAKE_BASIS([0,12], 7, 1) 
%    could be used for the monthly temperature/precipitation data to define
%    a 'fourier' basis over an interval of 12 months containing 7 basis
%    functions (the 3rd argument species the basis to be periodic.)
%    This argument is necessary, and there is no default value.
%    
%    The following arguments are optional:
%
%  NDERIV   ... A non-negative integer specifying the order of derivative
%   whose size is to be controlled by the roughness penalty
%       LAMBDA \int [D^NDERIV x(t)]^2 dt
%   NDERIV may also be the string 'h', in which the harmonic acceleration
%   operator is set up with period equal to the range specified in
%   BASISOBJ.
%   The default value is 4.
%
%  LAMBDA   ... A nonegative real number specifying the weight placed on 
%    the size of the derivative.
%   The default value is 0.0.
%
%  FDNAMES  ... A cell array object of length 3 with each cell containing a 
%    single string specifying a name for a dimension of Y.
%       1. argument domain, such as the string 'Time'
%       2. replications or cases
%       3. variables
%    For example, for the daily temperature data, 
%    fdnames{1} = 'Day';
%    fdnames{2} = 'Station';
%    fdnames{3} = 'Temperature (deg C)';
%    By default, the string 'time', 'reps' and 'values' are used.
%
%    These optional arguments can be supplied in two ways:
%    1. If only a fourth argument is supplied, and it is a cell object,
%       then it is taken to be FDNAMES.  In this case, no roughness
%       penalty will be used.  This is the argument sequence in the 
%       original release of this function.  On the other hand, if it is
%       a non-negative integer, then it is taken to be  NDERIV and 
%       LAMBDA and FDNAMES are set to their default values
%    2. Up to three arguments are supplied, and are assumed to be in the
%       order that follows: first NDERIV, second LAMBDA, and third FDNAMES.
%       If any argument is to be set to its default value, but a following
%       argument is required, it's position should be filled by the empty object []
%
%  DATA2FD Returns these objects: 
%
%  FDOBJ ... A functional data object containing the curves that fit the
%    data in a least squares sense.
%  DF    ... A real number specifying the equivalent of degrees of freedom
%    for the fitted curves.
%  GCV   ... N values of the GCV criterion associated with the fit, one
%    value per replication.  SUM(GCV) can be used as a selector of the
%    value of LAMBDA by searching for its minimum.
%  COEF: ... The array of coefficients for the expansions of the fitted
%    curves.  The first dimension corresponds to the basis functions, the
%    second dimension to the replications, and the third dimension if 
%    required to the variables for a multivariate functional data object.
%    That is, it is an NBASIS by N (by NVAR) matrix or array.
%  SSE   ... N values of the error sum of squares.
%
%  DATA2FD is intended for more casual smoothing not requiring a great deal  
%    sophistication in defining the functional data object.  It uses
%    function SMOOTH_BASIS to compute the functional data object. Indeed,
%    in the simplest and most common situation, DATA2FD consists of 
%    However, for more advanced applications requiring more smoothing 
%    control than is possible by setting the number of basis functions in
%    BASIS, function SMOOTH_BASIS should be used.  Or, alternatively,
%    DATA2FD may first be used with a generous number of basis functions,
%    followed by smoothing using function SMOOTH_FD.

%  Tests have now been installed to detect that Y and ARGVALS have been
%  supplied in reverse order, so that users can employ the same order as
%  that used in the other smoothing functions, namely ARGVALS before Y.

%  Last modified:  30 April 2009

 %  check the BASIS argument
 
 if nargin < 3
    error('Argument BASISOBJ is not supplied.');
 end
 if ~isa_basis(basisobj)
    error('Argument BASISOBJ is not a basis object.');
 end
 
 %  check that the sizes of the first dimensions match
 
 if size(y,1) ~= size(argvals,1)
     error('First dimensions of Y and ARGVALS do not match.');
 end
     
 %  check whether ARGVALS and Y should be swapped
 
 [y, argvals] = argvalsy_swap(argvals, y, basisobj);
 
 %  check the dimensions of Y
 
 yd   = size(y);
 
 if yd(1) == 1
     if length(yd) == 2
        y  = y';
        yd = size(y);
     else
        error(['Y is an array and length of ', ...
               'first dimension = 1.']);
     end
 end
 ydim = length(yd);
 if ydim > 3
    error('Too many dimensions for argument Y.');
 end
 %  set number of sampled values
 n    = yd(1);
 if n == 1
    error('Only one value in ARGVALS not allowed.');
 end
 
 %  check argument value array

 argd  = size(argvals);
 if length(argd) > 2
    error('Argument ARGVALS has too many dimensions.');
 end
 if min(diff(argvals)) < 0
     error('Values in ARGVALS are not in increasing order.');
 end
 %  Test if arguments are outside of the range in the basis.
 rangeval = getbasisrange(basisobj);
 if min(argvals) < rangeval(1) || max(argvals) > rangeval(2)
    error('Some arguments values are outside of the range in BASISOBJ.');
 end
 
 %  set default values for optional parameters
 
 defaultnames = cell(3,1);
 defaultnames{1} = 'argument values';
 defaultnames{2} = 'replications';
 defaultnames{3} = 'variables';
 
 if nargin > 4 || nargin < 4
     if nargin < 6,  fdnames = defaultnames;  end
     if nargin < 5,  lambda  = 0.0;           end
     if nargin < 4,  nderiv  = 4;             end
 else
     if iscell(nderiv)  
         fdnames = nderiv; 
         nderiv  = 4;
         lambda  = 0.0;
     end
 end
 
 %  check values of optional arguments
 
 if ~isnumeric(nderiv) && ~ischar(nderiv)
     error('NDERIV is neither a number nor a character.');
 end
 if isnumeric(nderiv)
     if nderiv < 0
         error('NDERIV is negative.');
     end
     if floor(nderiv) ~= nderiv
         error('NDERIV is not an integer.');
     end
 else
     if strcmp(nderiv, 'h')
         %  ----  set up the harmonic acceleration operator  -------
         Lbasis  = create_constant_basis(rangeval);
         Lcoef   = [0,(2*pi/(rangeval(2)-rangeval(1)))^2,0];
         wfd     = fd(Lcoef,Lbasis);
         wfdcell = fd2cell(wfd);
         nderiv  = Lfd(3, wfdcell);
         if nargin < 6,  fdnames = defaultnames;  end
         if nargin < 5,  lambda  = 0.0;           end
     else
         error('NDERIV is not an allowed character.');
     end
 end
 
 if ~isnumeric(lambda)
     error('LAMBDA is not a number.');
 end
 if lambda < 0.0
     error('LAMBDA is negative.');
 end
 
 if ~iscell(fdnames)
     error('FDNAMES is not a cell object');
 end
 if length(fdnames) < 3
     fdnames{3} = 'variables';
 end
 if length(fdnames) < 2
     fdnames{2} = 'replications';
 end
 
 %  do the smoothing
 
 fdParobj = fdPar(basisobj, nderiv, lambda);
 
 [fdobj, df, gcv, coef, SSE] = ...
                    smooth_basis(argvals, y, fdParobj, [], fdnames);

    
    




