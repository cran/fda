function fdobj = data2fd(y, argvals, basisobj, fdnames)
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
%  The data values used by DATA2FD to make the functional data object FDOBJ
%    are in array Y.  If each observation consists of a single curve, 
%    Y is two-dimensional If each observations has multiple variables, 
%    Y will be three-dimensional.  See below for further details.
%
%  DATA2FD assumes that each observation is evaluated at a set of argument 
%    specified in ARGVALS. These argument values may be common to all curves,
%    in which case ARGVALS is a one-dimensional vector.  If argument values
%    vary from observation to observation, ARGVALS will be two dimensional.
%    See below for further details.
%
%  Arguments for this function are as follows.  The first three are necessary
%    and the fourth is optional.
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
%  The values in Y may be missing, indicated by NaN.  The presence of missing
%     values will slow down the computation somewhat since each observation
%     must then be processed separately.
%  Example:  For monthly temperature data for 35 weather stations, 
%     Y will be 12 by 35.  For both temperature and precipitation observations, 
%     Y will be 12 by 35 by 2.  For temperature/precipitation data at Montreal 
%     only, Y will be 12 by 1 by 2.
%  This argument is necessary, and there is no default value.
%
%  ARGVALS  ... (necessary)
%  A set of argument values.  In most situations, these will be common to all
%    observations, and ARGVALS will be a one-dimensional vector, with one
%    element per observation.  These values need not be increasing.
%    In the weather station example for monthly data, ARGVALS is 
%    a vector of length 12, with values 0.5, 1.5,..., 11.5.
%    However, ARGVALS may also be a matrix if the argument values vary from
%    observation to observation.  In this case, the second dimension is 
%    the same as that of Y.  If the number of arguments also varies from
%    from observation to observation, the second dimension of ARGVALS should
%    equal the the largest number of argument values, and the elements in 
%    each show should be padded out with NaN.
%    Argument values falling outside of the range specified in the
%    BASIS and their corresponding values in Y will not be used, 
%    but if this happens, a warning message is displayed.
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
%    Earlier releases of DATA2FD supplied additional arguments for constructing
%    a default basis, but these have been eliminated in favor of using new
%    function MAKE_BASIS.
%    
%  FDNAMES  ... (optional)
%    A cell object of length 3 with each element being a string specifying
%    a name for each dimension of Y.
%       1. argument domain, such as the string 'Time'
%       2. replications or cases
%       3. variables
%    For example, for the monthly temperature/precipitation data, 
%    fdnames{1} = 'Month';
%    fdnames{2} = 'Station';
%    fdnames{3} = 'Variable';
%    By default, the string 'time', 'reps' and 'values' are used.
%
%  DATA2FD Returns the object FDOBJ of functional data class containing 
%    coefficients for the expansion and the functional data basis object 
%    BASIS.
%
%  DATA2FD is intended for more casual use not requiring a great deal of 
%    control over the smoothness of the functional data object.  It uses
%    function PROJECT_BASIS to compute the functional data object. Indeed,
%    in the simplest and most common situation, DATA2FD consists of 
%             coef  = project_basis(y, argvals, basis);
%             fdobj = fd(coef,basis,fdnames);
%    However, for more advanced applications requiring more smoothing 
%    control than is possible by setting the number of basis functions in
%    BASIS, function SMOOTH_BASIS should be used.  Or, alternatively,
%    DATA2FD may first be used with a generous number of basis functions,
%    followed by smoothing using function SMOOTH.

%  Last modified:  4 September 2007

%  set default names for dimensions of Y in FDNAMES

 if nargin < 4
    fdnames{1} = 'arguments';
    fdnames{2} = 'replications';
    fdnames{3} = 'functions';
 end
 
 %  check the BASIS argument
 
 if nargin < 3
    error('Argument BASIS must be supplied in this release.');
 end
 if ~isa_basis(basisobj)
    error('Argument BASISOBJ is not a basis obj.');
 end
 
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
 %  set number of replications
 if (ydim > 1)
    nrep = yd(2);
 else
    nrep = 1;
 end
 %  set number of variables
 if (ydim > 2)
    nvar = yd(3);
 else
    nvar = 1;
 end

 %  check argument value array

 argd  = size(argvals);
 if length(argd) > 2
    error('Argument ARGVALS has too many dimensions.');
 end
 if argd(1) == 1, argvals = argvals'; argd = size(argvals); end
 nargd = length(argd);
 if argd(nargd) == 1, nargd = nargd - 1; end
 if (argd(1) ~= n)
    error('Number of arg. values not equal to 1st dim. of Y.');
 end
 if nargd == 2 && argd(2) ~= 1 && argd(2) ~= nrep
    error(['Matrix ARGVALS must have same number of columns', ...
           'as the number of replicates.']);
 end
 %  Issue a warning if arguments are outside of the range in the basis.
 rangeval = getbasisrange(basisobj);
 if nargd == 1
    temp = argvals;
 else
    temp = reshape(argvals, n*nrep, 1);
 end
 temp = temp(~isnan(argvals));
 if min(temp) < rangeval(1) || max(temp) > rangeval(2)
    warning('Wid1:range', ...
            'Some arguments values are outside of the range in BASIS.');
    if nargd == 1
       index = argvals < rangeval(1) || argvals > rangeval(2);
       argvals(index) = NaN;
    else
       for irep=1:nrep
          index = argvals(:,irep) < rangeval(1) || ...
                  argvals(:,irep) > rangeval(2);
          argvals(index,irep) = NaN;
       end
    end
 end
 
 %  Determine which of three cases applies:
 %  First case:  ARGVALS a vector, and there are no missing values in Y.
 %  Second case: ARGVALS a vector, but missing values in ARGVALS and/or Y. 
 %  Third case:  ARGVALS is a matrix

 if nargd == 1 || (nargd == 2 && argd(2) == 1) 
    if ~any(isnan(argvals)) && ~any(isnan(y(:)))
      %  ------------------------------------------------
      %  First case:  no missing values, ARGVALS a vector
      %  ------------------------------------------------
      nbasis = getnbasis(basisobj);
      if nbasis <= n
        coef = project_basis(y, argvals, basisobj);
      else
        coef = project_basis(y, argvals, basisobj, 1);
      end

   else
      %  ---------------------------------------------------------
      %  Second case:  ARGVALS a vector, but missing data present. 
      %  ---------------------------------------------------------
      nbasis   = getnbasis(basisobj);
      coefd    = yd;
      coefd(1) = nbasis;
      coef     = zeros(coefd);
      % set up penalty and basis matrix
      index    = ~isnan(argvals);
      basismat = getbasismatrix(argvals, basisobj);
      penmat   = eval_penalty(basisobj);
      penmat   = penmat + 1e-10 .* max(max(penmat)) .* eye(nbasis);
      lambda1  = 0.0001 .* sum(sum(basismat(index,:).^2))./sum(sum(diag(penmat)));
      penmat1  = lambda1 .* penmat;
      %  process each observation in turn
      if(length(coefd) == 2)
        for j=1:nrep
          yy    = y(:,j);
          index = ~isnan(yy) && ~isnan(argvals);
          if length(yy(index)) < 2
             error(['Less than 2 data values available for curve ',...
                    num2str(j),'.']);
          end
          temp  = basismat(index,:);
          Cmat  = temp' * temp + penmat1;
          Dmat  = temp' * yy(index);
          coef(:,j) = symsolve(Cmat, Dmat);
        end
      else
        for j=1:nrep
          for k=1:nvar
            yy = y(:,j,k);
            index = ~isnan(yy) && ~isnan(argvals);
            if length(yy(index)) < 2
               error(['Less than 2 data values available for curve ',...
                      num2str(j),' and variable ',num2str(k),'.']);
            end
            temp  = basismat(index,:);
            Cmat = temp' * temp + penmat1;
            Dmat = temp' * yy(index);
            coef(:,j,k) = symsolve(Cmat,Dmat);
          end
        end
      end
    end
 else
    %  ------------------------------------------------
    %  Third case:  ARGVALS a matrix. 
    %  ------------------------------------------------
    nbasis   = getnbasis(basisobj);
    coefd    = yd;
    coefd(1) = nbasis;
    coef     = zeros(coefd);
    argv     = reshape(argvals,n*nrep,1);
    index    = ~isnan(argv);
    argv     = unique(argv(index));
    basismat = getbasismatrix(argv, basisobj);
    penmat   = eval_penalty(basisobj);
    penmat   = penmat + 1e-10 .* max(max(penmat)) .* eye(nbasis);
    lambda1  = 0.0001 .* sum(sum(basismat.^2))./sum(sum(diag(penmat)));
    penmat1  = lambda1 .* penmat;
    if length(coefd) == 2 
      for j=1:nrep
        yy    = y(:,j);
        argv  = argvals(:,j);
        index = (~isnan(yy) & ~isnan(argv));
        if length(yy(index)) < 2
           error(['Less than 2 data values available for curve ',...
                   num2str(j),'.']);
        end
        coef(:,j) =  ...
            project_basis(yy(index), argv(index), basisobj, penmat1);
      end
    else
      for j=1:nrep
        for k=1:nvar
          yy    = y(:,j,k);
          argv  = argvals(:,j);
          index = (~isnan(yy) & ~isnan(argv));
          if length(yy(index)) < 2
              error(['Less than 2 data values available for curve ',...
                      num2str(j),' and variable ',num2str(k),'.']);
          end
          coef(:,j,k) = ...
            project_basis(yy(index), argv(index), basisobj, penmat1);
        end
      end
    end
 end
 
 % set up the functional data object
 
 fdobj = fd(coef, basisobj, fdnames);


