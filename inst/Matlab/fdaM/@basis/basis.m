function basisobj = basis(basistype, rangeval, nbasis, params, ...
                          dropind, quadvals, values, basisvalues)
%  BASIS  Creates a functional data basis.
%  Arguments or slots:
%  BASISTYPE ... a string indicating the type of basis.  This may be one of
%               'Bspline', 'bspline', 'Bsp', 'bsp',
%               'con', 'const', 'constant'
%               'exp', 'exponen', 'exponential'
%               'Fourier', 'fourier', 'Fou', 'fou',
%               'mon', 'monom', 'monomial',
%               'polyg' 'polygon', 'polygonal'
%               'pol', 'poly', 'polynomial'
%               'power', 'pow'
%               'QW'
%               'QWM'
%               'QS'
%               'slide'
%               'fd'
%               'FEM'
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values
%  NBASIS   ... the number of basis functions
%  PARAMS   ... If the basis is 'fourier', this is a single number indicating
%                 the period.  That is, the basis functions are periodic on
%                 the interval (0,PARAMS) or any translation of it.
%               If the basis is 'bspline', the values are interior points at
%                 which the piecewise polynomials join.
%                 Note that the number of basis functions NBASIS is equal
%                 to the order of the Bspline functions plus the number of
%                 interior knots, that is the length of PARAMS.
%               This means that NBASIS must be at least 1 larger than the
%                 length of PARAMS.
%  DROPIND ... A set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];  Note
%                that argument NBASIS is reduced by the number of indices,
%                and the derivative matrices in VALUES are also clipped.
%  QUADVALS .. A NQUAD by 2 matrix.  The first column contains quadrature
%                points to be used in a fixed point quadrature.  The second
%                contains quadrature weights.  For example, for Simpson's
%                rule for NQUAD = 7, the points are equally spaced and the
%                weights are delta.*[1, 4, 2, 4, 2, 4, 1]/3.  DELTA is the
%                spacing between quadrature points.  The default is [].
%  VALUES  ... A cell array, with entries containing the values of
%                the basis function derivatives starting with 0 and
%                going up to the highest derivative needed.  The values
%                correspond to quadrature points in QUADVALS and it is
%                up to the user to decide whether or not to multiply
%                the derivative values by the square roots of the
%                quadrature weights so as to make numerical integration
%                a simple matrix multiplication.
%                Values are checked against QUADVALS to ensure the correct
%                number of rows, and against NBASIS to ensure the correct
%                number of columns.
%                The default is VALUES{1} = [];
%  BASISVALUES ... A cell array.  The cell array must be 2-dimensional,
%                with a variable number of rows and two or more columns.
%                This field is designed to avoid evaluation of a
%                basis system repeatedly at a set of argument values.
%                Each row corresponds to a specific set of argument values.
%                The first  cell in that row contains the argument values.
%                The second cell in that row contains a matrix of values of
%                the basis functions.
%                The third and subsequent cells contain matrices of values
%                their derivatives up to a maximum derivative order.
%                Whenever function getbasismatrix is called, it checks
%                the first cell in each row to see, first, if the number of
%                argument values corresponds to the size of the first dimension,
%                and if this test succeeds, checks that all of the argument
%                values match.  This takes time, of course, but is much
%                faster than re-evaluation of the basis system.  Even this
%                time can be avoided by direct retrieval of the desired
%                array.
%
%  Returns
%  BASISOBJ  ... a basis object with slots
%         type
%         rangeval
%         nbasis
%         params
%         dropind
%         quadvals
%         values
%         basisvalues
%  Slot VALUES contains values of basis functions and derivatives at
%   quadrature points weighted by square root of quadrature weights.
%   These values are only generated as required, and only if slot
%   quadvals is not empty.
%
%  An alternative name for this function is CREATE_BASIS, but PARAMS argument
%     must be supplied.
%  Specific types of bases may be set up more conveniently using functions
%  CREATE_BSPLINE_BASIS    ...  creates a b-spline basis
%  CREATE_CONSTANT_BASIS   ...  creates a constant basis
%  CREATE_FOURIER_BASIS    ...  creates a fourier basis
%  CREATE_MONOM_BASIS      ...  creates a monomial basis
%  CREATE_POLYGON_BASIS    ...  creates a polygonal basis
%  CREATE_POLYNOMIAL_BASIS ...  creates a polynomial basis
%  CREATE_POWER_BASIS      ...  creates a polygonal basis
%  CREATE_QW_BASIS         ...  creates a Weibull W basis
%  CREATE_QWM_BASIS        ...  creates a modified Weibull W basis
%  CREATE_SLIDE_BASIS      ...  creates a slide basis
%  CREATE_FD_BASIS         ...  creates a functional data object basis
%  CREATE_FEM_BASIS        ...  creates a finite element basis

%  Last modified 25 May 2010

%  Set up default basis if there are no arguments

if nargin==0
    basisobj.type        = 'bspline';
    basisobj.rangeval    = [0,1];
    basisobj.nbasis      = 2;
    basisobj.params      = [];
    basisobj.dropind     = [];
    basisobj.quadvals    = [];
    basisobj.values      = {};
    basisobj.basisvalues = {};
    basisobj = class(basisobj, 'basis');
    return;
end

%  If arguments are supplied, at least the first four must be supplied.

if nargin < 4
    error('Less than four arguments found.');
end

%  if first argument is a basis object, return

if isa(basistype, 'basis')
    basisobj = basistype;
    return;
end

%  check basistype

basistype = use_proper_basis(basistype);
if strcmp(basistype,'unknown')
    error ('TYPE unrecognizable.');
end

%  check if QUADVALS is present, and set to default if not

if nargin < 6
    quadvals = [];
else
    if ~isempty(quadvals)
        [nquad,ncol] = size(quadvals);
        if nquad == 2 && ncol > 2
            quadvals = quadvals';
            [nquad,ncol] = size(quadvals);
        end
        if nquad < 2
            error('Less than two quadrature points are supplied.');
        end
        if ncol ~= 2
            error('QUADVALS does not have two columns.');
        end
    end
end

%  check VALUES if present, and set to a single empty cell
%  if not.

if nargin < 7
    values = {};
else
    if ~iscell(values)
        error('VALUES argument is not a cell array.');
    end
    if ~isempty(values)
        if ~isempty(values{1})
            nvalues = length(values);
            for ivalue=1:nvalues
                [n,k] = size(full(values{ivalue}));
                if n ~= nquad
                    error(['Number of rows in VALUES not equal to ', ...
                        'number of quadrature points.']);
                end
                if k ~= nbasis
                    error(['Number of columns in VALUES not equal to ', ...
                        'number of basis functions.']);
                end
            end
        end
    end
end

%  check BASISVALUES are present, and if not set to empty cell array

if nargin < 8
    basisvalues = {};
else
    if ~isempty(basisvalues)
        %  check BASISVALUES
        if ~iscell(basisvalues)
            error('BASISVALUES is not a cell object.')
        end
        sizevec = size(basisvalues);
        if length(sizevec) > 2
            error('BASISVALUES is not 2-dimensional.')
        end
        for i=1:sizevec(1)
            if length(basisvalues{i,1}) ~= size(basisvalues{i,2},1)
                error(['Number of argument values not equal number ', ...
                       'of values.']);
            end
        end
    end
end

%  check if DROPIND is present, and set to default if not

if nargin < 5
    dropind = [];
else
    if ~isempty(dropind)
        %  check DROPIND
        ndrop = length(dropind);
        if ndrop >= nbasis
            error('Too many index values in DROPIND.');
        end
        dropind = sort(dropind);
        if ndrop > 1
            if any(diff(dropind)) == 0
                error('Multiple index values in DROPIND.');
            end
        end
        for i=1:ndrop;
            if dropind(i) < 1 || dropind(i) > nbasis
                error('An index value is out of range.');
            end
        end
        %  drop columns from VALUES cells if present
        droppad = [dropind,zeros(1,nbasis-ndrop)];
        keepind = (1:nbasis) ~= droppad;
        if ~isempty(values) && ~isempty(values{1})
            for ivalue=1:nvalues
                derivvals = values{ivalue};
                derivvals = derivvals(:,keepind);
                values{ivalue} = derivvals;
            end
        end
    end
end

%  select the appropriate type and process

switch basistype
    case 'fourier'
        period     = params(1);
        if (period <= 0)
            error ('Period must be positive for a Fourier basis');
        end
        params = period;
        if (2*floor(nbasis/2) == nbasis)
            nbasis = nbasis + 1;
        end

    case 'bspline'
        if ~isempty(params)
            nparams  = length(params);
            if (params(1) <= rangeval(1))
                error('Smallest value in BREAKS not within RANGEVAL');
            end
            if (params(nparams) >= rangeval(2))
                error('Largest value in BREAKS not within RANGEVAL');
            end
        end

    case 'expon'
        if (length(params) ~= nbasis)
            error(['No. of parameters not equal to no. of basis fns ',  ...
                    'for exponential basis.']);
        end

    case 'polyg'
        if (length(params) ~= nbasis)
            error(...
                'No. of parameters not equal to no. of basis fns for polygonal basis.');
        end

    case 'power'
        if length(params) ~= nbasis
            error(...
                'No. of parameters not equal to no. of basis fns for power basis.');
        end

    case 'const'
        params = 0;

    case 'monom'
        if length(params) ~= nbasis
            error(['No. of parameters not equal to no. of basis fns', ...
                   ' for monomial basis.']);
        end

    case 'polynom'
        if length(params) > 1
            error('More than one parameter for a polynomial basis.');
        end

    case 'QW'
        if ~isempty(params)
            error('More than zero parameters for a QW basis.');
        end
    case 'QWM'
        if ~isempty(params)
            error('More than zero parameters for a QWM basis.');
        end
    case 'QS'
        if ~isempty(params)
            error('More than zero parameters for a QWM basis.');
        end
    case 'slide'
        if length(params) ~= 2*nbasis-1
            error('Number of parameters not correct for a slide basis.');
        end
    case 'fd'
        if ~strcmp(class(params), 'fd')
            error('Parameter not a functional data object')
        end
    case 'FEM'
        if ~strcmp(class(params), 'struct')
            error('Parameter not a struct object')
        end
    otherwise
        error('Unrecognizable basis');
end

basisobj.type        = basistype;
basisobj.rangeval    = rangeval;
basisobj.nbasis      = nbasis;
basisobj.params      = params;
basisobj.dropind     = dropind;
basisobj.quadvals    = quadvals;
basisobj.values      = values;
basisobj.basisvalues = basisvalues;

basisobj = class(basisobj, 'basis');

%  ------------------------------------------------------------------------

function fdtype = use_proper_basis(fdtype)
%  USE_PROPER_BASIS recognizes type of basis by use of several variant spellings

switch fdtype

    %  B-spline basis

    case 'bspline'
        fdtype = 'bspline';
    case 'Bspline'
        fdtype = 'bspline';
    case 'Bsp'
        fdtype = 'bspline';
    case 'bsp'
        fdtype = 'bspline';

    %  constant basis

    case 'con'
        fdtype = 'const';
    case 'const'
        fdtype = 'const';
    case 'const'
        fdtype = 'const';

    %  exponential basis

    case 'exp'
        fdtype = 'expon';
    case 'expon'
        fdtype = 'expon';
    case 'exponential'
        fdtype = 'expon';

    % Fourier basis

    case 'Fourier'
        fdtype = 'fourier';
    case 'fourier'
        fdtype = 'fourier';
    case 'Fou'
        fdtype = 'fourier';
    case 'fou'
        fdtype = 'fourier';

    %  Monomial basis

    case 'mon'
        fdtype = 'monom';
    case 'monom'
        fdtype = 'monom';
    case 'monomial'
        fdtype = 'monom';

    %  Polygonal basis

    case 'polyg'
        fdtype = 'polyg';
    case 'polygon'
        fdtype = 'polyg';
    case 'polygonal'
        fdtype = 'polyg';

    %  Polynomial basis

    case 'poly'
        fdtype = 'polynom';
    case 'polynom'
        fdtype = 'polynom';
    case 'polynomial'
        fdtype = 'polynom';

    %  Power basis

    case 'power'
        fdtype = 'power';
    case 'pow'
        fdtype = 'power';

    %  Weibull W basis

    case 'QW'
        fdtype = 'QW';

    %  modified Weibull W basis

    case 'QWM'
        fdtype = 'QWM';

    case 'QS'
        fdtype = 'QS';
        
    %  slide basis    
        
    case 'slide'
        fdtype = 'slide';

    %  functional data object basis    
        
    case 'fd'
        fdtype = 'fd';

    %  FEM basis    
        
    case 'FEM'
        fdtype = 'FEM';

    %  Not a recognizable basis

    otherwise
        fdtype = 'unknown';

end
