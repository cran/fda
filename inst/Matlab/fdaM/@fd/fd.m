function fdobj = fd(coef, basisobj, fdnames)
%  FD   Creates a functional data object.
%  A functional data object consists of a basis for expanding a functional
%    observation and a set of coefficients defining this expansion.
%    The basis is contained=a 'basis' object; that is, a realization
%    of the 'basis' class.
%
%  Arguments:
%  COEF     ... An array containing coefficient values for the expansion
%               of each set of function values=terms of a set of basis
%               functions.
%               If COEF is a three-way array, then the first dimension
%               corresponds to basis functions, the second to replications,
%               and the third to variables.
%               If COEF is a matrix, it is assumed that there is only
%               one variable per replication, and then:
%                 rows    correspond to basis functions
%                 columns correspond to replications
%               If COEF is a vector, it is assumed that there is only one
%               replication and one variable.
%               If COEF is empty, it is left that way.
%  BASISOBJ ... a functional data basis object.
%               If this argument is missing, a B-spline basis is created
%               with the number of basis functions equal to the size
%               of the first dimension of coef.
%  FDNAMES  ... A cell array of length 3 with members containing
%               1. a name for the argument domain, such as 'Time'
%               2. a name for the replications or cases
%               3. a name for the function
%               If this argument is not supplied, the strings
%               'arguments', 'replications' and 'functions' are used
%               Each of the cells may itself be a cell array of length 2
%               in which case the first cell contains the name as above
%               for the dimension of the data, and the second cell 
%               contains a character matrix of names for each index
%               value.  Note that the rows must be of the same length,
%               so that variable-length names must be padded out with
%               blanks as required
%               For example, here is the setup for the multivariate
%               juggling data where the variables have names 
%               'X', 'Y', and 'Z':
%               jugglenames = cell(1,3);
%               jugglenames{1} = 'Seconds';
%               jugglenames{2} = 'Record';
%               varnames = cell(1,2);
%               varnames{1} = 'Coordinate';
%               varnames{2} = ['X'; 'Y'; 'Z'];
%               jugglenames{3} = varnames;

%
%  An alternative argument list:
%  The argument COEF can be dropped, so that BASISOBJ is the
%  leading argument.  In this case, COEF is set to [], an empty
%  array.  For many purposes, and especially for functional parameter
%  objects (fdPar), the coefficient array is either not needed, or
%  supplied later.
%
%  Returns:
%  FD ... a functional data object

%  last modified 24 December 2011

%  Set default FDNAMES

    defaultfdnames{1} = 'arguments';
    defaultfdnames{2} = 'replications';
    defaultfdnames{3} = 'functions';

%  Define the default FD object

if nargin == 0

    %  set up default object if there are no arguments:
    %    a single function with a constant basis and coefficient zero

    fdobj.coef     = zeros(1,1);
    fdobj.basisobj = create_constant_basis([0,1]);
    fdobj.fdnames  = defaultfdnames;
    fdobj = class(fdobj, 'fd');
    return;

end

%  set default arguments

if nargin < 3,  fdnames  = defaultfdnames;  end
if nargin < 2
    error('Less than two arguments.');         
end

%  check leading argument COEF

if ~strcmp(class(coef), 'double')
    error('Leading argument is not of class double.');
end

%  check dimensions of COEF

if ~isempty(coef)

    %  check dimensions of coefficient array

    coefd = size(coef);
    ndim  = length(coefd);
    if (ndim > 3)
        error('Coefficient array has more than three dimensions.');
    end
else
    error('Coefficient matrix is empty.');
end

%  check BASISOBJ

if ~strcmp(class(basisobj), 'basis')
    error('Argument BASISOBJ is not of basis class.');
end

nbasis  = getnbasis(basisobj);

if ~isempty(coef) && ~strcmp(getbasistype(basisobj),'fdVariance')
    if coefd(1) ~= nbasis
        error(['First dimension of coefficient array is ', ...
                'not equal to number of basis functions.']);
    end
end

%  check fdnames

if ~iscell(fdnames)
    error('FDNAMES not a cell object.');
end

if length(fdnames) > 3
    error('FDNAMES has length greater than 3.');
end

%  fill missing cells with default labels

if length(fdnames) == 2
    fdnames{3} = 'functions';
end
if length(fdnames) == 1
    fdnames{2} = 'replications';
    fdnames{3} = 'functions';
end

%  set up object

fdobj.coef     = coef;
fdobj.basisobj = basisobj;
fdobj.fdnames  = fdnames;

fdobj = class(fdobj, 'fd');

