function basisobj = create_bspline_basis(rangeval, nbasis, norder, ...
                                         breaks, dropind, quadvals, values)
%CREATE_SPLINE_BASIS Creates a bspline functional data basis.
%  This function is identical to BSPLINE_BASIS.
%  Arguments ...
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a 
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS   ... the number of basis functions
%  NORDER   ... order of b-splines (one higher than their degree).  The
%                 default of 4 gives cubic splines.
%  BREAKS   ... also called knots, these are a strictly increasing sequence
%               of junction points between piecewise polynomial segments.
%               They must satisfy BREAKS(1) = RANGEVAL(1) and
%               BREAKS(NBREAKS) = RANGEVAL(2), where NBREAKS is the total
%               number of BREAKS.  
%  There is a potential for inconsistency among arguments NBASIS, NORDER, and
%  BREAKS.  It is resolved as follows:
%     If BREAKS is supplied, NBREAKS = length(BREAKS), and
%     NBASIS = NBREAKS + NORDER - 2, no matter what value for NBASIS is
%     supplied.
%     If BREAKS is not supplied but NBASIS is, NBREAKS = NBASIS - NORDER + 2,
%        and if this turns out to be less than 3, an error message results.
%     If neither BREAKS nor NBASIS is supplied, NBREAKS is set to 21.
%  DROPIND ... a set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];
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
%
%  Returns
%  BASISOBJ  ... a functional data basis object

%  A B-spline basis may also be constructed using CREATE_EASY_BASIS 
%    CREATE_BASIS or MAKE_BASIS.

%  last modified 20 July 2006

%  check RANGEVAL

if length(rangeval) == 1
    if rangeval <= 0
        error('RANGEVAL a single value that is not positive.');
    end
    rangeval = [0,rangeval];
end

if rangechk(rangeval) ~= 1
    error('RANGEVAL is not a legitimate range.');
end

%  set some default values

%  number of basis functions
if nargin < 2
    nbasis = 1;
end
%  order of splines
if nargin < 3
    norder = min(4,nbasis);
end
%  knots
if nargin < 4
    breaks = [];
else
    if size(breaks,1) > 1, breaks = breaks';  end
    if size(breaks,1) > 1
        error('BREAKS must be a vector.');
    end
end
%  indices of basis functions to be dropped
if nargin < 5
    dropind = [];
end
%  quadrature points and weights
if nargin < 6
    quadvals = [];
end
%  basis derivative values at quadrature points
if nargin < 7
    values{1} = [];
end


type = 'bspline';

%  Determine what to do if some arguments are empty

% If both NBASIS and BREAKS are missing, but NORDER IS NOT,
%   use 21 equally spaced knots, and determinine NBASIS.
%   by NBASIS = NORDER + 19.
if isempty(nbasis) && isempty(breaks)
    nbreaks = 21;
    nbasis  = 19 + norder;
    breaks  = linspace(rangeval(1), rangeval(2), nbreaks);
end
% If NBASIS is empty but BREAKS are supplied, determine
%   NBASIS = NORDER + number of interior knots.
if isempty(nbasis) && ~isempty(breaks)
    nbreaks = length(breaks);
    nbasis  = nbreaks + norder - 2;
end
% If NBASIS and NORDER are present but no BREAKS supplied,
%   set up NBASIS - NORDER + 2 equally spaced breaks.
if ~isempty(nbasis) && isempty(breaks)
    nbreaks = nbasis - norder + 2;
    breaks  = linspace(rangeval(1), rangeval(2), nbreaks);
end
nbreaks = length(breaks);

%  check the compatibility of NBASIS, NBREAKS and RANGEVAL

if (nbreaks < 2)
    error ('Number of values in BREAKS less than 2.');
end
if (nbasis < nbreaks-1)
    error ('NBASIS is less than number of values=BREAKS.');
end
if (breaks(1) ~= rangeval(1))
    error('Smallest value in BREAKS not equal to RANGEVAL(1).');
end
if (breaks(nbreaks) ~= rangeval(2))
    error('Largest  value in BREAKS not equal to RANGEVAL(2).');
end

%  The PARAMS field contains only the interior knots

if nbreaks > 2
    params   = breaks(2:(nbreaks-1));
else
    params = [];
end

%  set default values

%  construct basis object

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values);
