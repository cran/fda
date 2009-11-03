function basisobj = create_nspline_basis(rangeval, nbasis, norder, ...
                                         breaks, dropind)
%CREATE_NPLINE_BASIS Creates a natural spline functional data basis.
%
%  Arguments ...
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS   ... the number of basis functions
%  NORDER   ... order of n-splines (one higher than their degree).  The
%                 default of 4 gives (natural) cubic splines.
%  BREAKS   ... also called knots, these are a strictly increasing sequence
%               of junction points between piecewise polynomial segments.
%               They must satisfy BREAKS(1) = RANGEVAL(1) and
%               BREAKS(NBREAKS) = RANGEVAL(2), where NBREAKS is the total
%               number of BREAKS.
%  There is a potential for inconsistency among arguments NBASIS, NORDER, and
%  BREAKS.  It is resolved as follows:
%     If BREAKS is supplied, NBREAKS = length(BREAKS), and
%     NBASIS = NBREAKS + NORDER - 4, no matter what value for NBASIS is
%     supplied.
%     If BREAKS is not supplied but NBASIS is, NBREAKS = NBASIS - NORDER + 2,
%        and if this turns out to be less than 3, an error message results.
%     If neither BREAKS nor NBASIS is supplied, NBREAKS is set to 21.
%  DROPIND ... a set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];
%  Returns
%  BASISOBJ  ... a functional data basis object

%  added by Kris Villez in August 2011 based on bsplinepen file in the
%  FDA toolbox by Jim Ramsay

%  last modified 31 October 2011 by Kris Villez: changed terminology


%  Default basis for missing arguments

if nargin==0
    type        = 'nspline';
    rangeval    = [0,1];
    nbasis      = 2;
    params      = [];
    dropind     = [];
    quadvals    = [];
    values      = {};
    basisvalues = {};

    basisobj = basis(type, rangeval, nbasis, params, ...
                     dropind, quadvals, values, basisvalues);
    return
end

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

if nargin < 2, nbasis = 4               ;               end
if nargin < 3, norder = min(4,nbasis)   ;               end
if nargin < 4
    breaks = [];
else
    if size(breaks,1) > 1, breaks = breaks';  end
    if size(breaks,1) > 1
        error('BREAKS must be a vector.');
    end
end
if nargin < 5, dropind = [];           end

%  Determine what to do if some arguments are empty

% case of empty NBASIS and empty BREAKS: set up splines with a
%  single interior knot

if isempty(nbasis) && isempty(breaks)
    nbasis  = 5;
    norder  = 4;
    nbreaks = 3;
    breaks  = linspace(rangeval(1), rangeval(2), nbreaks);
end

% If NBASIS is empty but BREAKS are supplied, determine
%   NBASIS = NORDER + number of interior knots.
if isempty(nbasis) && ~isempty(breaks)
    nbreaks = length(breaks);
    nbasis  = nbreaks + norder - 4;
end
% If NBASIS and NORDER are present but no BREAKS supplied,
%   set up NBASIS - NORDER + 2 equally spaced breaks.
if ~isempty(nbasis) && isempty(breaks)
    nbreaks = nbasis - norder + 4; 
    breaks  = linspace(rangeval(1), rangeval(2), nbreaks);
end

%  Special argument configurations taken care of.  
%  Now go ahead and set up the basis

nbreaks = length(breaks);

%  check the compatibility of NBASIS, NBREAKS and RANGEVAL

if (nbreaks < 2)
    error ('Number of values in BREAKS less than 2.');
end

if (nbasis < nbreaks-1)
    error ('NBASIS is less than number of values=BREAKS.');
end
if breaks(1) ~= rangeval(1)
    error('Smallest value in BREAKS not equal to RANGEVAL(1).');
end
if breaks(nbreaks) ~= rangeval(2)
    error('Largest  value in BREAKS not equal to RANGEVAL(2).');
end

%  The PARAMS field contains only the interior knots; drop end breaks

if nbreaks > 2
    params   = breaks(2:(nbreaks-1));
else
    params = [];
end

%  check DROPIND

if length(dropind) > 0
    if length(dropind) >= nbasis
        error('Too many index values in DROPIND.');
    end
    dropind = sort(dropind);
    if length(dropind) > 1
        if any(diff(dropind)) == 0
            error('Multiple index values in DROPIND.');
        end
    end
    for i=1:length(dropind);
        if dropind(i) < 1 || dropind(i) > nbasis
            error('An index value is out of range.');
        end
    end
end

%  construct basis object

type        = 'nspline';
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);
