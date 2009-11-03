function basisobj = create_polynomial_basis(rangeval, nbasis, ctr, dropind)
%  CREATE_POLY_BASIS  Creates a monomial basis:, 1, x, ..., x^{nbasis-1}
%  Arguments:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a 
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS   ... the number of basis functions
%  CTR      ... If 1, center the range before evaluating.
%  DROPIND ... a set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];
%  Return:
%  BASIS.FD  ... a functional data basis object of type 'constant'

%  last modified 20 July 20064

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

if nargin < 3, ctr = 0;          end
if nargin < 2, nbasis = 2;       end
if nargin < 1, rangeval = [0,1]; end
type    = 'polynom';
params  = ctr;

%  check DROPIND

if nargin < 5
    dropind = [];
end

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

dropind   = [];
quadvals  = [];
values{1} = [];

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values);

