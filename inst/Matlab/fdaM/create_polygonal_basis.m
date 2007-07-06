function basisobj = create_polygonal_basis(rangeval, argvals, dropind)
%  CREATE_POLYGON_BASIS Creates a polygonal basis
%  Argument:
%  ARGVALS  ... strictly increasing argument values
%  Return:
%  BASIS_FD ... a functional data basis object of type 'polygon'

%  last modified 3 January 2008

%  default RANGEVAL

if nargin < 1, rangeval = [0,1];  end

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

%  set up default arguments

if nargin < 2, argvals = rangeval; end
if nargin < 3, dropind = [];       end

%  check ARGVALS

if min(diff(argvals)) <= 0
    error('ARGVALS are not strictly increasing.');
end

if min(argvals) < rangeval(1) || max(argvals) > rangeval(2)
    error('ARGVALS are out of range.');
end

%  make sure range is included

argvals = sort(unique([rangeval(:); argvals(:)]));
nbasis  = length(argvals);

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

type        = 'polyg';
nbasis      = length(argvals);
params      = argvals';
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);
