function basisobj = create_polygon_basis(argvals)
%  CREATE_POLYGON_BASIS Creates a polygonal basis
%  Argument:
%  ARGVALS  ... strictly increasing argument values
%  Return:
%  BASIS_FD ... a functional data basis object of type 'polygon'

%  last modified 18 June 2007

%  check that argument values are strictly increasing

if min(diff(argvals)) <= 0
    error('ARGVALS are not strictly increasing.');
end

type     = 'polyg';
nbasis   = length(argvals);
rangeval = [min(argvals), max(argvals)];
params   = argvals;

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

dropind     = [];
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);
