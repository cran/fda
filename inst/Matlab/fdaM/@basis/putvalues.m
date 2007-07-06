function newbasisobj = putvalues(basisobj, values)
%  PUTVALUES   Enters a cell array of values of basis functions
%  and a number of their derivatives into basis object BASISOBJ 

%  last modified 18 June 2007

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

%  check values

if isempty(basisobj.quadvals)
    error('QUADVALS is empty.');
else
    nquad = size(basisobj.quadvals,1);
end

if ~iscell(values)
    error('VALUES argument is not a cell array.');
end
if ~isempty(values{1})
    nvalues = length(values);
    for ivalue=1:nvalues
        [n,k] = size(full(values{ivalue}));
        if n ~= nquad
            error(['Number of rows in VALUES not equal to ', ...
                    'number of quadrature points.']);
        end
        if k ~= basisobj.nbasis
            error(['Number of columns in VALUES not equal to ', ...
                    'number of basis functions.']);
        end
    end
end

newbasisobj.type        = basisobj.type;
newbasisobj.rangeval    = basisobj.rangeval;
newbasisobj.nbasis      = basisobj.nbasis;
newbasisobj.params      = basisobj.params;
newbasisobj.dropind     = basisobj.dropind;
newbasisobj.quadvals    = basisobj.quadvals;
newbasisobj.values      = values;
newbasisobj.basisvalues = basisobj.basisvalues;

newbasisobj = class(newbasisobj, 'basis');