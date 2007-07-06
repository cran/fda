function newbasisobj = putbasisvalues(basisobj, basisvalues)
%  PUTVALUES   Enters a cell array of values of basis functions
%  and a number of their derivatives into basis object BASISOBJ 

%  last modified 18 June 2007

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

%  check values

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
            error(['Number of argument valuesdoes not equal ', ...
                   'number of values.']);
        end
    end
end

newbasisobj.type        = basisobj.type;
newbasisobj.rangeval    = basisobj.rangeval;
newbasisobj.nbasis      = basisobj.nbasis;
newbasisobj.params      = basisobj.params;
newbasisobj.dropind     = basisobj.dropind;
newbasisobj.quadvals    = basisobj.quadvals;
newbasisobj.values      = basisobj.values;
newbasisobj.basisvalues = basisvalues;

newbasisobj = class(newbasisobj, 'basis');