function newbasisobj = putdropind(basisobj, dropind)
%  PUQUADVALS   Enters drop indices for
%     basis object BASISOBJ into slot basisobj.dropind

%  last modified 18 June 2007

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

%  check DROPIND

if length(dropind) >= basisobj.nbasis
    error('Too many index values in DROPIND.');
end
dropind = sort(dropind);
if length(dropind) > 1
    if any(diff(dropind)) == 0
        error('Multiple index values in DROPIND.');
    end
end
for i=1:length(dropind);
    if dropind(i) < 1 || dropind(i) > basisobj.nbasis
        error('An index value is out of range.');
    end
end

newbasisobj.type        = basisobj.type;
newbasisobj.rangeval    = basisobj.rangeval;
newbasisobj.nbasis      = basisobj.nbasis;
newbasisobj.params      = basisobj.params;
newbasisobj.dropind     = dropind;
newbasisobj.quadvals    = basisobj.quadvals;
newbasisobj.values      = basisobj.values;
newbasisobj.basisvalues = basisobj.basisvalues;

newbasisobj = class(newbasisobj, 'basis');

