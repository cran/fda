function isabasis = isa_basis(basisobj)
%  ISA_BASIS  checks a struct object for fields for basis objects

%  last modified 12 September 2004

isabasis = 1;
if ~strcmp(class(basisobj),'basis')
    isabasis = 0;
end

