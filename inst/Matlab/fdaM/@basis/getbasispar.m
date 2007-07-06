function params = getbasispar(basisobj)
%  GETBASISPAR   Extracts the basis parameters from basis object BASISOBJ.

%  last modified 16 March 1999

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

params = basisobj.params;
