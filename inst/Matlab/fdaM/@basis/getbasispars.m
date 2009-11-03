function params = getbasispars(basisobj)
%  GETBASISPARS  Extracts the basis parameters from basis object BASISOBJ.
%  It is an exact copy of function GETBASISPAR.

%  last modified 16 March 1999

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

params = basisobj.params;
