function dropind = getdropind(basisobj)
%  GETDROPIND   Extracts the indices of basis functions to be dropped
%     from basis object BASISOBJ.

%  last modified 6 April 2004

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

dropind = basisobj.dropind;
