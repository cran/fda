function nbasis = getnbasis(basisobj)
%  GETNBASIS   Extracts the type of basis from basis object BASISOBJ.

%  last modified 3 January 2008

  if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
  end

  nbasis = basisobj.nbasis;
