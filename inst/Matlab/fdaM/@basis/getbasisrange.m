function rangeval = getbasisrange(basisobj)
%  GETBASISRANGE   Extracts the range from basis object BASISOBJ.

%  last modified 30 June 1998

  if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
  end

  rangeval = basisobj.rangeval;
