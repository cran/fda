function coefmat = getcoef(fdobj)
%  GETCOEF   Extracts the coefficient array from FDOBJ.
%    If FDOBJ is a basis object, it assigns the identity matrix as the
%    coefficient array.

%  last modified 20 July 2006

  if isa_fd(fdobj) || isa_basis(fdobj)
    fnames = fieldnames(fdobj);
    switch fnames{1}
      case 'coef'
        coefmat = fdobj.coef;
      case 'type'
        basisobj = fdobj.basisobj;
        nbasis = getnbasis(basisobj);
        coefmat = eye(nbasis);
      otherwise
        error('Structure is neither of fd type or of basis type');
    end
  else
    error(['Argument is neither a functional data object', ...
           ' nor a functional basis object.']);
  end

