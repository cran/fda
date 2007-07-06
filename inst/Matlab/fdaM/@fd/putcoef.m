function newfdobj = putcoef(fdobj, coef)
%  PUTCOEF   Replaces the coefficient array from FDOBJ.

%  last modified 1 July 1998

  if isa_fd(fdobj)
    newfdobj.coef     = coef;
    newfdobj.basisobj = fdobj.basisobj;
    newfdobj.fdnames  = fdobj.fdnames;
    newfdobj = class(newfdobj, 'fd');
  else
    error('Argument FDOBJ is not a functional data object.');
  end

