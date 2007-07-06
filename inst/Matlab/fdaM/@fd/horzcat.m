function horzcatfd = horzcat(xfd, yfd)
%   HORZCAT  concatenates two fd objects.
%   It is assumed that the objects have the same basisobj objects,
%   and that all the coef arrays have the same first dimension,
%   and, if there are more than one function, the same third dimension.
%   Returns CATFD, the concatenated functional data object

%   last modified 20 July 2006

  if ~(isa_fd(xfd) && isa_fd(yfd))
    error('Both arguments must be of functional data objects');
  end

  coefs    = getcoef(xfd);
  basisobj = getbasis(xfd);
  dimcoef  = size(coefs);
  ndimcoef = length(dimcoef);
  if yfd.basisobj ~= basisobj
    error('Objects must all have the same basis');
  end
  if length(size(getcoef(yfd))) ~= ndimcoef
    error('Objects must all have the same number of multiple functions');
  end
  if ndimcoef == 2
    coefs = [coefs, getcoef(yfd)];
  else
    coefs = [coefs, permute(getcoef(yfd), [1, 3, 2])];
    coefs = reshape(coefs,[dimcoef(1), dimcoef(3), length(coefs) ...
               /(dimcoef(1) * dimcoef(3))]);
    coefs = permute(coefs, [1, 3, 2]);
  end
  horzcatfd = fd(coefs, basisobj);

