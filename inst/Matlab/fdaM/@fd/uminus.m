function minusfd = uminus(fdobj)
% Unary minus or negative of functional data object.

%  last modified 24 April 2003

if ~isa_fd(fdobj)
    error('FD is not a functional data object.');
end
coef     = getcoef(fdobj);
basisobj = getbasis(fdobj);
fdnames  = getnames(fdobj);
minusfd  = fd(-coef, basisobj, fdnames);

