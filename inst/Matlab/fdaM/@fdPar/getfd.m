function fdobj = getfd(fdParobj)
%GETFD   Extracts the functional data object from 
%   a functional parameter object.

%  last modified 9 September 2003

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object.');
end

fdobj = fdParobj.fd;


