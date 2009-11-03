function fdParobj = putfd(fdParobj, fdobj)
%  Replace the fd parameter and reset penmat parameter to []

%  Last modified 17 May 2004

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end

fdParobj.fd     = fdobj;
fdParobj.penmat = [];