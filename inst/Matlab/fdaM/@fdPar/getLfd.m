function Lfdobj = getLfd(fdParobj)
%GETFD   Extracts the linear differential operator object from 
%   a functional parameter object.

%  last modified 6 May 2003

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end

Lfdobj = fdParobj.Lfd;


