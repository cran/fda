function penmat = getpenmat(fdParobj)
%GETPENMAT   Extracts the penmat value from 
%   a functional parameter object.

%  last modified 17 May 2003

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end
penmat = fdParobj.penmat;


