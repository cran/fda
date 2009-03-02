function lambda = getlambda(fdParobj)
%GETLAMBDA   Extracts the lambda value from 
%   a functional parameter object.

%  last modified 6 May 2003

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end

lambda = fdParobj.lambda;


