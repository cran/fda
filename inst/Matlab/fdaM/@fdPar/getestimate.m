function estimate = getestimate(fdParobj)
%GETESTIMATE   Extracts the estimate value from 
%   a functional parameter object.

%  last modified 6 May 2003

if ~isa_fdPar(fdParobj)
    error('Argument is not a functional parameter object');
end

estimate = fdParobj.estimate;


