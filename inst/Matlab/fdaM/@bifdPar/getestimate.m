function estimate = getestimate(bifdParobj)
%GETESTIMATE   Extracts the estimate value from 
%   a bivariate functional parameter object.

%  last modified 14 October 2003

if ~isa_bifdPar(bifdParobj)
    error('Argument is not a functional parameter object');
end

estimate = bifdParobj.estimate;


