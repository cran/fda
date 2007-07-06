function [lambdas, lambdat] = getlambda(bifdParobj)
%GETLAMBDA   Extracts the lambda value from 
%   a bivariate functional parameter object.

%  last modified 14 October 2003

if ~isa_bifdPar(bifdParobj)
    error('Argument is not a bivariate functional parameter object');
end

lambdas = bifdParobj.lambdas;
lambdat = bifdParobj.lambdat;


