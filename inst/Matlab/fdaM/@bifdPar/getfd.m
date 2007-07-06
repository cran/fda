function bifdobj = getfd(bifdParobj)
%GETFD   Extracts the bivariate functional data object from 
%   a bivariate functional parameter object.

%  last modified 14 October 2003

if ~isa_bifdPar(bifdParobj)
    error('Argument is not a bivariate functional parameter object.');
end

bifdobj = bifdParobj.bifd;


