function Lfdcell = getLfd(bifdParobj)
%GETLFD   Extracts the linear differential operator object from 
%   a bivariate functional parameter object.

%  last modified 14 October 2003

if ~isa_bifdPar(bifdParobj)
    error('Argument is not a bivariate functional parameter object');
end

Lfdcell{1} = bifdParobj.Lfds;
Lfdcell{2} = bifdParobj.Lfdt;


