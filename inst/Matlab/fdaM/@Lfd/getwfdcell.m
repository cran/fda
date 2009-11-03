function wfdcell = getwfdcell(Lfdobj)
%  GETWFDCELL   Extracts the weight function cell object from LFDOBJ.

%  last modified 18 November 2003

if ~isa_Lfd(Lfdobj)
    error('Argument is not a linear differential operator object');
end

wfdcell = Lfdobj.bwtcell;


