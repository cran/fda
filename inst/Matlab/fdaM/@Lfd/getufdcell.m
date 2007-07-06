function ufdcell = getufdcell(Lfdobj)
%  GETUFD   Extracts the forcing function from LFDOBJ.

%  last modified 5 November 2003

if ~isa_Lfd(Lfdobj)
    error('Argument is not a linear differential operator object');
end

ufdcell = Lfdobj.ufdcell;


