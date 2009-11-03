function nderiv = getnderiv(Lfdobj)
%  GETNDERIV   Extracts the order of the operator from LFDOBJ.

%  last modified 3 January 2008

if ~isa_Lfd(Lfdobj)
    error('Argument is not a linear differential operator object');
end

nderiv = Lfdobj.nderiv;


