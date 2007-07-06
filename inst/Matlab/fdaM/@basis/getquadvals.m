function quadvals = getquadvals(basisobj)
%  GETQUADVALS   Extracts the quadrature points and weights
%     from basis object BASISOBJ.

%  last modified 30 March 2006

if ~isa_basis(basisobj)
    error('Argument is not a functional basis object.');
end

%  check if values are available

if isempty(basisobj.quadvals)
    error('No quadrature points and weights available.');
end

quadvals = basisobj.quadvals;
