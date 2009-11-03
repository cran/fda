function basisvalues = getbasisvalues(basisobj, argvals, nderiv)
%  Return values of the derivative of order NDERIV of basis 
%    functions at quadrature points multiplied by the square 
%    root of the quadrature weights.
%  If only basisobj is supplied, the cell array containing
%  basis values is returned.
%  Arguments:
%  BASISOBJ ... a basis object
%  ARGVALS  ... a vector of argument values
%  NDERIV   ... the order of the derivative to be retrieved.
%               This can be from 0 to the highest order of 
%               derivative that is stored in BASISOBJ.VALUES.   
%               If this is not present, then the entire
%               cell array containing all of the available
%               derivative matrices is returned.

%  Last modified 10 March  2011 

%  check BASISOBJ

if nargin < 1
    error('No basis object supplied.');
end

%  check NDERIV

if nargin < 3
    nderiv = 0;
end

if nderiv < 0
    error('Order of derivative is negative.');
end

%  check if values are available

if isempty(basisobj.basisvalues)
    basisvalues = [];
    return;
end

%  if no arguments are supplied, return the cell array itself

if nargin < 2
    basisvalues = basisobj.basisvalues;
    return;
end

%  check that requested derivative is available

sizevec = size(basisobj.basisvalues);
if nderiv + 2 > sizevec(2)
    basisvalues = [];
    return;
end

%  search for argvals match

N = length(argvals);
for i=1:sizevec(1)
    evalarg = basisobj.basisvalues{i,1};
    if ~isempty(basisobj.basisvalues{i,nderiv+2})
        if N == length(evalarg)
            if all(argvals(:) == evalarg(:))
                basisvalues = basisobj.basisvalues{i,nderiv+2};
                return;
            end
        end
    end
end

basisvalues = [];