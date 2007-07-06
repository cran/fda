function values = getvalues(basisobj, nderiv)
%  Return values of the derivative of order NDERIV of basis 
%    functions at quadrature points multiplied by the square 
%    root of the quadrature weights.
%  Arguments:
%  BASISOBJ ... a basis object
%  NDERIV   ... the order of the derivative to be retrieved.
%               This can be from 0 to the highest order of 
%               derivative that is stored in BASISOBJ.VALUES.   
%               If this is not present, then the entire
%               cell array containing all of the available
%               derivative matrices is returned.

%  Last modified 30 March 2006

%  check BASISOBJ

if nargin < 1
    error('No basis object supplied.');
end

%  check NDERIV

if nargin == 2
    if nderiv < 0
        error('Order of derivative is negative.');
    end
end

%  check if values are available

if isempty(basisobj.values)
    error('No basis function values available.');
end

%  
if nargin < 2
    values = basisobj.values{1};
else
    nderivcal = length(basisobj.values);
    if nderivcal >= nderiv+1
        values = basisobj.values{nderiv+1};
    else
        error('Derivative not contained in VALUES.');
    end
end
