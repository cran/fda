function basisobj = create_default_basis(argvals, nresol, nderiv, periodic)
%  CREATE_DEfAULT_BASIS takes a vector or matrix argvals and creates a 
%   default basis to be used for data observed on these arguments.
%
%  ARGVALS  ... A vector or matrix of argument values. Missing values 
%               allowed.
%  NRESOL   ... A number that specifies the number of the finest features 
%               or events that are of interest that can occur within the
%               range of the argument values. By feature or event is
%               meant things like peaks, valleys, zero crossings,
%               plateaus, linear slopes, and so on.  NRESOL specifies
%               the amount of resolution required of the functional
%               data object.
%  NDERIV   ... A natural number, 0, 1, 2, ..., specifying the number
%               of derivatives that the functional data object must
%               have.  The default is 2.
%  PERIODIC ... If T, functions are treated as periodic, and=the
%               case of vector ARGVALS the
%               argument domain is extended below by one value to become
%                (ARGVALS(1) - (ARGVALS(N)-ARGVALS(1))/(N-1), ARGALS(N).
%               The default is F.
%
%  Returns an object of class BASISOBJ, a functional data basis object
%
%  See also BASIS, CREATE_BSPLINE_BASIS,  CREATE_CONSTANT_BASIS, 
%  CREATE_EXPONENTIAL_BASIS, CREATE_FD_BASIS, CREATE_FDVARIANCE_BASIS,
%  CREATE_FOURIER_BASIS, CREATE_MONOMIAL_BASIS, CREATE_POLYGONAL_BASIS, 
%  CREATE_POLYNOMIAL_BASIS, CREATE_POWER_BASIS, CREATE_POLYGONAL_BASIS, 
%  CREATE_FEM_BASIS, CREATE_PRODUCT_BASIS, CREATE_TP_BASIS

%  Last modified 3 October 2011

if nargin < 4, periodic = 0; end
if nargin < 3, nderiv = 2;   end

%  Check values used to set up basis

argvals = argvals(:);

n = length(argvals);

rangeval = [min(argvals), max(argvals)];

%  check NRESOL

nresol = round(nresol);
if(nresol < 1 || nresol > n)
    error('NRESOL is not between 1 and N.');
end

%  check NDERIV

nderiv = round(nderiv);
if (nderiv < 0 || nderiv > n - 1)
    error('NDERIV is not between 0 and N-1.');
end

%  Set up basis object.

if(periodic)
    rangeval(1) = rangeval(1) - diff(rangeval)/(n - 1);
    basisobj    = create_fourier_basis(rangeval, nresol);
else
    if(nresol == 1)
        basisobj = create_constant_basis(rangeval);
    else
        if(nderiv == 0 && nresol == n && length(argvalsd) ~= 2)
            basisobj = create_polygon_basis(argvals);
        else
            basisobj = create_bspline_basis(rangeval, nresol, nderiv + 2);
        end
    end
end


