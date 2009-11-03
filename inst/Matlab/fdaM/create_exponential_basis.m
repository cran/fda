function basisobj = create_exponential_basis(rangeval, nbasis, ...
                                             ratevec, dropind)
%  CREATE_EXPONENTIAL_BASIS  Creates a exponential basis:
%            exp[RATEVEC(1)*x], exp[RATEVEC(2)*x], ...
%  Argument:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%  NBASIS    ... number of basis functions
%  RATEVEC   ... an array of NBASIS rate values
%  DROPIND   ... a set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];
%  Return:
%  BASIS     ... a functional data basis object of type 'expon'
%
%  See also BASIS, CREATE_BSPLINE_BASIS, CREATE_FD_BASIS, 
%  CREATE_FDVARIANCE_BASIS, CREATE_FOURIER_BASIS, CREATE_MONOMIAL_BASIS, 
%  CREATE_POLYGONAL_BASIS,  CREATE_POLYNOMIAL_BASIS, CREATE_POWER_BASIS, 
%  CREATE_POLYGONAL_BASIS,  CREATE_FEM_BASIS, CREATE_PRODUCT_BASIS, 
%  CREATE_TP_BASIS

%  Last modified 3 October 2011

%  default RANGEVAL

if nargin < 1, rangeval = [0,1];  end

%  check RANGEVAL

if length(rangeval) == 1
    if rangeval <= 0
        error('RANGEVAL a single value that is not positive.');
    end
    rangeval = [0,rangeval];
end

if rangechk(rangeval) ~= 1
    error('RANGEVAL is not a legitimate range.');
end

%  Set up some default arguments

if nargin < 2, nbasis  = 1;            end
if nargin < 3, ratevec = 0:(nbasis-1); end
if nargin < 4, dropind = [];           end

% check if there are duplicate values in ratevec

if min(diff(sort(ratevec))) <= 0
    error('There are duplicate ratevec.');
end

%  check DROPIND

if length(dropind) > 0
    if length(dropind) >= nbasis
        error('Too many index values in DROPIND.');
    end
    dropind = sort(dropind);
    if length(dropind) > 1
        if any(diff(dropind)) == 0
            error('Multiple index values in DROPIND.');
        end
    end
    for i=1:length(dropind);
        if dropind(i) < 1 || dropind(i) > nbasis
            error('An index value is out of range.');
        end
    end
end

type        = 'expon';
params      = ratevec;
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);

