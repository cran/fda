function basisobj = create_constant_basis(rangeval)
%  CREATE_CONSTANT_BASIS  Creates a constant basis
%  Argument:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%               If no argument is present the range is set to [0,1];
%  Return:
%  BASIS.FD  ... a functional data basis object of type 'constant'
%
%  See also BASIS, CREATE_BSPLINE_BASIS, CREATE_EXPONENTIAL_BASIS, 
%  CREATE_FD_BASIS, CREATE_FDVARIANCE_BASIS, CREATE_FOURIER_BASIS, 
%  CREATE_MONOMIAL_BASIS, CREATE_POLYGONAL_BASIS, CREATE_POLYNOMIAL_BASIS, 
%  CREATE_POWER_BASIS, CREATE_POLYGONAL_BASIS, CREATE_FEM_BASIS, 
%  CREATE_PRODUCT_BASIS, CREATE_TP_BASIS

%  Last modified 3 October 2011

%  default rangeval

if nargin < 1, rangeval = [0,1];  end

%  check RANGEVAL

if length(rangeval) == 1
    if rangeval <= 0
        error('RANGEVAL a single value that is not postive.');
    end
    rangeval = [0,rangeval];
end

if rangechk(rangeval) ~= 1
    error('RANGEVAL is not a legitimate range.');
end

if nargin < 1
    rangeval = [0,1];
end

type        = 'const';
params      = [];
nbasis      = 1;
dropind     = [];
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);

