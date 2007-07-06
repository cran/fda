function basisobj = create_power_basis(rangeval, nbasis, exponents, ...
                                       dropind)
%  CREATE_POWER_BASIS  Creates a power basis:, x.^exponents(1), x.^exponents(2), ...
%  Argument:
%  RANGEVAL ... an array of length 2 containing the lower and upper
%               boundaries for the rangeval of argument values.  If a
%               single value is input, it must be positive and the lower
%               limit of the range is set to 0.
%               For the power basis, the lower limit must not be negative.
%  NBASIS    ... number of basis functions
%  EXPONENTS ... an array of NBASIS exponents
%                by default this is 0:(NBASIS-1)
%  DROPIND  ... a set of indices in 1:NBASIS of basis functions to drop
%                when basis objects are arguments.  Default is [];
%  Return:
%  BASISOBJ ... a functional data basis object of type 'power'

%  Last modified 3 January 2008

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

if rangeval(1) < 0
    error('Lower limit in RANGEVAL is negative.');
end

%  set default argument values

if nargin < 2, nbasis = 2;               end
if nargin < 3, exponents = 0:(nbasis-1); end

%  check whether exponents are negative, 
%  and if so whether the range includes nonpostive values

if any(exponents < 0) && rangeval(1) <= 0
    error('An exponent is negative and range contains 0.');
end

% check if there are duplicate exponents

if min(diff(sort(exponents))) <= 0
    error('There are duplicate exponents.');
end

type        = 'power';
params      = sort(exponents);
dropind     = [];
quadvals    = [];
values      = {};
basisvalues = {};

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values, basisvalues);
