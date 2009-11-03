function basisobj = create_power_basis(rangeval, nbasis, exponents)
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
%  Return:
%  BASIS  ... a functional data basis object of type 'power'

%  last modified 20 July 2006

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

if nargin < 2, nbasis = 2;    end
if nargin < 3, exponents = 0:(nbasis-1); end

%  check whether exponents are negative, and if so whether any x are zero

if any(exponents < 0) && any(x == 0)
    error('An exponent is negative and an argument is equal to zero.')
end

% check if there are duplicate exponents

if min(diff(sort(exponents))) <= 0
    error('There are duplicate exponents.');
end

type   = 'power';
params = exponents;

dropind   = [];
quadvals  = [];
values{1} = [];

basisobj = basis(type, rangeval, nbasis, params, ...
                 dropind, quadvals, values);
